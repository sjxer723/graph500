#include "generate_graph.h"

double GenerateGraph(int SCALE, uint64_t seed1, uint64_t seed2, int64_t nglobalverts, tuple_graph* tg)
{
	double make_graph_start = MPI_Wtime();
	if( !tg -> data_in_file || tg -> write_file )
	{
		/* Spread the two 64-bit numbers into five nonzero values in the correct
		 * range. */
		uint_fast32_t seed[5];
		make_mrg_seed(seed1, seed2, seed);

		/* As the graph is being generated, also keep a bitmap of vertices with
		 * incident edges.  We keep a grid of processes, each row of which has a
		 * separate copy of the bitmap (distributed among the processes in the
		 * row), and then do an allreduce at the end.  This scheme is used to avoid
		 * non-local communication and reading the file separately just to find BFS
		 * roots. */
		MPI_Offset nchunks_in_file = (tg -> nglobaledges + FILE_CHUNKSIZE - 1) / FILE_CHUNKSIZE;
		int64_t bitmap_size_in_bytes = int64_min(BITMAPSIZE, (nglobalverts + CHAR_BIT - 1) / CHAR_BIT);
		if (bitmap_size_in_bytes * size * CHAR_BIT < nglobalverts) {
			bitmap_size_in_bytes = (nglobalverts + size * CHAR_BIT - 1) / (size * CHAR_BIT);
		}
		int ranks_per_row = tg -> data_in_file ? ((nglobalverts + CHAR_BIT - 1) / CHAR_BIT + bitmap_size_in_bytes - 1) / bitmap_size_in_bytes : 1;
		int nrows = size / ranks_per_row;
		int my_row = -1, my_col = -1;
		MPI_Comm cart_comm;
		{
			int dims[2] = {size / ranks_per_row, ranks_per_row};
			int periods[2] = {0, 0};
			MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);
		}
		int in_generating_rectangle = 0;
		if (cart_comm != MPI_COMM_NULL) {
			in_generating_rectangle = 1;
			{
				int dims[2], periods[2], coords[2];
				MPI_Cart_get(cart_comm, 2, dims, periods, coords);
				my_row = coords[0];
				my_col = coords[1];
			}
			MPI_Comm this_col;
			MPI_Comm_split(cart_comm, my_col, my_row, &this_col);
			MPI_Comm_free(&cart_comm);
			/* Every rank in a given row creates the same vertices (for updating the
			 * bitmap); only one writes them to the file (or final memory buffer). */
			packed_edge* buf = (packed_edge*)xmalloc(FILE_CHUNKSIZE * sizeof(packed_edge));
#ifdef SSSP
			float* wbuf = (float*)xmalloc(FILE_CHUNKSIZE*sizeof(float));
#endif
			MPI_Offset block_limit = (nchunks_in_file + nrows - 1) / nrows;
			/* fprintf(stderr, "%d: nchunks_in_file = %" PRId64 ", block_limit = %" PRId64 " in grid of %d rows, %d cols\n", rank, (int64_t)nchunks_in_file, (int64_t)block_limit, nrows, ranks_per_row); */
			if (tg -> data_in_file) {
				tg -> edgememory_size = 0;
				tg -> edgememory = NULL;
			} else {
				int my_pos = my_row + my_col * nrows;
				int last_pos = (tg -> nglobaledges % ((int64_t)FILE_CHUNKSIZE * nrows * ranks_per_row) != 0) ?
					(tg -> nglobaledges / FILE_CHUNKSIZE) % (nrows * ranks_per_row) :
					-1;
				int64_t edges_left = tg -> nglobaledges % FILE_CHUNKSIZE;
				int64_t nedges = FILE_CHUNKSIZE * (tg -> nglobaledges / ((int64_t)FILE_CHUNKSIZE * nrows * ranks_per_row)) +
					FILE_CHUNKSIZE * (my_pos < (tg -> nglobaledges / FILE_CHUNKSIZE) % (nrows * ranks_per_row)) +
					(my_pos == last_pos ? edges_left : 0);
				/* fprintf(stderr, "%d: nedges = %" PRId64 " of %" PRId64 "\n", rank, (int64_t)nedges, (int64_t)tg -> nglobaledges); */
				tg -> edgememory_size = nedges;
				tg -> edgememory = (packed_edge*)xmalloc(nedges * sizeof(packed_edge));
#ifdef SSSP
				tg -> weightmemory = (float*)xmalloc(nedges*sizeof(float));
#endif
			}
			MPI_Offset block_idx;
			for (block_idx = 0; block_idx < block_limit; ++block_idx) {
				/* fprintf(stderr, "%d: On block %d of %d\n", rank, (int)block_idx, (int)block_limit); */
				MPI_Offset start_edge_index = int64_min(FILE_CHUNKSIZE * (block_idx * nrows + my_row), tg -> nglobaledges);
				MPI_Offset edge_count = int64_min(tg -> nglobaledges - start_edge_index, FILE_CHUNKSIZE);
				packed_edge* actual_buf = (!tg -> data_in_file && block_idx % ranks_per_row == my_col) ?
					tg -> edgememory + FILE_CHUNKSIZE * (block_idx / ranks_per_row) :
					buf;
#ifdef SSSP
				float* actual_wbuf = (!tg -> data_in_file && block_idx % ranks_per_row == my_col) ?
					tg -> weightmemory + FILE_CHUNKSIZE * (block_idx / ranks_per_row) :
					wbuf;
#endif
				/* fprintf(stderr, "%d: My range is [%" PRId64 ", %" PRId64 ") %swriting into index %" PRId64 "\n", rank, (int64_t)start_edge_index, (int64_t)(start_edge_index + edge_count), (my_col == (block_idx % ranks_per_row)) ? "" : "not ", (int64_t)(FILE_CHUNKSIZE * (block_idx / ranks_per_row))); */
				if (!tg -> data_in_file && block_idx % ranks_per_row == my_col) {
					assert (FILE_CHUNKSIZE * (block_idx / ranks_per_row) + edge_count <= tg -> edgememory_size);
				}
				if (tg -> write_file) {
					generate_kronecker_range(seed, SCALE, start_edge_index, start_edge_index + edge_count, actual_buf
#ifdef SSSP
							,actual_wbuf
#endif
							);
					if (tg -> data_in_file && my_col == (block_idx % ranks_per_row)) { /* Try to spread writes among ranks */
						MPI_File_write_at(tg -> edgefile, start_edge_index, actual_buf, edge_count, packed_edge_mpi_type, MPI_STATUS_IGNORE);
#ifdef SSSP
						MPI_File_write_at(tg -> weightfile, start_edge_index, actual_wbuf, edge_count, MPI_FLOAT, MPI_STATUS_IGNORE);
#endif
					}
				} else {
					/* All read rather than syncing up for a row broadcast. */
					MPI_File_read_at(tg -> edgefile, start_edge_index, actual_buf, edge_count, packed_edge_mpi_type, MPI_STATUS_IGNORE);
#ifdef SSSP
					MPI_File_read_at(tg -> weightfile, start_edge_index, actual_wbuf, edge_count, MPI_FLOAT, MPI_STATUS_IGNORE);
#endif
				}
			}
			free(buf);
			MPI_Comm_free(&this_col);
		} else {
			tg -> edgememory = NULL;
			tg -> edgememory_size = 0;
#ifdef SSSP
			tg -> weightmemory = NULL;
#endif
		}
		MPI_Allreduce(&tg -> edgememory_size, &tg -> max_edgememory_size, 1, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);
		if (tg -> data_in_file && tg -> write_file) {
			MPI_File_sync(tg -> edgefile);
#ifdef SSSP
			MPI_File_sync(tg -> weightfile);
#endif
		}
	}

	double make_graph_stop = MPI_Wtime();
	double make_graph_time = make_graph_stop - make_graph_start;
	if (rank == 0) { /* Not an official part of the results */
		fprintf(stderr, "graph_generation:               %f s\n", make_graph_time);
	}

	return make_graph_time;
}
