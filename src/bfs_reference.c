/* Copyright (c) 2011-2017 Graph500 Steering Committee
   All rights reserved.
   Developed by:		Anton Korzh anton@korzh.us
				Graph500 Steering Committee
				http://www.graph500.org
   New code under University of Illinois/NCSA Open Source License
   see license.txt or https://opensource.org/licenses/NCSA
*/

// Graph500: Kernel 2: BFS
// Simple level-synchronized BFS with visits as Active Messages

#include "common.h"
#include "aml.h"
#include "statistics.h"
#include "csr_reference.h"
#include "bitmap_reference.h"
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <stdint.h>

#ifdef DEBUGSTATS
extern int64_t nbytes_sent,nbytes_rcvd;
#endif
// two arrays holding visited VERTEX_LOCALs for current and next level
// we swap pointers each time
int *q1,*q2;
int qc,q2c; //pointer to first free element

//VISITED bitmap parameters
unsigned long *visited;
int64_t visited_size;

//global variables of CSR graph to be used inside of AM-handlers
int64_t *column;
int64_t *pred_glob;
unsigned int * rowstarts;

oned_csr_graph g;

typedef struct visitmsg {
	//both vertexes are VERTEX_LOCAL components as we know src and dest PEs to reconstruct VERTEX_GLOBAL
	int vloc;
	int vfrom;
} visitmsg;

//AM-handler for check&visit
void visithndl(int from,void* data,int sz) {
	visitmsg *m = data;
	if (!TEST_VISITEDLOC(m->vloc)) {
		SET_VISITEDLOC(m->vloc);
		q2[q2c++] = m->vloc;
		pred_glob[m->vloc] = VERTEX_TO_GLOBAL(from,m->vfrom);
	}
}

void send_visit(int64_t glob, int from) {
	visitmsg m = {VERTEX_LOCAL(glob),from};
	aml_send(&m,1,sizeof(visitmsg),VERTEX_OWNER(glob));
}

double make_graph_data_structure(const tuple_graph* const tg) {
	double data_struct_start = MPI_Wtime();
	
	int i,j,k;
	convert_graph_to_oned_csr(tg, &g);
	column=g.column;
	rowstarts=g.rowstarts;
	visited_size = (g.nlocalverts + ulong_bits - 1) / ulong_bits;
	aml_register_handler(visithndl,1);
	q1 = xmalloc(g.nlocalverts*sizeof(int)); //100% of vertexes
	q2 = xmalloc(g.nlocalverts*sizeof(int));
	for(i=0;i<g.nlocalverts;i++) q1[i]=0,q2[i]=0; //touch memory
	visited = xmalloc(visited_size*sizeof(unsigned long));

	double data_struct_end = MPI_Wtime();
	
	return data_struct_end - data_struct_start;
}

void run_bfs(int64_t root, int64_t* pred) {
	int64_t nvisited;
	long sum;
	unsigned int i,j,k,lvl=1;
	pred_glob=pred;
	aml_register_handler(visithndl,1);

	CLEAN_VISITED();

	qc=0; sum=1; q2c=0;

	nvisited=1;
	if(VERTEX_OWNER(root) == rank) {
		pred[VERTEX_LOCAL(root)]=root;
		SET_VISITED(root);
		q1[0]=VERTEX_LOCAL(root);
		qc=1;
	} 

	// While there are vertices in current level
	while(sum) {
#ifdef DEBUGSTATS
		double t0=aml_time();
		nbytes_sent=0; nbytes_rcvd=0;
#endif
		//for all vertices in current level send visit AMs to all neighbours
		for(i=0;i<qc;i++)
			for(j=rowstarts[q1[i]];j<rowstarts[q1[i]+1];j++)
				send_visit(COLUMN(j),q1[i]);
				
		aml_barrier();

		qc=q2c;int *tmp=q1;q1=q2;q2=tmp;
		sum=qc;
		aml_long_allsum(&sum);

		nvisited+=sum;

		q2c=0;
#ifdef DEBUGSTATS
		aml_long_allsum(&nbytes_sent);
		t0-=aml_time();
		if(!my_pe()) printf (" --lvl%d : %lld(%lld,%3.2f) visited in %5.2fs, network aggr %5.2fGb/s\n",lvl++,sum,nvisited,((double)nvisited/(double)g.notisolated)*100.0,-t0,-(double)nbytes_sent*8.0/(1.e9*t0));
#endif
	}
	aml_barrier();

}

void run_bfs_for_roots(tuple_graph *tg, int num_bfs_roots, int64_t *bfs_roots,
					stat_t *statis, uint64_t nlocalverts, int64_t *pred, float *shortest)
{
	int bfs_root_idx,i;

	statis -> validation_passed = 1;
	statis -> bfs_times = (double*)xmalloc(num_bfs_roots * sizeof(double));
	statis -> validate_times = (double*)xmalloc(num_bfs_roots * sizeof(double));
	statis -> edge_counts = (double*)xmalloc(num_bfs_roots * sizeof(double));

	clean_pred(&pred[0]); //user-provided function from bfs_implementation.c
	run_bfs(bfs_roots[0], &pred[0]); //warm-up
#ifdef ENERGYLOOP_BFS
            int eloop;
            if(!my_pe()) printf("starting energy loop BFS\n");
            for(eloop=0;eloop<1000000;eloop++)
                    for (bfs_root_idx = 0; bfs_root_idx < num_bfs_roots; ++bfs_root_idx) {
            clean_pred(&pred[0]);
            run_bfs(bfs_roots[bfs_root_idx], &pred[0]);
            }
            if(!my_pe()) printf("finished energy loop BFS\n");
#endif
	if (!getenv("SKIP_VALIDATION")) {
		int64_t nedges=0;
		validate_result(1, tg, nlocalverts, bfs_roots[0], pred,shortest,NULL);
	}

	for (bfs_root_idx = 0; bfs_root_idx < num_bfs_roots; ++bfs_root_idx) {
		int64_t root = bfs_roots[bfs_root_idx];

		//if (rank == 0) fprintf(stderr, "Running BFS %d\n", bfs_root_idx);

		clean_pred(&pred[0]); //user-provided function from bfs_implementation.c
		/* Do the actual BFS. */
		double bfs_start = MPI_Wtime();
		run_bfs(root, &pred[0]);
		double bfs_stop = MPI_Wtime();
		statis -> bfs_times[bfs_root_idx] = bfs_stop - bfs_start;
		//if (rank == 0) fprintf(stderr, "Time for BFS %d is %f\n", bfs_root_idx, bfs_times[bfs_root_idx]);
		int64_t edge_visit_count=0;
		get_edge_count_for_teps(&edge_visit_count);
		statis -> edge_counts[bfs_root_idx] = (double)edge_visit_count;
		//if (rank == 0) fprintf(stderr, "TEPS for BFS %d is %g\n", bfs_root_idx, edge_visit_count / bfs_times[bfs_root_idx]);

		/* Validate result. */
		if (!getenv("SKIP_VALIDATION")) {
			//if (rank == 0) fprintf(stderr, "Validating BFS %d\n", bfs_root_idx);

			double validate_start = MPI_Wtime();
			int validation_passed_one = validate_result(1, tg, nlocalverts, root, pred,shortest,&edge_visit_count);
			double validate_stop = MPI_Wtime();

			statis -> validate_times[bfs_root_idx] = validate_stop - validate_start;
			//if (rank == 0) fprintf(stderr, "Validate time for BFS %d is %f\n", bfs_root_idx, validate_times[bfs_root_idx]);

			if (!validation_passed_one) {
				statis -> validation_passed = 0;
				//if (rank == 0) fprintf(stderr, "Validation failed for this BFS root; skipping rest.\n");
				break;
			}
		} else
			statis -> validate_times[bfs_root_idx] = -1;
	}
}

//we need edge count to calculate teps. Validation will check if this count is correct
void get_edge_count_for_teps(int64_t* edge_visit_count) {
	long i,j;
	long edge_count=0;
	for(i=0;i<g.nlocalverts;i++)
		if(pred_glob[i]!=-1) {
			for(j=rowstarts[i];j<rowstarts[i+1];j++)
				if(COLUMN(j)<=VERTEX_TO_GLOBAL(my_pe(),i))
					edge_count++;

		}

	aml_long_allsum(&edge_count);
	*edge_visit_count=edge_count;
}

void clean_pred(int64_t* pred) {
	int i;
	for(i=0;i<g.nlocalverts;i++) pred[i]=-1;
}
void free_graph_data_structure(void) {
	int i; 
	free_oned_csr_graph(&g);
	free(q1); free(q2); free(visited);
}

size_t get_nlocalverts_for_pred(void) {
	return g.nlocalverts;
}
