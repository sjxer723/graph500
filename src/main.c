/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */
/*           Anton Korzh                                                   */


/* These need to be before any possible inclusions of stdint.h or inttypes.h.
 * */
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include "../generator/make_graph.h"
#include "../generator/utils.h"
#include "aml.h"
#include "generate_graph.h"
#include "csr_reference.h"
#include "common.h"
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>

//int isisolated(int64_t v);
static int compare_doubles(const void* a, const void* b) {
	double aa = *(const double*)a;
	double bb = *(const double*)b;
	return (aa < bb) ? -1 : (aa == bb) ? 0 : 1;
}

enum {s_minimum, s_firstquartile, s_median, s_thirdquartile, s_maximum, s_mean, s_std, s_LAST};
void get_statistics(const double x[], int n, volatile double r[s_LAST]) {
	double temp;
	int i;
	/* Compute mean. */
	temp = 0.0;
	for (i = 0; i < n; ++i) temp += x[i];
	temp /= n;
	r[s_mean] = temp;
	double mean = temp;
	/* Compute std. dev. */
	temp = 0;
	for (i = 0; i < n; ++i) temp += (x[i] - mean) * (x[i] - mean);
	temp /= n - 1;
	r[s_std] = sqrt(temp);
	/* Sort x. */
	double* xx = (double*)xmalloc(n * sizeof(double));
	memcpy(xx, x, n * sizeof(double));
	qsort(xx, n, sizeof(double), compare_doubles);
	/* Get order statistics. */
	r[s_minimum] = xx[0];
	r[s_firstquartile] = (xx[(n - 1) / 4] + xx[n / 4]) * .5;
	r[s_median] = (xx[(n - 1) / 2] + xx[n / 2]) * .5;
	r[s_thirdquartile] = (xx[n - 1 - (n - 1) / 4] + xx[n - 1 - n / 4]) * .5;
	r[s_maximum] = xx[n - 1];
	/* Clean up. */
	free(xx);
}

int main(int argc, char** argv) {
	int SCALE = 16;
	int edgefactor = 16; /* nedges / nvertices, i.e., 2*avg. degree */
	int num_bfs_roots = 64;
	int64_t* bfs_roots = NULL;
	uint64_t seed1 = 2, seed2 = 3;
	double make_graph_time = 0;
	
	/* Parse arguments. */
	aml_init(&argc,&argv); //includes MPI_Init inside
	setup_globals();
	if (argc >= 2) SCALE = atoi(argv[1]);
	if (argc >= 3) edgefactor = atoi(argv[2]);
	if (argc <= 1 || argc >= 4 || SCALE == 0 || edgefactor == 0) {
		if (rank == 0) {
			fprintf(stderr, "Usage: %s SCALE edgefactor\n  SCALE = log_2(# vertices) [integer, required]\n  edgefactor = (# edges) / (# vertices) = .5 * (average vertex degree) [integer, defaults to 16]\n(Random number seed and Kronecker initiator are in main.c)\n", argv[0]);
		}
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	/* Initial the tuple graph tg */
	tuple_graph tg;
	InitTupleGraph(SCALE, seed1, seed2, edgefactor, &tg);
	
	/* Make graph according to SCALE and nglobalverts */
	int64_t nglobalverts = (int64_t)(1) << SCALE;
	make_graph_time = GenerateGraph(SCALE, seed1, seed2, nglobalverts, &tg);

	/* Make user's graph data structure. */
	double data_struct_start = MPI_Wtime();
	make_graph_data_structure(&tg);
	double data_struct_stop = MPI_Wtime();
	double data_struct_time = data_struct_stop - data_struct_start;
	if (rank == 0) { /* Not an official part of the results */
		fprintf(stderr, "construction_time:              %f s\n", data_struct_time);
	}

	/* Make the raw graph edges. */
	/* Get roots for BFS runs, plus maximum vertex with non-zero degree (used by
	 * validator). */
	//generate non-isolated roots
	num_bfs_roots = GenerateNonIsolatedRoots(seed1, seed2, nglobalverts, num_bfs_roots, &bfs_roots); 
	
	/* Number of edges visited in each BFS; a double so get_statistics can be
	 * used directly. */
	double* edge_counts = (double*)xmalloc(num_bfs_roots * sizeof(double));

	/* Run BFS. */
	int validation_passed = 1;
	double* bfs_times = (double*)xmalloc(num_bfs_roots * sizeof(double));
	double* validate_times = (double*)xmalloc(num_bfs_roots * sizeof(double));
	uint64_t nlocalverts = get_nlocalverts_for_pred();
	int64_t* pred = (int64_t*)xMPI_Alloc_mem(nlocalverts * sizeof(int64_t));
	float* shortest = (float*)xMPI_Alloc_mem(nlocalverts * sizeof(float));


	int bfs_root_idx,i;
	if (!getenv("SKIP_BFS")) {
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
			validate_result(1,&tg, nlocalverts, bfs_roots[0], pred,shortest,NULL);
		}

		for (bfs_root_idx = 0; bfs_root_idx < num_bfs_roots; ++bfs_root_idx) {
			int64_t root = bfs_roots[bfs_root_idx];

			if (rank == 0) fprintf(stderr, "Running BFS %d\n", bfs_root_idx);

			clean_pred(&pred[0]); //user-provided function from bfs_implementation.c
			/* Do the actual BFS. */
			double bfs_start = MPI_Wtime();
			run_bfs(root, &pred[0]);
			double bfs_stop = MPI_Wtime();
			bfs_times[bfs_root_idx] = bfs_stop - bfs_start;
			if (rank == 0) fprintf(stderr, "Time for BFS %d is %f\n", bfs_root_idx, bfs_times[bfs_root_idx]);
			int64_t edge_visit_count=0;
			get_edge_count_for_teps(&edge_visit_count);
			edge_counts[bfs_root_idx] = (double)edge_visit_count;
			if (rank == 0) fprintf(stderr, "TEPS for BFS %d is %g\n", bfs_root_idx, edge_visit_count / bfs_times[bfs_root_idx]);

			/* Validate result. */
			if (!getenv("SKIP_VALIDATION")) {
				if (rank == 0) fprintf(stderr, "Validating BFS %d\n", bfs_root_idx);

				double validate_start = MPI_Wtime();
				int validation_passed_one = validate_result(1,&tg, nlocalverts, root, pred,shortest,&edge_visit_count);
				double validate_stop = MPI_Wtime();

				validate_times[bfs_root_idx] = validate_stop - validate_start;
				if (rank == 0) fprintf(stderr, "Validate time for BFS %d is %f\n", bfs_root_idx, validate_times[bfs_root_idx]);

				if (!validation_passed_one) {
					validation_passed = 0;
					if (rank == 0) fprintf(stderr, "Validation failed for this BFS root; skipping rest.\n");
					break;
				}
			} else
				validate_times[bfs_root_idx] = -1;
		}

	}
#ifdef SSSP
	double* sssp_times = (double*)xmalloc(num_bfs_roots * sizeof(double));
	double* validate_times2 = (double*)xmalloc(num_bfs_roots * sizeof(double));

	clean_shortest(shortest);
	clean_pred(pred);
	run_sssp(bfs_roots[0], &pred[0],shortest); //warm-up
#ifdef ENERGYLOOP_SSSP
		int eloop;
		if(!my_pe()) printf("starting energy loop SSSP\n");
		for(eloop=0;eloop<1000000;eloop++)
			for (bfs_root_idx = 0; bfs_root_idx < num_bfs_roots; ++bfs_root_idx) {
				clean_shortest(shortest);
				clean_pred(&pred[0]);
				run_sssp(bfs_roots[bfs_root_idx], &pred[0],shortest);
		}
		if(!my_pe()) printf("finished energy loop SSSP\n");
#endif
	for (bfs_root_idx = 0; bfs_root_idx < num_bfs_roots; ++bfs_root_idx) {
		int64_t root = bfs_roots[bfs_root_idx];

		if (rank == 0) fprintf(stderr, "Running SSSP %d\n", bfs_root_idx);

		clean_pred(&pred[0]);
		clean_shortest(shortest);

		/* Do the actual SSSP. */
		double sssp_start = MPI_Wtime();
		run_sssp(root, &pred[0],shortest);
		double sssp_stop = MPI_Wtime();
		sssp_times[bfs_root_idx] = sssp_stop - sssp_start;
		int64_t edge_visit_count=0;
		get_edge_count_for_teps(&edge_visit_count);
		edge_counts[bfs_root_idx] = (double)edge_visit_count;
		if (rank == 0) fprintf(stderr, "Time for SSSP %d is %f\n", bfs_root_idx, sssp_times[bfs_root_idx]);
		if (rank == 0) fprintf(stderr, "TEPS for SSSP %d is %g\n", bfs_root_idx, edge_counts[bfs_root_idx] / sssp_times[bfs_root_idx]);

		/* Validate result. */
		if (!getenv("SKIP_VALIDATION")) {
			if (rank == 0) fprintf(stderr, "Validating SSSP %d\n", bfs_root_idx);

			double validate_start = MPI_Wtime();
			int validation_passed_one = validate_result(0,&tg, nlocalverts, root, pred, shortest,&edge_visit_count);
			double validate_stop = MPI_Wtime();

			validate_times2[bfs_root_idx] = validate_stop - validate_start;
			if (rank == 0) fprintf(stderr, "Validate time for SSSP %d is %f\n", bfs_root_idx, validate_times2[bfs_root_idx]);

			if (!validation_passed_one) {
				validation_passed = 0;
				if (rank == 0) fprintf(stderr, "Validation failed for this SSSP root; skipping rest.\n");
				break;
			}
		} else {
			validate_times2[bfs_root_idx] = -1;
		}
	}

#endif
	MPI_Free_mem(pred);
#ifdef SSSP
	MPI_Free_mem(shortest);
#endif
	free(bfs_roots);
	free_graph_data_structure();

	if (tg.data_in_file) {
		MPI_File_close(&tg.edgefile);
#ifdef SSSP
		MPI_File_close(&tg.weightfile);
#endif    
	} else {
		free(tg.edgememory); tg.edgememory = NULL;
#ifdef SSSP
		free(tg.weightmemory); tg.weightmemory = NULL;
#endif
	}

	/* Print results. */
	if (rank == 0) {
		if (!validation_passed) {
			fprintf(stdout, "No results printed for invalid run.\n");
		} else {
			int i;
			//for (i = 0; i < num_bfs_roots; ++i) printf(" %g \n",edge_counts[i]);
			fprintf(stdout, "SCALE:                          %d\n", SCALE);
			fprintf(stdout, "edgefactor:                     %d\n", edgefactor);
			fprintf(stdout, "NBFS:                           %d\n", num_bfs_roots);
			fprintf(stdout, "graph_generation:               %g\n", make_graph_time);
			fprintf(stdout, "num_mpi_processes:              %d\n", size);
			fprintf(stdout, "construction_time:              %g\n", data_struct_time);
			volatile double stats[s_LAST];
			get_statistics(bfs_times, num_bfs_roots, stats);
			fprintf(stdout, "bfs  min_time:                  %g\n", stats[s_minimum]);
			fprintf(stdout, "bfs  firstquartile_time:        %g\n", stats[s_firstquartile]);
			fprintf(stdout, "bfs  median_time:               %g\n", stats[s_median]);
			fprintf(stdout, "bfs  thirdquartile_time:        %g\n", stats[s_thirdquartile]);
			fprintf(stdout, "bfs  max_time:                  %g\n", stats[s_maximum]);
			fprintf(stdout, "bfs  mean_time:                 %g\n", stats[s_mean]);
			fprintf(stdout, "bfs  stddev_time:               %g\n", stats[s_std]);
#ifdef SSSP
			get_statistics(sssp_times, num_bfs_roots, stats);
			fprintf(stdout, "sssp min_time:                  %g\n", stats[s_minimum]);
			fprintf(stdout, "sssp firstquartile_time:        %g\n", stats[s_firstquartile]);
			fprintf(stdout, "sssp median_time:               %g\n", stats[s_median]);
			fprintf(stdout, "sssp thirdquartile_time:        %g\n", stats[s_thirdquartile]);
			fprintf(stdout, "sssp max_time:                  %g\n", stats[s_maximum]);
			fprintf(stdout, "sssp mean_time:                 %g\n", stats[s_mean]);
			fprintf(stdout, "sssp stddev_time:               %g\n", stats[s_std]);
#endif
			get_statistics(edge_counts, num_bfs_roots, stats);
			fprintf(stdout, "min_nedge:                      %.11g\n", stats[s_minimum]);
			fprintf(stdout, "firstquartile_nedge:            %.11g\n", stats[s_firstquartile]);
			fprintf(stdout, "median_nedge:                   %.11g\n", stats[s_median]);
			fprintf(stdout, "thirdquartile_nedge:            %.11g\n", stats[s_thirdquartile]);
			fprintf(stdout, "max_nedge:                      %.11g\n", stats[s_maximum]);
			fprintf(stdout, "mean_nedge:                     %.11g\n", stats[s_mean]);
			fprintf(stdout, "stddev_nedge:                   %.11g\n", stats[s_std]);
			double* secs_per_edge = (double*)xmalloc(num_bfs_roots * sizeof(double));
			for (i = 0; i < num_bfs_roots; ++i) secs_per_edge[i] = bfs_times[i] / edge_counts[i];
			get_statistics(secs_per_edge, num_bfs_roots, stats);
			fprintf(stdout, "bfs  min_TEPS:                  %g\n", 1. / stats[s_maximum]);
			fprintf(stdout, "bfs  firstquartile_TEPS:        %g\n", 1. / stats[s_thirdquartile]);
			fprintf(stdout, "bfs  median_TEPS:               %g\n", 1. / stats[s_median]);
			fprintf(stdout, "bfs  thirdquartile_TEPS:        %g\n", 1. / stats[s_firstquartile]);
			fprintf(stdout, "bfs  max_TEPS:                  %g\n", 1. / stats[s_minimum]);
			fprintf(stdout, "bfs  harmonic_mean_TEPS:     !  %g\n", 1. / stats[s_mean]);
			/* Formula from:
			 * Title: The Standard Errors of the Geometric and Harmonic Means and
			 *        Their Application to Index Numbers
			 * Author(s): Nilan Norris
			 * Source: The Annals of Mathematical Statistics, Vol. 11, No. 4 (Dec., 1940), pp. 445-448
			 * Publisher(s): Institute of Mathematical Statistics
			 * Stable URL: http://www.jstor.org/stable/2235723
			 * (same source as in specification). */
			fprintf(stdout, "bfs  harmonic_stddev_TEPS:      %g\n", stats[s_std] / (stats[s_mean] * stats[s_mean] * sqrt(num_bfs_roots - 1)));
#ifdef SSSP
			for (i = 0; i < num_bfs_roots; ++i) secs_per_edge[i] = sssp_times[i] / edge_counts[i];
			get_statistics(secs_per_edge, num_bfs_roots, stats);
			fprintf(stdout, "sssp min_TEPS:                  %g\n", 1. / stats[s_maximum]);
			fprintf(stdout, "sssp firstquartile_TEPS:        %g\n", 1. / stats[s_thirdquartile]);
			fprintf(stdout, "sssp median_TEPS:               %g\n", 1. / stats[s_median]);
			fprintf(stdout, "sssp thirdquartile_TEPS:        %g\n", 1. / stats[s_firstquartile]);
			fprintf(stdout, "sssp max_TEPS:                  %g\n", 1. / stats[s_minimum]);
			fprintf(stdout, "sssp harmonic_mean_TEPS:     !  %g\n", 1. / stats[s_mean]);
			fprintf(stdout, "sssp harmonic_stddev_TEPS:      %g\n", stats[s_std] / (stats[s_mean] * stats[s_mean] * sqrt(num_bfs_roots - 1)));
#endif
			free(secs_per_edge); secs_per_edge = NULL;
			free(edge_counts); edge_counts = NULL;
			get_statistics(validate_times, num_bfs_roots, stats);
			fprintf(stdout, "bfs  min_validate:              %g\n", stats[s_minimum]);
			fprintf(stdout, "bfs  firstquartile_validate:    %g\n", stats[s_firstquartile]);
			fprintf(stdout, "bfs  median_validate:           %g\n", stats[s_median]);
			fprintf(stdout, "bfs  thirdquartile_validate:    %g\n", stats[s_thirdquartile]);
			fprintf(stdout, "bfs  max_validate:              %g\n", stats[s_maximum]);
			fprintf(stdout, "bfs  mean_validate:             %g\n", stats[s_mean]);
			fprintf(stdout, "bfs  stddev_validate:           %g\n", stats[s_std]);
#ifdef SSSP
			get_statistics(validate_times2, num_bfs_roots, stats);
			fprintf(stdout, "sssp min_validate:              %g\n", stats[s_minimum]);
			fprintf(stdout, "sssp firstquartile_validate:    %g\n", stats[s_firstquartile]);
			fprintf(stdout, "sssp median_validate:           %g\n", stats[s_median]);
			fprintf(stdout, "sssp thirdquartile_validate:    %g\n", stats[s_thirdquartile]);
			fprintf(stdout, "sssp max_validate:              %g\n", stats[s_maximum]);
			fprintf(stdout, "sssp mean_validate:             %g\n", stats[s_mean]);
			fprintf(stdout, "sssp stddev_validate:           %g\n", stats[s_std]);
#endif
#if 0
			for (i = 0; i < num_bfs_roots; ++i) {
				fprintf(stdout, "Run %3d:                        %g s, validation %g s\n", i + 1, bfs_times[i], validate_times[i]);
				fprintf(stdout, "Run %3d:                        %g s, validation %g s\n", i + 1, sssp_times[i], validate_times2[i]);
			}
#endif


		}
	}
	free(bfs_times);
	free(validate_times);
#ifdef SSSP
	free(sssp_times);
	free(validate_times2);
#endif
	cleanup_globals();
	aml_finalize(); //includes MPI_Finalize()
	return 0;
}
