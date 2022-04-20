#include "../generator/utils.h"
#include "statistics.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>

#define OUT(fp) (fp == NULL)? stdout: fp

static int compare_doubles(const void* a, const void* b) 
{
	double aa = *(const double*)a;
	double bb = *(const double*)b;
	return (aa < bb) ? -1 : (aa == bb) ? 0 : 1;
}

void get_statistics(const double x[], int n, volatile double r[s_LAST]) 
{
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

void get_all_statistics(config_t graph500_cfg, stat_t stat, FILE* fp)
{
    volatile double stats[s_LAST];
    get_statistics(stat.bfs_times, graph500_cfg.num_bfs_roots, stats);
	fprintf(OUT(fp), "bfs  min_time:                  %g\n", stats[s_minimum]);
	fprintf(OUT(fp), "bfs  firstquartile_time:        %g\n", stats[s_firstquartile]);
	fprintf(OUT(fp), "bfs  median_time:               %g\n", stats[s_median]);
	fprintf(OUT(fp), "bfs  thirdquartile_time:        %g\n", stats[s_thirdquartile]);
	fprintf(OUT(fp), "bfs  max_time:                  %g\n", stats[s_maximum]);
	fprintf(OUT(fp), "bfs  mean_time:                 %g\n", stats[s_mean]);
	fprintf(OUT(fp), "bfs  stddev_time:               %g\n", stats[s_std]);
    
    get_statistics(stat.edge_counts, graph500_cfg.num_bfs_roots, stats);
	fprintf(OUT(fp), "min_nedge:                      %.11g\n", stats[s_minimum]);
	fprintf(OUT(fp), "firstquartile_nedge:            %.11g\n", stats[s_firstquartile]);
	fprintf(OUT(fp), "median_nedge:                   %.11g\n", stats[s_median]);
	fprintf(OUT(fp), "thirdquartile_nedge:            %.11g\n", stats[s_thirdquartile]);
	fprintf(OUT(fp), "max_nedge:                      %.11g\n", stats[s_maximum]);
	fprintf(OUT(fp), "mean_nedge:                     %.11g\n", stats[s_mean]);
	fprintf(OUT(fp), "stddev_nedge:                   %.11g\n", stats[s_std]);
	
    double* secs_per_edge = (double*)xmalloc(graph500_cfg.num_bfs_roots * sizeof(double));
	for (int i = 0; i < graph500_cfg.num_bfs_roots; ++i) secs_per_edge[i] = stat.bfs_times[i] / stat.edge_counts[i];
	get_statistics(secs_per_edge, graph500_cfg.num_bfs_roots, stats);
	fprintf(OUT(fp), "bfs  min_TEPS:                  %g\n", 1. / stats[s_maximum]);
	fprintf(OUT(fp), "bfs  firstquartile_TEPS:        %g\n", 1. / stats[s_thirdquartile]);
	fprintf(OUT(fp), "bfs  median_TEPS:               %g\n", 1. / stats[s_median]);
	fprintf(OUT(fp), "bfs  thirdquartile_TEPS:        %g\n", 1. / stats[s_firstquartile]);
	fprintf(OUT(fp), "bfs  max_TEPS:                  %g\n", 1. / stats[s_minimum]);
	fprintf(OUT(fp), "bfs  harmonic_mean_TEPS:     !  %g\n", 1. / stats[s_mean]);
	fprintf(OUT(fp), "bfs  harmonic_stddev_TEPS:      %g\n", stats[s_std] / (stats[s_mean] * stats[s_mean] * sqrt(graph500_cfg.num_bfs_roots - 1)));

    free(secs_per_edge); secs_per_edge = NULL;		
    get_statistics(stat.validate_times, graph500_cfg.num_bfs_roots, stats);
	fprintf(OUT(fp), "bfs  min_validate:              %g\n", stats[s_minimum]);
	fprintf(OUT(fp), "bfs  firstquartile_validate:    %g\n", stats[s_firstquartile]);
	fprintf(OUT(fp), "bfs  median_validate:           %g\n", stats[s_median]);
	fprintf(OUT(fp), "bfs  thirdquartile_validate:    %g\n", stats[s_thirdquartile]);
	fprintf(OUT(fp), "bfs  max_validate:              %g\n", stats[s_maximum]);
	fprintf(OUT(fp), "bfs  mean_validate:             %g\n", stats[s_mean]);
	fprintf(OUT(fp), "bfs  stddev_validate:           %g\n", stats[s_std]);
}