#ifndef STATISTICS_H
#define STATISTICS_H

#include <stdint.h>
#include <stdio.h>

enum {
    s_minimum,
    s_firstquartile,
    s_median,
    s_thirdquartile,
    s_maximum,
    s_mean,
    s_std,
    s_LAST
};

typedef struct graph500_configure {
    int SCALE;
    int edgefactor; /* nedges / nvertices, i.e., 2*avg. degree */
	int num_bfs_roots;
    uint64_t seed1;
    uint64_t seed2;
} config_t;

typedef struct statistics {
    int validation_passed;
    double* edge_counts;
    double* bfs_times;
    double* validate_times;
} stat_t;

inline void init_graph500_cfg(config_t *cfg, int SCALE, int edgefactor, int num_bfs_roots, uint64_t seed1, uint64_t seed2)
{
    cfg -> SCALE = SCALE;
    cfg -> edgefactor = edgefactor;
    cfg -> num_bfs_roots = num_bfs_roots;
    cfg -> seed1 = seed1;
    cfg -> seed2 = seed2;
}

void get_statistics(const double x[], int n, volatile double r[s_LAST]);

void get_all_statistics(config_t graph500_cfg, stat_t stat, FILE* fp);

#endif
