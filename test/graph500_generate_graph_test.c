#include "generate_graph.h"
#include "common.h"
#include <stdio.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>
#include <stdint.h>
#include <inttypes.h>

/* Test cases that generate a graph. */
static void graph500_generate_graph1(void **state) {
    int SCALE = 1;
	int edgefactor = 3; /* nedges / nvertices, i.e., 2*avg. degree */
	uint64_t seed1 = 2, seed2 = 3;
    double make_graph_time = 0;
    tuple_graph tg;
	
    tg.nglobaledges = (int64_t)(edgefactor) << SCALE;
    tg.write_file = 0;
    int64_t nglobalverts = (int64_t)(1) << SCALE;

    make_graph_time = GenerateGraph(SCALE, seed1, seed2, nglobalverts, &tg);
    printf("Generate graph 500 with SCALE = %d, edgefactor = %d costs %f\n", SCALE, edgefactor, make_graph_time);
}

static void graph500_generate_graph2(void **state) {
    int SCALE = 16;
	int edgefactor = 16; /* nedges / nvertices, i.e., 2*avg. degree */
	uint64_t seed1 = 2, seed2 = 3;
    double make_graph_time = 0;
    tuple_graph tg;
	
    tg.nglobaledges = (int64_t)(edgefactor) << SCALE;
    tg.write_file = 0;
    int64_t nglobalverts = (int64_t)(1) << SCALE;

    make_graph_time = GenerateGraph(SCALE, seed1, seed2, nglobalverts, &tg);
    printf("Generate graph 500 with SCALE = %d, edgefactor = %d costs %f\n", SCALE, edgefactor, make_graph_time);
}

int main()
{
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(graph500_generate_graph1),
        cmocka_unit_test(graph500_generate_graph2),
    };

    return cmocka_run_group_tests(tests, NULL, NULL);
}