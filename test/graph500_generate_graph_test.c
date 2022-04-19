#include "generate_graph.h"
#include "common.h"
#include "utils.h"
#include <stdio.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>
#include <stdint.h>
#include <inttypes.h>

/* Test cases that generate a graph. */
static void graph500_generate_graph1(void **state) {
    int SCALE = SCALES[0];
	int edgefactor = edgefactors[0]; /* nedges / nvertices, i.e., 2*avg. degree */
	uint64_t seed1 = 2, seed2 = 3;
    double make_graph_time = 0;
    
    make_graph_time = GenerateGraph_WithoutWritingToFile(SCALE, seed1, seed2, edgefactor);
    printf("Generate graph 500 with SCALE = %d, edgefactor = %d costs %f\n", SCALE, edgefactor, make_graph_time);
}

static void graph500_generate_graph2(void **state) {
    int SCALE = SCALES[1];
	int edgefactor = edgefactors[1]; /* nedges / nvertices, i.e., 2*avg. degree */
	uint64_t seed1 = 2, seed2 = 3;
    double make_graph_time = 0;
    
    make_graph_time = GenerateGraph_WithoutWritingToFile(SCALE, seed1, seed2, edgefactor);
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