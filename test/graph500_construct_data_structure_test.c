#include "generate_graph.h"
#include "common.h"
#include "aml.h"
#include "utils.h"
#include <stdio.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>
#include <stdint.h>
#include <inttypes.h>

static void graph500_construct_data_structure_test1()
{   
    int SCALE = SCALES[0];
    uint64_t edgefactor = edgefactors[0];
    uint64_t seed1 = 2, seed2 = 3;
    int64_t nglobalverts = (int64_t)(1) << SCALE;
    tuple_graph tg;

    InitTupleGraph(SCALE, seed1, seed2, edgefactor, &tg);
    GenerateGraph(SCALE, seed1, seed2, nglobalverts, &tg);

    double data_struct_start = MPI_Wtime();
	make_graph_data_structure(&tg);
	double data_struct_stop = MPI_Wtime();
	double data_struct_time = data_struct_stop - data_struct_start;
	fprintf(stderr, "construction_time:              %f s\n", data_struct_time);
}

static void graph500_construct_data_structure_test2()
{   
    int SCALE = SCALES[1];
    uint64_t edgefactor = edgefactors[1];
    uint64_t seed1 = 2, seed2 = 3;
    int64_t nglobalverts = (int64_t)(1) << SCALE;
    tuple_graph tg;

    InitTupleGraph(SCALE, seed1, seed2, edgefactor, &tg);
    GenerateGraph(SCALE, seed1, seed2, nglobalverts, &tg);

    double data_struct_start = MPI_Wtime();
	make_graph_data_structure(&tg);
	double data_struct_stop = MPI_Wtime();
	double data_struct_time = data_struct_stop - data_struct_start;
	fprintf(stderr, "construction_time:              %f s\n", data_struct_time);
}

int main()
{   
    int res = -1; 

    aml_init(NULL, NULL);    
    setup_globals();

    const struct CMUnitTest tests[] = {
        cmocka_unit_test(graph500_construct_data_structure_test1),
        cmocka_unit_test(graph500_construct_data_structure_test2),
    };

    res = cmocka_run_group_tests(tests, NULL, NULL);
    cleanup_globals();
	aml_finalize(); //includes MPI_Finalize()
   
    return res;
}