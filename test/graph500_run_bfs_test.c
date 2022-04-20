#include "generate_graph.h"
#include "common.h"
#include "aml.h"
#include "csr_reference.h"
#include "statistics.h"
#include "utils.h"
#include <stdio.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <assert.h>

static void graph500_generate_nonisolated_roots_test1()
{
    int64_t nglobalverts = (int64_t)(1) << SCALES[0];
	int64_t* bfs_roots = NULL;
	uint64_t seed1 = 2, seed2 = 3;
	double make_graph_time = 0;
	double data_struct_time = 0;
    char filename[64];
    tuple_graph tg;
    config_t cfg;

    init_graph500_cfg(&cfg, SCALES[0], edgefactors[0], 64, seed1, seed2);
    InitTupleGraph(cfg, &tg);
   
	make_graph_time = GenerateGraph(SCALES[0], seed1, seed2, nglobalverts, &tg);
    data_struct_time = make_graph_data_structure(&tg);
	assert_int_equal(GenerateNonIsolatedRoots(seed1, seed2, nglobalverts, cfg.num_bfs_roots, &bfs_roots), cfg.num_bfs_roots);    
	
    
	/* Run BFS. */
    stat_t bfs_statis; 
	uint64_t nlocalverts = get_nlocalverts_for_pred();
	float* shortest = (float*)xMPI_Alloc_mem(nlocalverts * sizeof(float));
    int64_t* pred = (int64_t*)xMPI_Alloc_mem(nlocalverts * sizeof(int64_t));
	run_bfs_for_roots(&tg, cfg.num_bfs_roots, bfs_roots, &bfs_statis, nlocalverts, pred, shortest);    
    
    sprintf(filename, "out/run_bfs_test_scale_%d_edgefactor_%d.log", SCALES[0], edgefactors[0]);
    FILE *fp = fopen(filename, "w+");
    get_all_statistics(cfg, bfs_statis, fp);
}

static void graph500_generate_nonisolated_roots_test2()
{

}

int main()
{    
    int res = -1; 

    aml_init(NULL, NULL);    
    setup_globals();    
    
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(graph500_generate_nonisolated_roots_test1),
        cmocka_unit_test(graph500_generate_nonisolated_roots_test2),
    };

    res = cmocka_run_group_tests(tests, NULL, NULL);
    cleanup_globals();
	aml_finalize(); //includes MPI_Finalize()

    return 0;
}