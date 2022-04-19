#ifndef GENERATE_GRAPH_H
#define GENERATE_GRAPH_H

#include "../generator/make_graph.h"
#include "../generator/utils.h"
#include "common.h"

void InitTupleGraph (int SCALE, uint64_t seed1, uint64_t seed2,  int edgefactor, tuple_graph *tg);

double GenerateGraph(int SCALE, uint64_t seed1, uint64_t seed2, int64_t nglobalverts, tuple_graph* tg);

double GenerateGraph_WithoutWritingToFile(int SCALE, uint64_t seed1, uint64_t seed2, int edgefactor);

#endif

