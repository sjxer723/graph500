#ifndef GENERATE_GRAPH_H
#define GENERATE_GRAPH_H

#include "../generator/make_graph.h"
#include "../generator/utils.h"
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


double GenerateGraph(int SCALE, uint64_t seed1, uint64_t seed2, int64_t nglobalverts, tuple_graph* tg);

#endif

