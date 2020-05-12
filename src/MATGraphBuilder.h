/* vim: set sts=0 ts=4 sw=4 noet : */

#ifndef _FILL_MAT_GRAPH_BUILDER_H_
#define _FILL_MAT_GRAPH_BUILDER_H_

#include "MATGraph.h"

void buildMATGraphWithCGAL(MATGraph& mat, const Paths & inputPoly);
void buildMATGraphWithBOOST(MATGraph& mat, const Paths & inputPoly);

#endif // _FILL_MAT_GRAPH_BUILDER_H_
