/* vim: set sts=0 ts=4 sw=4 noet : */

#ifndef _FILL_MAT_GRAPH_BUILDER_H_
#define _FILL_MAT_GRAPH_BUILDER_H_

#include "fill_config.h"

#include "MATGraph.h"

#if FILL_MA_CGAL
void buildMATGraphWithCGAL(MATGraph& mat, const Paths & inputPoly);
#endif
#if FILL_MA_BOOST
void buildMATGraphWithBOOST(MATGraph& mat, const Paths & inputPoly, const double units_per_mm);
#endif

#endif // _FILL_MAT_GRAPH_BUILDER_H_
