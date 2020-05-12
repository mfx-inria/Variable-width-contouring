/* vim: set sts=0 ts=4 sw=4 noet : */

#ifndef _CAVITIES_CAIRO_SURFACE_H_
#define _CAVITIES_CAIRO_SURFACE_H_

#include "MATGraph.h"

#include <string>
#include <vector>

#include "utils.h"
#include "Sampling.h"

using namespace std;

class CairoSurface {
	public:
	
	CairoSurface(const BBox & bbox, const string & filename, double min_radius, double max_radius,
			bool tight, bool recording=false);
	~CairoSurface();
	
	static void fuse(const BBox & bbox,
			CairoSurface & contours,
			CairoSurface & fused,
			MATGraph & mat,
			MATGraph::Components & comps,
			const vector<vector<Sample>> & pts,
			double min_radius,
			double max_radius,
			bool drawBisector
			);
	
	void clear();
	void polygonPath(const Paths & paths);
	void drawPolygons(const Paths & paths);
	void setDrawGrid(bool b) { drawGrid_ = b; }
	void setDrawVertices(bool b) { drawVertices_ = b; }
	void setUnderfill(bool b) { voidArea_ = b; }
	void setDrawMAVertices(bool b) { drawMAVertices_ = b; }
	void setLineWidthMultiplier(double m) { lineWidthMultiplier_ = m; }

	void drawTools(const BBox & bbox);

	void drawSmoothPath(const SmoothPath & sp);
	void drawSmoothPaths(const SmoothPaths & sp);
	void drawMA(const MATGraph& mat, bool noShave = false);

	void drawContour(const MATGraph::Component &);
	void drawContour(const SmoothPaths &);
	void drawDisks(const vector<Disk> & disks);
	void drawPaths(const vector<vector<Sample>> & pts);
	void step(int i = -1) { if( i == -1 ) { step_ = 1 - step_; } else { step_ = (abs(i)%2); } }
	void drawGrid(double grid_step = 10.0, double line_width = 0.3, bool skip10 = false);
	//protected:
	double min_radius, max_radius;
	int step_;
	cairo_t * context_;
	cairo_surface_t * surface_;
	double canvasWidth_, canvasHeight_, scale_;
	BBox bbox_;
	bool drawVertices_;
	bool drawMAVertices_;
	bool drawGrid_;
	bool tight_, voidArea_;
	double lineWidthMultiplier_;
};

#endif // _CAVITIES_CAIRO_SURFACE_H_
