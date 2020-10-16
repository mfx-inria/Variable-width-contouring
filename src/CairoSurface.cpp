/* vim: set sts=0 ts=4 sw=4 noet : */

#include "CairoSurface.h"

extern bool gInnerUnderfill, gOuterUnderfill, gUnderfill;
extern int gTotalSteps;
extern string gToolPos;

void pathColor(int step_, bool voidArea_, cairo_t * context_) {
	if( voidArea_ ) {
		cairo_set_source_rgb(context_,1,1,1);
		return;
	}
	if( 1 == step_ ) {
		cairo_set_source_rgb(context_, 1.0, 1.0, 1.0);
	} else {
		cairo_set_source_rgb(context_, 0.4, 0.4, 0.4);
	}
}

void errorColor(bool voidArea_, cairo_t * context_) {
	if( voidArea_ )
		cairo_set_source_rgb(context_, 0, 0, 0);
	else {
		cairo_set_source_rgb(context_, 1.0, 170.0/255.0, 201.0/255.0);
		//cairo_set_source_rgb(context_, 1, 0.3, 0.3);
	}
}

void trackColor(int step_, bool voidArea_, cairo_t * context_) {
	if( voidArea_ ) {
		cairo_set_source_rgb(context_,1,1,1);
		return;
	}
	if( step_ == 1 ) {
		cairo_set_source_rgb(context_, 0.70, 0.70, 1.0);
	} else {
		cairo_set_source_rgb(context_, 0.85, 0.85, 1.0);
	}
}

CairoSurface::
CairoSurface(const BBox & bbox, const string & filename, double min_radius, double max_radius,
		bool tight, bool recording)
: min_radius(min_radius)
, max_radius(max_radius)
, bbox_(bbox)
{
	const double W = 180 /* mm */  * 2.8346457; // 1 mm = 2.8346457 pt
	canvasWidth_  = bbox.max_.x() - bbox.min_.x();
	canvasHeight_ = bbox.max_.y() - bbox.min_.y();
	const double margin = tight ? (gInnerUnderfill ? 0.0 : 0.05) : (0.05*canvasWidth_);
	canvasWidth_  += 2*margin;
	canvasHeight_ += 2*margin;
	const double H = W * canvasHeight_ / canvasWidth_; // in pt

	if( recording ) {
		cairo_rectangle_t rect({0,0,W,H});
		surface_ = cairo_recording_surface_create(CAIRO_CONTENT_COLOR_ALPHA, & rect);
	} else {
		surface_ = cairo_pdf_surface_create(filename.c_str(), W, H);
	}
	context_ = cairo_create(surface_);

	scale_ = W/canvasWidth_;
	const double s = scale_;
	cairo_translate(context_, (margin-bbox.min_.x())*s, (margin+bbox.max_.y())*s);
	cairo_scale(context_ , s, -s);

	cairo_set_line_cap(context_, CAIRO_LINE_CAP_BUTT);
	cairo_set_line_join(context_, CAIRO_LINE_JOIN_ROUND);
	cairo_set_line_width(context_, 1.0/s);
	cairo_set_source_rgb(context_,0,0,0);
	step_ = 0;
	drawVertices_ = false;
	drawMAVertices_ = false;
	tight_ = tight;
	lineWidthMultiplier_ = 1.0;
}

CairoSurface::
~CairoSurface() {
	cairo_destroy(context_);
	cairo_surface_destroy(surface_);
}

void
CairoSurface::
drawGrid(double grid_step, double lw, bool skip10) {
	double grid_x_min = floor(bbox_.min_.x() / grid_step) * grid_step;
	double grid_y_min = floor(bbox_.min_.y() / grid_step) * grid_step;
	double grid_x_max = ceil(bbox_.max_.x() / grid_step) * grid_step;
	double grid_y_max = ceil(bbox_.max_.y() / grid_step) * grid_step;
	cairo_save(context_);
	cairo_set_source_rgba(context_, 0, 0, 0, 0.5);
	cairo_set_line_width(context_, lw/scale_);
	for( double x = grid_x_min; x < grid_x_max+grid_step; x += grid_step ) {
		if( skip10 && (((int)x) % 10 == 0) ) continue;
		cairo_move_to(context_, x, grid_y_min);
		cairo_line_to(context_, x, grid_y_max);
	}
	for( double y = grid_y_min; y < grid_y_max+grid_step; y += grid_step ) {
		if( skip10 && (((int)y) % 10 == 0) ) continue;
		cairo_move_to(context_, grid_x_min, y);
		cairo_line_to(context_, grid_x_max, y);
	}
	cairo_stroke(context_);
	cairo_restore(context_);
}

void
CairoSurface::
fuse(const BBox & bbox,
		CairoSurface & contours,
		CairoSurface & fused,
		MATGraph & mat,
		MATGraph::Components & comps,
		const vector<vector<Sample>> & pts,
		double min_radius,
		double max_radius,
		bool drawBisector,
		int currentStep
		)
{
	CairoSurface ma(bbox, "", min_radius, max_radius, fused.tight_, true);
	ma.setDrawMAVertices(fused.drawMAVertices_);
	ma.setLineWidthMultiplier(fused.lineWidthMultiplier_);
	cairo_t * cr = cairo_create(fused.surface_);
	cairo_save (cr);
	cairo_set_source_rgba(cr, 0, 0, 0, 0);
	cairo_set_operator (cr, CAIRO_OPERATOR_SOURCE);
	cairo_paint (cr);
	cairo_restore (cr);
	
	if( fused.drawGrid_ ) {
		//ma.drawGrid(.1, 0.03);
		ma.drawGrid(1.0, 0.1, true);
		ma.drawGrid(10.0, 0.3);
	}
	ma.drawMA(mat);

	for( const MATGraph::Component & comp : comps ) {
		contours.drawContour(comp);
	}
	if( drawBisector) {
		contours.drawPaths(pts);
	}
	contours.step();
	if( (! gUnderfill) || (currentStep==gTotalSteps+1) ) {
		cairo_set_source_surface(cr, contours.surface_, 0, 0);
		cairo_paint(cr);
		cairo_set_source_surface(cr, ma.surface_, 0.0, 0.0);
		cairo_paint(cr);
		fused.drawTools(bbox);
		cairo_show_page(fused.context_);
	}
	cairo_destroy(cr);
}

void
CairoSurface::
clear() {
}

void
CairoSurface::
polygonPath(const Paths & paths) {
	for( const auto & path : paths ) {
		if( path.size() < 3 ) continue;
		const size_t N(path.size());
		cairo_move_to(context_, path[N-1].x(), path[N-1].y());
		for( size_t i(0); i < N; ++i ) {
			cairo_line_to(context_,path[i].x(), path[i].y());
		}
		cairo_close_path(context_);
	}
}

void
CairoSurface::
drawPolygons(const Paths & paths) {
	if( voidArea_ ) {
		drawVertices_ = false;
		if( gOuterUnderfill ) {
			cairo_set_source_rgb(context_, 0., 0., 0.);
		} else {
			return;
		}
	} else {
		cairo_set_source_rgb(context_, 0.8, 0.8, 0.8);
	}
	cairo_save(context_);
	cairo_set_line_width(context_, lineWidthMultiplier_/scale_);
	polygonPath(paths);
	cairo_set_fill_rule(context_, CAIRO_FILL_RULE_EVEN_ODD);
	cairo_fill(context_);
	if( drawVertices_ ) {
		cairo_set_line_cap(context_, CAIRO_LINE_CAP_ROUND);
		cairo_set_line_width(context_, lineWidthMultiplier_ * bbox_.diagonal() / 50.0);
		cairo_set_source_rgb(context_, 1.0, 0.1, 0.1);
		for( const auto & path : paths ) {
			if( path.size() < 3 ) continue;
			const size_t N(path.size());
			for( size_t i(0); i < N; ++i ) {
				cairo_move_to(context_, path[i].x(), path[i].y());
				cairo_line_to(context_, path[i].x(), path[i].y());
			}
		}
		cairo_stroke(context_);
	}
	cairo_restore(context_);
}

void
CairoSurface::
drawTools(const BBox & bbox) {
	double minR = min_radius;
	double maxR = max_radius;
	double x, y;
	double W = 2*(minR + maxR);
	double H = 2*maxR;
	x = 0.5 * (bbox.max_.x() + bbox.min_.x()) - W/2;
	y = 0.5 * (bbox.max_.y() + bbox.min_.y()) - H/2;
	for( char c : gToolPos ) {
		if( c == 'n' ) return;
		if( c == 'l' ) x = bbox_.min_.x();
		if( c == 'r' ) x = bbox_.max_.x()-W;
		if( c == 'b' ) {
			y = bbox_.min_.y();
			if( ! tight_ ) {
				y -= H;
			}
		}
		if( c == 't' ) y = bbox_.max_.y()-H;
	}
	x += minR;
	y += H-minR;
	cairo_set_line_width(context_, 2.0 * minR);
	cairo_move_to(context_, x, y); cairo_line_to(context_, x, y);
	cairo_set_line_cap(context_, CAIRO_LINE_CAP_ROUND);
	cairo_stroke(context_);
	x += minR + maxR;
	y += - maxR + minR;
	cairo_set_line_width(context_, 2.0 * maxR);
	cairo_move_to(context_, x, y);
	cairo_line_to(context_, x, y);
	cairo_stroke(context_);
}

void
CairoSurface::
drawSmoothPaths(const SmoothPaths & sps) {
	for( const auto & sp : sps )
		drawSmoothPath(sp);
}

void
CairoSurface::
drawSmoothPath(const SmoothPath & sp) {
	size_t N = sp.size();
	if( N == 0 ) return;
	if( N == 1 ) {
		const Vec2d c = sp[0].center();
		const double r =  sp[0].radius();
		cairo_move_to(context_, c.x()+r, c.y());
			cairo_arc(context_, c.x(), c.y(), r, 0.0, 2.0*M_PI);
		return;
	}
	Vec2d P, Q;
	BitangentComputer bmake(min_radius);
	bmake(sp[N-1], sp[0], P, Q);
	cairo_move_to(context_, P.x(), P.y());
	cairo_line_to(context_, Q.x(), Q.y());
	for( size_t i(0); i < N; ++i ) {
		size_t j = (i+1) % N;
		Vec2d p0, p1;
		bmake(sp[i], sp[j], p0, p1);
		// draw MATedge from Q to p0
		const Vec2d & c = sp[i].center();
		if( sp[i].radius() < 1e-5 ) {
			cairo_line_to(context_, c.x(), c.y());
		} else {
			Vec2d v = Q - c;
			double angleQ = atan2(v.y(), v.x());
			v = p0 - c;
			double angleP = atan2(v.y(), v.x());
			if( sp[i].passage_ == ToTheRight ) {
				while( angleP < angleQ ) angleP += 2.0*M_PI;
				if( angleP - angleQ > 2.0*M_PI-1e-3 ) {
					//cerr << "KABOUM at " << c << endl;
					std::swap(angleQ, angleP);
				}
				cairo_arc(context_, c.x(), c.y(), sp[i].radius(), angleQ, angleP);
			} else {
				while( angleQ < angleP ) angleQ += 2.0*M_PI;
				if( angleQ >= angleP + M_PI ) { // inversion due to numerical error in BitangentMaker
					angleQ -= 2.0*M_PI;
					std::swap(angleQ, angleP);
				}
				cairo_arc_negative(context_, c.x(), c.y(), sp[i].radius(), angleQ, angleP);
			}
		}
		Q = p1;
	}
}

void
CairoSurface::
drawMA(const MATGraph& mat, bool noShave) {
	cairo_save(context_);
	double lw = lineWidthMultiplier_*bbox_.diagonal()/500.0;
	const bool onlyTrimmed = true;
	const bool allOthers = false;
	mat.drawMedialAxis(context_, lw, drawMAVertices_, onlyTrimmed, noShave);
	mat.drawMedialAxis(context_, lw, drawMAVertices_, allOthers, noShave);
	cairo_restore(context_);
}

void
CairoSurface::
drawContour(const MATGraph::Component & comp) {
	cairo_save(context_);
	cairo_set_fill_rule(context_, CAIRO_FILL_RULE_EVEN_ODD);
	const SmoothPaths * outer = & comp.outerSmoothPaths;
	const SmoothPaths * inner = & comp.innerSmoothPaths;
	if( ! gOuterUnderfill ) {
		if( ! comp.voids.empty() ) {
			errorColor(voidArea_, context_);
			inner = & comp.voids;
			drawSmoothPaths(comp.voids);
			drawSmoothPaths(comp.innerSmoothPaths);
			cairo_fill(context_);
		}
		if( ! comp.simpleOuterSmoothPaths.empty() ) {
			errorColor(voidArea_, context_);
			outer = & comp.simpleOuterSmoothPaths;
			drawSmoothPaths(comp.simpleOuterSmoothPaths);
			drawSmoothPaths(comp.outerSmoothPaths);
			cairo_fill(context_);
		}
	}
	trackColor(step_, voidArea_, context_);
	drawSmoothPaths(*outer);
	if( ! gOuterUnderfill ) {
		drawSmoothPaths(*inner);
	}
	cairo_fill(context_);
	cairo_restore(context_);
}

void
CairoSurface::
drawContour(const SmoothPaths & sps) {
	trackColor(step_, voidArea_, context_);
	cairo_save(context_);
	drawSmoothPaths(sps);
	cairo_set_fill_rule(context_, CAIRO_FILL_RULE_EVEN_ODD);
	cairo_fill(context_);
	cairo_restore(context_);
}

void
CairoSurface::
drawDisks(const vector<Disk> & disks) {
	cairo_save(context_);
	cairo_set_line_cap(context_, CAIRO_LINE_CAP_ROUND);
	cairo_set_source_rgba(context_, 1, 0.6, 0.6, 0.3);
	for( const auto & d : disks ) {
		double x(d.center_.x());
		double y(d.center_.y());
		double r(d.radius_);
		cairo_move_to(context_, x+r, y);
		cairo_arc(context_, x, y, r, 0, 2*M_PI);
	}
	cairo_fill(context_);
	cairo_restore(context_);
}

void
CairoSurface::
drawPaths(const vector<vector<Sample>> & pts) {
	cairo_save(context_);
	double lw = lineWidthMultiplier_*bbox_.diagonal()/700.0;
	cairo_set_line_width(context_, lw);
	cairo_set_line_cap(context_, CAIRO_LINE_CAP_BUTT);
	for( auto & path : pts ) {
		const size_t N(path.size());
		if( N < 2 ) continue;
#if 0
		cairo_set_source_rgb(context_, 0.6, 0.3, 0.6);
		for( size_t i(0); i < N; i+=2 ) {
			cairo_move_to(context_, path[i+0].pos.x(), path[i+0].pos.y());
			cairo_line_to(context_, path[(i+1)%N].pos.x(), path[(i+1)%N].pos.y());
		}
		cairo_stroke(context_);
		cairo_set_source_rgb(context_, 0.5, 0.7, 0.1);
		for( size_t i(1); i < N; i+=2 ) {
			cairo_move_to(context_, path[i+0].pos.x(), path[i+0].pos.y());
			cairo_line_to(context_, path[(i+1)%N].pos.x(), path[(i+1)%N].pos.y());
		}
		cairo_stroke(context_);
		cairo_set_source_rgb(context_, 0,0,0);
#else
		//static const double dashed[] = {0.2};
		//static int dashlen = sizeof(dashed) / sizeof(dashed[0]);
		//cairo_set_dash(context_, dashed, dashlen, 0);
		pathColor(step_, voidArea_, context_);
		cairo_set_line_cap(context_, CAIRO_LINE_CAP_BUTT);
		cairo_set_line_join(context_, CAIRO_LINE_JOIN_ROUND);
		cairo_move_to(context_, path[0].pos.x(), path[0].pos.y());
		for( size_t i(1); i < N; ++i ) {
			cairo_line_to(context_, path[i].pos.x(), path[i].pos.y());
		}
		cairo_close_path(context_);
		cairo_stroke(context_);
#endif
#if 0
		// highlight first vertex
		cairo_set_line_cap(context_, CAIRO_LINE_CAP_ROUND);
		cairo_set_source_rgb(context_, 1.0, 0.8, 1.0);
		cairo_set_line_width(context_, lw * 2);
		cairo_move_to(context_, path[0].pos.x(), path[0].pos.y());
		cairo_line_to(context_, path[0].pos.x(), path[0].pos.y());
		cairo_stroke(context_);
		cairo_set_line_width(context_, lw);
#endif
	}
	cairo_restore(context_);
}
