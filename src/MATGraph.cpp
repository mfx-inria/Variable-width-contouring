/* vim: set sts=0 ts=4 sw=4 noet : */

#include "MATGraph.h"

#include <iterator> // for std::distance
#include <algorithm> // rotate
#include <stack>

void
CollapsedAxis::
jump() {
	if( subAxes.empty() ) return;
	curAP = farAP;
	//cerr << "Jump " << farAP.axis() << '/' << farAP.pos() << endl;
}

void
CollapsedAxis::
reorder() {
	assert( 1 == subAxes.size() );
	SubAxis & sub = subAxes[0];
	if( 0 >= farAP.pos() ) { // farPos is 0 or -1
		//cerr << "Reorder: overlap segments.\n";
		sub.points.push_back(sub.points[1]); // first and last segments are the same
		return;
	}
	assert( nullptr == sub.clippingVertex );
	//cerr << "Reordering (0, " << farPos << ", " << (sub.points.size()-1) << ')' << endl;
	sub.points.pop_back(); // open the loop
	std::rotate(sub.points.begin(), sub.points.begin()+farAP.pos(), sub.points.end());
	sub.points.push_back(sub.points[0]); // close the loop
	sub.points.push_back(sub.points[1]); // first and last segments are the same
	farAP.set(0,-1);
	curAP.set(0,-1);
}

bool
CollapsedAxis::
incr(AxisPos & ap) const {
	const SubAxis & axis = subAxes[ap.axis()];
	++ap.pos();
	if( ap.pos() == axis.points.size() - 1 ) {
		ap.pos() = -1;
		++ap.axis();
		if( ap.axis() == subAxes.size() ) { // go back to starting |ap|
			--ap.axis();
			ap.pos() = axis.points.size() - 1;
			return false;
		}
	}
	return true;
}

void
CollapsedAxis::
computeBoundaryCircles(SmoothPath & sp, Samples::const_iterator it, Samples::const_iterator end)
	const {
	//size_t nbSamples = std::distance(it, end);
	//cerr << nbSamples << " samples along the collapses axis." << endl;

	Samples::const_iterator lastDisk = end;
	Samples::const_iterator prevIt = it++;
	++it;
	for( ; it != end; ++it ) {
		AxisPos ap = it->axisPos;
		if( ap.axis() >= subAxes.size() ) {
			cerr << "Bad axis " << ap.axis() << " instead of <" << subAxes.size() << endl;
		}
		const SubAxis & axis = subAxes[ap.axis()];
		if( (ap.pos() >= 0) && (ap.pos() >= axis.points.size()) ) {
			cerr << "Bad pos " << ap.pos() << " instead of <" << axis.points.size() << endl;
		}
		AxisPos prevAP = prevIt->axisPos;
		if( prevAP < ap ) {
			if( -1 == prevAP.pos() || -1 == ap.pos() ) {
				if( nullptr != axis.clippingVertex ) {
					if( lastDisk != prevIt )
						sp.emplace_back(ToTheLeft, prevIt->pos, prevIt->radius);
					sp.emplace_back(ToTheLeft, it->pos, it->radius);
					lastDisk = it;
				}
			} else {
				if( prevAP.axis() == ap.axis() ) {
					int i = prevAP.pos() + 1;
					bool onlyLeftTurns = true;
					while( i <= ap.pos() ) {
						Vec2d a = axis.points[i] - axis.points[i-1];
						Vec2d b = axis.points[i+1] - axis.points[i];
						//if( det(a, b) == 0.0 ) cerr << "U-turn at " << axis.points[i] << endl;
						if( det(a, b) < 0.0 )
							onlyLeftTurns = false;
						++i;
					}
					if( onlyLeftTurns ) {
						while( prevAP != ap ) {
							incr(prevAP);
							sp.emplace_back(ToTheLeft, subAxes[prevAP.axis()].points[prevAP.pos()], 0);
							lastDisk = end;
						}
					} else {
						//sp.emplace_back(ToTheLeft, subAxes[prevAP.axis()].points[prevAP.pos()], 0);
						if( lastDisk != prevIt )
							sp.emplace_back(ToTheLeft, prevIt->pos, prevIt->radius);
						sp.emplace_back(ToTheLeft, it->pos, it->radius);
						lastDisk = it;
						//sp.emplace_back(ToTheLeft, subAxes[ap.axis()].points[ap.pos()+1], 0);
					}
				} else {
					const SubAxis & curAxis = subAxes[prevAP.axis()];
					cerr << "Changing axis! from " << prevAP << " to " << ap << " (at"
						<< curAxis.clippingVertex->pos() <<  " --> "
						<< axis.clippingVertex->pos() <<  ")" << endl;
				}
			}
		}
		prevIt = it;
	}
}

Sample
CollapsedAxis::
operator()(const Sample & sample) {
	bool bad(false), nothingYet(true);
	Sample best;
	AxisPos bestAP;
	auto update = [&]( const Sample & candidate, AxisPos axisPos) {
		if( bad ) {
			//cerr << '$' << pos << ' ';
			return;
		}
		if( nothingYet || (candidate.radius < best.radius) ) {
			best = candidate;
			bestAP = axisPos;
			// update to lexicographically largest
			if( farAP < axisPos ) farAP = axisPos;
		}
		nothingYet = false;
	};
	for( int ax(curAP.axis()); ax < subAxes.size(); ++ax ) {
		const SubAxis & subAxis = subAxes[ax];
		AxisPos ap;
		ap.axis() = ax;
		ap.pos() = (ax > curAP.axis()) ? -1 : curAP.pos();
		if( -1 == ap.pos() ) {
 			if( nullptr != subAxis.clippingVertex ) {
				Sample s = Sampling::moveTowardBisectorWithDisk(sample, 0, subAxis.clippingVertex->circumcircle, & bad);
				update(s, ap);
			}
			++ap.pos();
		}
		const vector<Vec2d> & points = subAxis.points;
		const int end = points.size() - 1;
		for( ; ap.pos() < end; ++ap.pos() ) {
			update(Sampling::moveTowardBisectorWithSegment(sample, 0, points[ap.pos()], points[ap.pos()+1], bad), ap);
		}
	}
	best.axisPos = bestAP;
	return best;
}

void
CollapsedAxis::
print() const {
	cerr << "{ AXIS }" << endl;
	for( const auto & sub : subAxes ) {
		if( nullptr != sub.clippingVertex ) {
			cerr << "Disk(" << sub.clippingVertex->pos() << ", " << sub.clippingVertex->radius() << "), ";
		} else {
			cerr << "No clipping disk, ";
		}
		cerr << sub.points.size() << " points: (";
		for( int i(0); i < sub.points.size(); ++i ) {
			cerr << i << ':' << sub.points[i] << ") (";
		}
	}
	cerr << "END)" << endl;
}

void
CollapsedAxis::
addPoint(const Vec2d & p) {
	assert( ! subAxes.empty());
	subAxes.back().points.push_back(p);
}

void
CollapsedAxis::
addSubAxis(MATvert * v) {
	subAxes.emplace_back();
	subAxes.back().clippingVertex = v;
}

void
CollapsedAxis::
clear() {
	subAxes.clear();
	farAP.set(0,-1);
	curAP.set(0,-1);
}

void
MATGraph::
computeConnectedComponents(Components & components) {
	components.clear();
	VertexSet visited;
	for ( MATvert & vert_ : verts ) {
		MATvert * vert = & vert_;
		if ( visited.count(vert) > 0 ) continue;
		components.emplace_back();
		VertexSet & component = components.back().vertexSet;
		std::vector< MATvert *> to_check;
		to_check.emplace_back(vert);
		component.emplace(vert);
		while ( ! to_check.empty() ) {
			MATvert * v = to_check.back(); to_check.pop_back();
			for ( int i = 0 ; i < 3 ; ++i ) {
				MATvert * neighbor = v->edge[i].to();
				if (v->has_neighbor(i) && component.count(neighbor) == 0) {
					component.emplace(neighbor);
					to_check.emplace_back(neighbor);
				}
			}
		}
		visited.insert(component.begin(), component.end());
	}
}

int
MATGraph::
statusDegree(const MATvert * vert, MatStatus status) const {
	MATedge * deadend, * freeway;
	return statusDegree(vert, status, freeway, deadend);
}

int
MATGraph::
statusDegree(const MATvert * vert, MatStatus status, MATedge *& freeway, MATedge *& deadend) const {
	int deg = 0;
	deadend = nullptr, freeway = nullptr;
	for( int i(0); i < 3; ++i ) {
		if( vert->has_neighbor(i) && (vert->edge[i].to()->status == status) ) {
			freeway = const_cast<MATedge *>(&vert->edge[i]);
			++deg;
		} else {
			deadend = const_cast<MATedge *>(&vert->edge[i]);
		}
	}
	switch( deg ) {
		case 1 : deadend = nullptr; break;
		case 2 : freeway = nullptr; break;
		default: deadend = freeway = nullptr;
	}
	assert(deg >= 0 && deg < 4);
	return deg;
}

int
MATGraph::
degree(const MATvert * vert, bool walk_on_fire) const {
	MATedge * freeway, * deadend;
	return degree(vert, freeway, deadend, walk_on_fire);
}

int
MATGraph::
degree(const MATvert * vert, MATedge *& freeway, MATedge *& deadend, bool walk_on_fire) const {
	int deg = 0;
	deadend = nullptr, freeway = nullptr;
	for( int i(0); i < 3; ++i ) {
		if( vert->has_neighbor(i) && (walk_on_fire || (vert->edge[i].to()->status == MatStatus::Normal)) ) {
			freeway = const_cast<MATedge *>(&vert->edge[i]);
			++deg;
		} else {
			deadend = const_cast<MATedge *>(&vert->edge[i]);
		}
	}
	switch( deg ) {
		case 1 : deadend = nullptr; break;
		case 2 : freeway = nullptr; break;
		default: deadend = freeway = nullptr;
	}
	assert(deg >= 0 && deg < 4);
	return deg;
}

bool
MATGraph::
componentHasCollapsedPart(const VertexSet & vs) const {
	std::unordered_set<const MATvert *> visited;
	for( const MATvert * vert : vs ) {
		if( visited.count(vert) > 0 ) continue;
		if( vert->is_collapsed() ) return true;
		std::stack<const MATvert *> jobs;
		jobs.push(vert);
		while( ! jobs.empty() ) {
			const MATvert * v = jobs.top(); jobs.pop();
			visited.insert(v);
			for( int i = 0; i < 3; ++i ) {
				if( ! v->edge[i].to() ) continue;
				if( v->edge[i].to()->is_collapsed() ) return true;
				if( visited.count(v->edge[i].to()) > 0 ) continue;
				jobs.push(v->edge[i].to());
			}
		}
	}
	return false;
}

MATvert&
MATGraph::
insert(Component * component, MATedge & edge, const Disk & circumcircle)
{
	debugCheck();
	MATedge * original_edge = & edge;
	MATedge * original_twin = original_edge->twin();
	verts.emplace_back(circumcircle);
	MATvert & vert = verts.back();
	if( nullptr != component ) {
		component->vertexSet.insert(&vert);
	}
	vert.edge[0].twin_ = original_edge;
	vert.edge[0].to_ = original_twin->to_;
	vert.edge[0].site_ = original_twin->site_;
	vert.edge[0].type = original_edge->type;

	vert.edge[1].twin_ = original_twin;
	vert.edge[1].to_ = original_edge->to_;
	vert.edge[1].site_ = original_edge->site_;
	vert.edge[1].type = original_twin->type;

	vert.edge[2].to_ = nullptr;
	vert.edge[2].twin_ = nullptr;

	original_edge->twin_ = & vert.edge[0];
	original_edge->to_ = &vert;
	original_twin->twin_ = & vert.edge[1];
	original_twin->to_ = &vert;

	debugCheck();
	return vert;
}

void
MATGraph::
removeVertexSet(const std::unordered_set<MATvert *> to_remove) {
	for ( auto vert_it = verts.begin(); vert_it != verts.end(); ) {
		MATvert * vert = &(*vert_it);
		if( to_remove.count(vert) > 0 ) {
			for (int i = 0; i < 3; ++i ) {
				if ( vert->edge[i].twin() ) {
					vert->edge[i].twin()->to_ = nullptr;
					vert->edge[i].twin()->twin_ = nullptr;
				}
			}
			vert_it = verts.erase(vert_it);
		} else {
			++vert_it;
		}
	}
	debugCheck();
}

void
MATGraph::
removeDestroyed() {
	for ( auto vert_it = verts.begin(); vert_it != verts.end(); ) {
		MATvert & vert = *vert_it;
		if (vert.status != MatStatus::Normal) {
			for (int i = 0; i < 3; ++i ) {
				if ( vert.edge[i].twin() ) {
					vert.edge[i].twin()->to_ = nullptr;
					vert.edge[i].twin()->twin_ = nullptr;
				}
			}
			vert_it = verts.erase(vert_it);
		} else {
			++vert_it;
		}
	}
	debugCheck();
}

MATedge*
MATedge::
cw() const
{
	assert(twin() && twin()->to());
	MATvert * from_ = from();
	int idx = this - static_cast<MATedge*>(from_->edge);
	for ( int i = 0; i < 2; ++i ) {
		idx--;
		if (idx == -1) idx = 2;
		if ( from_->edge[idx].to() ) return &from_->edge[idx];
	}
	return const_cast<MATedge*>(this);
}

MATedge*
MATedge::
ccw() const
{
	assert(twin() && twin()->to());
	MATvert * from_ = from();
	int idx = this - static_cast<MATedge*>(from_->edge);
	for ( int i = 0; i < 2; ++i ) {
		idx++;
		if (idx == 3) idx = 0;
		if ( from_->edge[idx].to() ) return &from_->edge[idx];
	}
	return const_cast<MATedge*>(this);
}

MATedge*
MATedge::
nextNotTrimmed() const {
	assert(to() && twin());
	MATedge* edge = twin()->ccw();
	for ( int i = 0; i < 2; ++i ) {
		if ( edge->to() && (!edge->to()->is_trimmed()) ) return edge;
		edge = edge->ccw();
	}
	return twin();
}

MATedge*
MATedge::
nextNormalOrTwin() const {
	assert(to() && twin());
	MATedge* edge = twin()->ccw();
	for ( int i = 0; i < 2; ++i ) {
		if ( edge->to() && (!edge->to()->is_destroyed()) ) return edge;
		edge = edge->ccw();
	}
	// or twin:
	return twin();
}

MATedge*
MATedge::
nextNormalOrNextAny() const {
	assert(to() && twin());
	MATedge* edge = twin()->ccw();
	for ( int i = 0; i < 3; ++i ) {
		if ( edge->to() && (!edge->to()->is_destroyed()) ) return edge;
		edge = edge->ccw();
	}
	// or next any:
	return twin()->ccw();
}

MATedge*
MATedge::
next(bool walk_on_fire) const {
	assert(to() && twin());
	assert(walk_on_fire || (!from()->is_destroyed() && !to()->is_destroyed()));
	MATedge* edge = twin()->ccw();
	for ( int i = 0; i < 2; ++i ) {
		if ( edge->to() && (walk_on_fire || !edge->to()->is_destroyed()) ) return edge;
		edge = edge->ccw();
	}
	return twin();
}

MATedge*
MATedge::
prev(bool walk_on_fire) const {
	assert(twin());
	assert(walk_on_fire || (!from()->is_destroyed() && !to()->is_destroyed()));
	MATedge* edge = cw();
	for ( int i = 0; i < 2; ++i ) {
		if ( edge->to() && (walk_on_fire || !edge->to()->is_destroyed()) ) return edge->twin();
		edge = edge->cw();
	}
	return twin();
}

Vec2d
MATedge::
getGeneratorLocation() const {
	if (site().is_point()) return site().point();
	assert(site().is_segment());
	assert(twin() && twin()->to());
	Vec2d p = from()->circumcircle.center_;
	Vec2d a = site().source();
	Vec2d b = site().destination();
	Vec2d ab = b - a;
	Vec2d ap = p - a;
	Vec2d ret = a + (ab.dot(ap) / ab.length()) * ab.normalized();
	return ret;
}

double
MATedge::
getRadiusGradientAtStart() const {
	Vec2d right = getGeneratorLocation();
	Vec2d left = twin()->next(true)->getGeneratorLocation();
	Vec2d p = from()->circumcircle.center_;
	Vec2d l = (left-p).normalized();
	Vec2d r = (p-right).normalized();
	double cos2a = l.dot(r);
	double squaredSina = max(0.0, (1.0-cos2a)/2.0);
	double sina = sqrt(squaredSina);
	if( det(l, r) < 0.0 )
		return sina;
	else
		return -sina;
}

#ifndef FILL_NO_CAIRO
// A Cairo helper function that draw a quadratic bezier curve with a cubic one.
void
MATGraph::
helper_quadratic_to(cairo_t *cr,
                    double x1, double y1,
                    double x2, double y2) const
{
	double x0, y0;
	cairo_get_current_point (cr, &x0, &y0);
	cairo_curve_to (cr,
	                2.0 / 3.0 * x1 + 1.0 / 3.0 * x0,
	                2.0 / 3.0 * y1 + 1.0 / 3.0 * y0,
	                2.0 / 3.0 * x1 + 1.0 / 3.0 * x2,
	                2.0 / 3.0 * y1 + 1.0 / 3.0 * y2,
	                x2, y2);
}

// draw medial axis to a cairo context
bool
MATGraph::
shaveConnectedToNormal(const MATvert * v, const MATvert * from) const {
	bool r = false;
	for( int i = 0; i < 3; ++i ) {
		const MATedge & edge = v->edge[i];
		if (!edge.to()) continue;
		const MATvert * neighbor = edge.to();
		if( neighbor->is_normal() ) return true;
		if( ! neighbor->is_shaved() ) continue;
		if( neighbor == from ) continue;
		r = r || shaveConnectedToNormal(neighbor, v);
	}
	return r;
}

static Vec4d gNormalColor(0.6, 1.0, 0.6, 1.0);
static Vec4d gTrimmedColor(0.0, 0.0, 0.0, 1.0);
static Vec4d gCollapsedColor(1.0, 0.0, 0.0, 1.0);
static Vec4d gShavedColor(1.0, 1.0, 0.2, 1.0);

// draw medial axis to a cairo context
void
MATGraph::
drawMedialAxis(cairo_t * context, double lw, bool drawVertices,
		double maxCollapseRadius, bool onlyTrimmed, bool noShave) const {
	cairo_set_line_width(context, lw);
	cairo_set_line_cap(context, CAIRO_LINE_CAP_BUTT);
	cairo_set_source_rgba(context,1,1,1,1);
	for( const MATvert & vert : verts ) {
		int nbArcs(0);
		for( int i = 0; i < 3; ++i ) {
			const MATedge & edge = vert.edge[i];
			if (!edge.to()) continue;
			++nbArcs;
			const MATvert& neighbor = *edge.to();
			if( neighbor < vert ) continue;
			// OK, draw the arc
			Vec2d a = vert.circumcircle.center_;
			Vec2d b = neighbor.circumcircle.center_;
			Vec4d rgba;
			bool trimmedEdge(false);
			auto colorEdge = [&](const MATedge & edge) {
				if( edge.to()->is_collapsed() || edge.from()->is_collapsed()) rgba = gCollapsedColor;
				if( edge.to()->is_shaved() || edge.from()->is_shaved()) {
					rgba = gShavedColor;
					const MATvert * shavedVert = edge.to();
					if( edge.from()->is_shaved() ) shavedVert = edge.from();
					if( shaveConnectedToNormal(shavedVert) )
						rgba = gCollapsedColor;
				}
				if( edge.to()->is_trimmed() || edge.from()->is_trimmed()) {
					rgba = gTrimmedColor;
					trimmedEdge = true;
				}
				if( (!noShave) && ((edge.to()->is_normal() && edge.from()->is_collapsed()) ||
						(edge.to()->is_collapsed() && edge.from()->is_normal())) )
					rgba = gCollapsedColor;//Vec4d(1,0.9,0.2,1); // INNER SHAVE
				cairo_set_source_rgba(context, rgba[0], rgba[1], rgba[2], rgba[3]);
			};
			rgba = gNormalColor;
			if( edge.type != EdgeType::VertEdge ) {
				colorEdge(edge);
				if( onlyTrimmed != trimmedEdge ) continue;
				cairo_move_to(context, a.x(), a.y());
				cairo_line_to(context, b.x(), b.y());
			} else {
				colorEdge(edge);
				if( onlyTrimmed != trimmedEdge ) continue;
				QuadraticArc qArc;
				buildBezierQuadraticOfArc(edge, qArc);
				cairo_move_to(context, a.x(), a.y());
				helper_quadratic_to(context, qArc.b.x(), qArc.b.y(), b.x(), b.y());
			}
			cairo_stroke(context);
		}
		if( 0 == nbArcs ) {
			Vec2d b = vert.circumcircle.center_;
			cairo_move_to(context, b.x(), b.y());
			cairo_line_to(context, b.x(), b.y());
			cairo_set_line_cap(context, CAIRO_LINE_CAP_ROUND);
			cairo_stroke(context);
			cairo_set_line_cap(context, CAIRO_LINE_CAP_BUTT);
		}
	}
	// draw medial axis vertices
	if( onlyTrimmed ) return;
	cairo_set_line_cap(context, CAIRO_LINE_CAP_ROUND);
	const double pointDiameter = lw * 2.5;
	cairo_set_line_width(context, pointDiameter);
	for( const MATvert & vert : verts ) {
		if( ! drawVertices ) { // draw only full-degree 0 Collapsed and Normal
			if( ! ( ( vert.is_normal() && 0 == degree(& vert) ) ||
						( vert.is_collapsed() && 0 == statusDegree(& vert, MatStatus::Collapsed) ) ) )
				continue;
		}
		const Disk & d = vert.circumcircle;
		Vec2d p = d.center_;
		Vec4d rgba;
		bool drawCircle(false);
		switch( vert.status ) {
			case MatStatus::Normal: rgba = gNormalColor; drawCircle = true; break;
			case MatStatus::Trimmed: rgba = gTrimmedColor; break;
			case MatStatus::Collapsed: rgba = gCollapsedColor; break;
			case MatStatus::Shaved: rgba = gShavedColor; break;
			default: break;
		}
		cairo_set_source_rgba(context, rgba[0], rgba[1], rgba[2], rgba[3]);
		cairo_move_to(context, p.x(), p.y());
		cairo_line_to(context, p.x(), p.y());
		cairo_stroke(context);
		if( false ) {//drawCircle ) {
			static const double dashed1[] = {0.5*lw, 0.8*lw};
			static int len1  = sizeof(dashed1) / sizeof(dashed1[0]);
			cairo_set_dash(context, dashed1, len1, 0);
			cairo_set_line_width(context, lw*0.15);
			cairo_set_source_rgba(context,0,0,0,1);
			cairo_move_to(context, p.x()+d.radius_, p.y());
			cairo_arc(context, p.x(), p.y(), d.radius_, 0.0, 2*M_PI);
			cairo_stroke(context);
			cairo_set_dash(context, nullptr, 0, 0);
			cairo_set_line_width(context, pointDiameter);
		}
	}
}
#endif

void
MATGraph::
buildBezierQuadraticOfArc(const MATedge & edge, QuadraticArc & qArc) const {
	// Get |p| and |q| the two elements whose bisector is the Voronoi
	// arc that we want to traverse
	const Site & p = edge.site();
	const Site & q = edge.twin()->site();
	bool isVV = edge.type == EdgeType::VertVert;
	bool isArcAStraightLine = isVV || ( p.is_segment() && q.is_segment() );
	qArc.a = edge.from()->circumcircle.center_;
	qArc.c = edge.to()->circumcircle.center_;
	// set middle of Quadratic to the mid-point of 'a' abd 'c'
	qArc.b = Vec2d(0.5*(qArc.a.x()+qArc.c.x()), 0.5*(qArc.a.y()+qArc.c.y()));
	if( isArcAStraightLine ) {
		return;
	}
	PointSegmentBisector psb;
	makePSBisector(edge, psb);
	// BEGIN Compute the antenna point for the quadratic arc
	double dx = psb.x0 - psb.x1;
	Vec2d v = (-dx*dx/4.0/psb.h) * psb.v;
	qArc.b = qArc.b + v;
	return;
}


// fills the psb struct. See doc of this struct at its definition (cgal.h)
void
MATGraph::
makePSBisector(const MATedge & edge, PointSegmentBisector & psb) const
{
	// Get |p| and |q| the two elements whose bisector is the Voronoi
	// arc that we want to traverse
	const Site * p = &edge.site();
	const Site * q = &edge.twin()->site();
	Vec2d c0 = edge.from()->circumcircle.center_;
	Vec2d c1 = edge.to()->circumcircle.center_;
	if( ! p->is_point() ) {
		std::swap(p, q);
	}
	// Here, p is a point and q is a segment.
	Vec2d P, S, T;
	P = p->from;
	S = q->from;
	T = q->to;
	Vec2d segVec = T - S;
	Utils::reduceVec(segVec); // this modifies segVec
	const double h = Utils::distanceToLine(P, S, segVec);
	if( h == 0.0 ) {
		std::cerr << std::endl << "BAD QUADRATIC BEZIER ARC" << std::endl;
		std::cerr << "P: " << P << std::endl;
		std::cerr << "S--: " << S << std::endl;
		std::cerr << "--T: " << T << std::endl;
		std::cerr << "segVec: " << segVec << std::endl;
		throw;//exit(-1); // return d01;
	}
	double segReducedLen = segVec.length();
	double x0 = (c0 - P).dot(segVec) / segReducedLen;
	double x1 = (c1 - P).dot(segVec) / segReducedLen;
	if( x0 > x1 ) {
		segVec = -1.0 * segVec;
		x0 = - x0;
		x1 = - x1;
	}
	// fill the data
	psb.point = P;
	psb.h = h;
	psb.x0 = x0;
	psb.x1 = x1;
	psb.u = segVec / segReducedLen;
	Vec2d v(-psb.u.y(), psb.u.x());
	psb.leftToRight = v.dot(P - S) > 0.0;
	if( psb.leftToRight ) {
		psb.v = v;
	} else {
		psb.v = -1.0 * v;
	}
}



bool
MATGraph::
is_endpoint_of_segment(const Site & point_site, const Site & segment_site) {
	assert(point_site.is_point());
	assert(segment_site.is_segment());
	return point_site.point() == segment_site.from || point_site.point() == segment_site.to;
}

// This can be used, for example, to cut a zero-area hole along the medial axis, inside a polygon.
void
MATGraph::
linearApproxOfMedialAxis(std::vector<std::vector<Sample>> & out, bool collapsedOnly) const {
	out.clear();
	std::unordered_set<const MATedge *> visited;
	Sample s;
	s.tangent = Vec2d(0,0);
	if( collapsedOnly ) {
		for( const MATvert & vert : verts ) {
			if( vert.status != MatStatus::Collapsed ) {
				for( int i = 0; i < 3; ++i ) {
					visited.insert(& vert.edge[i]);
					visited.insert(vert.edge[i].twin());
				}
			}
		}
	}
	for( const MATvert & vert_ : verts ) {
		const MATvert * vert = &vert_;
		MATedge * freeway, * deadend;
		if( collapsedOnly && ( ! vert->is_collapsed() ) ) continue;
		if( (collapsedOnly && statusDegree(vert, MatStatus::Collapsed, freeway, deadend) == 0)
				||
			   (!collapsedOnly && degree(vert, true) == 0) ) {
			out.emplace_back();
			s.pos = vert->circumcircle.center_;
			s.radius = vert->circumcircle.radius_;
			out.back().push_back(s);
			out.back().push_back(s);
		}
		for( int i = 0; i < 3; ++i ) {
			if ( ! vert->has_neighbor(i) ) continue;
			if ( visited.count(& vert->edge[i]) > 0 ) continue;
			const MATedge * arc = & vert->edge[i];
			out.emplace_back();
			s.pos = arc->from()->circumcircle.center_;
			s.radius = arc->from()->circumcircle.radius_;
			out.back().push_back(s);
			while( true ) {
				visited.insert(arc);
				visited.insert(arc->twin());
				const MATvert * vert = arc->from();
				const MATvert * neighbor = arc->to();
				// OK, approx the arc
				Vec2d b = neighbor->circumcircle.center_;
				//out.back().push_back(a);
				if( arc->type != EdgeType::VertEdge ) {
					s.pos = b;
					s.radius = neighbor->circumcircle.radius_;
					out.back().push_back(s);
				} else {
					QuadraticArc qArc;
					buildBezierQuadraticOfArc(*arc, qArc);
					Vec2d focus;
					if( arc->site().is_point() )
						focus = arc->site().point();
					else
						focus = arc->twin()->site().point();
					diceQuadratic(focus, qArc, out.back());
				}
				bool at_least_one = false;
				for( int i = 0; i < 3; ++i ) {
					if ( ! neighbor->has_neighbor(i) ) continue;
					if ( visited.count(& neighbor->edge[i]) > 0 ) continue;
					at_least_one = true;
					arc = & neighbor->edge[i];
					break;
				}
				if( ! at_least_one ) break;
			}
		}
	}
}

void
MATGraph::
debugCheck() const
{
#ifndef NDEBUG
	for ( const MATvert & vert : verts ) {
		for (int i = 0; i < 3 ; ++i ) {
			assert(!vert.edge[i].to() || vert.edge[i].twin() != nullptr);
			assert(!vert.edge[i].to() || vert.edge[i].twin()->twin() == & vert.edge[i]);
			assert(!vert.edge[i].to() || vert.edge[i].from() == & vert);
			assert(!vert.edge[i].to() || vert.edge[i].next(true)->prev(true) == & vert.edge[i]);
			assert(!vert.edge[i].to() || vert.edge[i].from()->is_destroyed() || vert.edge[i].to()->is_destroyed() || vert.edge[i].next(false)->prev(false) == & vert.edge[i]);
		}
	}
#endif
}
