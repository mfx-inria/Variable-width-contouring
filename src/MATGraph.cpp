/* vim: set sts=0 ts=4 sw=4 noet : */

#include "MATGraph.h"

#include <iterator> // for std::distance
#include <algorithm> // rotate
#include <stack>
#include <map>

// - - - - - - - - - - - - - - - - - - - - - COLLAPSEDAXIS

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
#if ICESL_PLUGIN==0
					cerr << "[Variable-width contouring] Changing axis! from " << prevAP << " to " << ap << " (at "
						<< curAxis.clippingVertex->pos() <<  " --> "
						<< axis.clippingVertex->pos() <<  ")" << endl;
#endif
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
	//best.tangent = Vec2d(0,1);
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
			update(Sampling::moveTowardBisectorWithSegment(sample, points[ap.pos()], points[ap.pos()+1], bad), ap);
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

// - - - - - - - - - - - - - - - - - - - - - MATGRAPH

void
MATGraph::
print() const {
	std::map<const MATvert *, std::string> vName;
	std::vector<const MATvert *> sorted;
	for( const MATvert & v : verts ) {
		sorted.push_back(&v);
	}
	std::sort(sorted.begin(), sorted.end(), [](const MATvert * l, const MATvert * r) {
			return l->circumcircle.center_ < r->circumcircle.center_;
			});
	int i(0);
	for( const MATvert * v : sorted ) {
		int j = i++;
		vName[v] = std::string("");
		do {
			vName[v] = (char)('A' + (j%26)) + vName[v];
			j /= 26;
		} while( j > 0 );
	}
	for( const MATvert * v : sorted ) {
		cerr << vName[v] << " : " << v->circumcircle.center_ << ", " << v->status << endl;
		v->iter_edges([&](ConstEdgeIterator e) {
			const MATvert * n = e->to();
			if( e->type == EdgeType::VertEdge )
				cerr << "\t~~> ";
			else
				cerr << "\t--> ";
			cerr	<< vName[n] << '(' << n->circumcircle.center_ << ')' << endl;
		});
	}
}

void
MATGraph::
computeConnectedComponents(Components & components) {
	components.clear();
	VertexSet visited;
	for ( MATvert & vert_ : verts ) {
		MATvert * vert = & vert_;
		if( visited.count(vert) > 0 ) continue;
		components.emplace_back();
		VertexSet & component = components.back().vertexSet;
		std::vector< MATvert *> to_check;
		to_check.emplace_back(vert);
		component.emplace(vert);
		while ( ! to_check.empty() ) {
			MATvert * v = to_check.back(); to_check.pop_back();
			v->iter_edges([&](EdgeIterator edge) {
				MATvert * neighbor = edge->to();
				if( component.count(neighbor) == 0 ) {
					component.emplace(neighbor);
					to_check.emplace_back(neighbor);
				}
			});
		}
		visited.insert(component.begin(), component.end());
	}
}

int
MATGraph::
statusDegree(const MATvert * vert, const MatStatus status) const {
	/* Returns the number of neighbors with given status.
	 */
	EdgeIterator freeway;
	return statusDegree(vert, status, freeway);
}

int
MATGraph::
statusDegree(const MATvert *vert, const MatStatus status, EdgeIterator & freeway) const {
	/* Returns the number of neighbors with given status.
	   If it is 1, then freeway is the edge leading to the neighbor with that status.
	 */
	int deg = 0;
	MATvert *non_const_v = const_cast<MATvert *>(vert);
	freeway = non_const_v->edges_.end();
	non_const_v->iter_edges([&](EdgeIterator edge) {
		if( edge->to()->status == status ) {
			freeway = edge;
			++deg;
		}
	});
	if( deg != 1 ) freeway = non_const_v->edges_.end();
	return deg;
}

int
MATGraph::
degree(const MATvert * vert, bool walk_on_fire) const {
	/* Returns the number of MatStatus::Normal neighbors (or any status if |walk_on_fire| is true).
	 */
	EdgeIterator freeway;
	return degree(vert, freeway, walk_on_fire);
}

int
MATGraph::
degree(const MATvert * vert, EdgeIterator & freeway, bool walk_on_fire) const {
	/* Returns the number of MatStatus::Normal neighbors (or any status if |walk_on_fire| is true).
	   If it is 1, then freeway is the edge leading to the neighbor with that status.
	 */
	int deg = 0;
	MATvert *non_const_v = const_cast<MATvert *>(vert);
	freeway = non_const_v->edges_.end();
	non_const_v->iter_edges([&](EdgeIterator edge) {
		if( walk_on_fire || (edge->to()->status == MatStatus::Normal) ) {
			freeway = edge;
			++deg;
		}
	});
	if( deg != 1 ) freeway = non_const_v->edges_.end();
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
			for( ConstEdgeIterator edge = v->edges_.begin(); edge != v->edges_.end(); ++edge ) {
				MATvert * neighbor = edge->to();
				if( neighbor->is_collapsed() ) return true;
				if( visited.count(neighbor) > 0 ) continue;
				jobs.push(neighbor);
			}
		}
	}
	return false;
}

MATvert&
MATGraph::
insert(Component * component, EdgeIterator edge, const Disk & circumcircle)
{
	debugCheck();
	EdgeIterator original_twin = edge->twin();
	verts.emplace_back(circumcircle);
	MATvert & vert = verts.back();
	if( nullptr != component ) {
		component->vertexSet.insert(&vert);
	}
	vert.edges_.emplace_front();
	EdgeIterator e0 = vert.edges_.begin();
	vert.edges_.emplace_front();
	EdgeIterator e1 = vert.edges_.begin();
	//cerr << "e0, e1, begin " << &(*e0) << ' ' << &(*e1) << ' ' << &(*vert.edges_.begin()) << endl;
	e0->twin_ = edge;
	e0->to_ = original_twin->to_;
	e0->site_ = original_twin->site_;
	e0->type = original_twin->type;

	e1->twin_ = original_twin;
	e1->to_ = edge->to_;
	e1->site_ = edge->site_;
	e1->type = edge->type;
	assert(edge->type == original_twin->type);

	edge->twin_ = e0;
	edge->to_ = &vert;
	original_twin->twin_ = e1;
	original_twin->to_ = &vert;

	debugCheck();
	return vert;
}

void
MATGraph::
removeVertexSet(const std::unordered_set<MATvert *> to_remove) {
	for( auto vert_it = verts.begin(); vert_it != verts.end(); ) {
		MATvert * vert = &(*vert_it);
		if( to_remove.count(vert) > 0 ) {
			vert->iter_edges([](EdgeIterator edge) {
					// kill the twin
					edge->to()->edges_.erase(edge->twin());
			});
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
		if (vert_it->status != MatStatus::Normal) {
			vert_it->iter_edges([](EdgeIterator edge) {
					edge->to()->edges_.erase(edge->twin());
			});
			vert_it = verts.erase(vert_it);
		} else {
			++vert_it;
		}
	}
	debugCheck();
}

// - - - - - - - - - - - - - - - - - - - - - MATEDGE

EdgeIterator
cw(EdgeIterator e)
{
	MATvert * v = e->from();
	if( e == v->edges_.begin() )
		return --(v->edges_.end());
	return --e;
}

EdgeIterator
ccw(EdgeIterator e)
{
	MATvert * v = e->from();
	e++;
	if( e == v->edges_.end() )
		return v->edges_.begin();
	return e;
}

EdgeIterator
MATedge::
nextNotTrimmed() const {
	EdgeIterator edge = ccw(twin());
	while( edge->to()->is_trimmed() && (edge != twin()) )
		edge = ccw(edge);
	return edge;
}

EdgeIterator
MATedge::
nextNormalOrTwin() const {
	EdgeIterator edge = ccw(twin());
	while( edge->to()->is_destroyed() && (edge != twin()) )
		edge = ccw(edge);
	return edge;
}

EdgeIterator
MATedge::
nextNormalOrNextAny() const {
	EdgeIterator edge, e0;
	e0 = edge = ccw(twin());
	do {
		if( ! e0->to()->is_destroyed() ) return e0;
		e0 = ccw(e0);
	} while( e0 != edge );
	return edge;
}

EdgeIterator
MATedge::
next(bool walk_on_fire) const {
	assert(walk_on_fire || (from()->is_normal() && to()->is_normal()));
	EdgeIterator edge = ccw(twin());
	if( walk_on_fire ) return edge;
	while( edge->to()->is_destroyed() && (edge != twin()) )
		edge = ccw(edge);
	return edge;
}

EdgeIterator
MATedge::
prev(bool walk_on_fire) const {
	assert(walk_on_fire || (from()->is_normal() && to()->is_normal()));
	EdgeIterator self = twin()->twin();
	assert(this == &(*self));
	EdgeIterator edge = cw(self);
	if( walk_on_fire ) return edge->twin();
	while( edge->to()->is_destroyed() && (edge != self) )
		edge = cw(edge);
	return edge->twin();
}

Vec2d
MATedge::
getGeneratorLocation() const {
	if (site().is_point()) return site().point();
	assert(site().is_segment());
	assert(twin()->to());
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

// used to decide which color to draw some arcs
bool
MATGraph::
isShavedVertexConnectedToNormalAxis(const MATvert * v, const MATvert * from) const {
	bool r = false;
	for( ConstEdgeIterator edge = v->edges_.begin(); edge != v->edges_.end(); ++edge ) {
		const MATvert * neighbor = edge->to();
		if( neighbor->is_normal() ) return true;
		if( ! neighbor->is_shaved() ) continue;
		if( neighbor == from ) continue;
		r = r || isShavedVertexConnectedToNormalAxis(neighbor, v);
		if( r ) break;
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
		bool onlyTrimmed, bool noShave) const {
	cairo_set_line_width(context, lw);
	cairo_set_line_cap(context, CAIRO_LINE_CAP_BUTT);
	cairo_set_source_rgba(context,1,1,1,1);
	for( const MATvert & vert : verts ) {
		int nbArcs(0);
		for( ConstEdgeIterator edge = vert.edges_.begin(); edge != vert.edges_.end(); ++edge ) {
			++nbArcs;
			const MATvert & neighbor = *edge->to();
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
					if( isShavedVertexConnectedToNormalAxis(shavedVert) )
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
			if( edge->type != EdgeType::VertEdge ) {
				colorEdge(*edge);
				if( onlyTrimmed != trimmedEdge ) continue;
				cairo_move_to(context, a.x(), a.y());
				cairo_line_to(context, b.x(), b.y());
			} else {
				colorEdge(*edge);
				if( onlyTrimmed != trimmedEdge ) continue;
				QuadraticArc qArc;
				buildBezierQuadraticOfArc(*edge, qArc);
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
#endif // FILL_NO_CAIRO

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
				vert.iter_edges([&](ConstEdgeIterator edge) {
					visited.insert(&(*edge));
					visited.insert(&(*edge->twin()));
				});
			}
		}
	}
	for( const MATvert & vert_ : verts ) {
		const MATvert * vert = &vert_;
		if( collapsedOnly && ( ! vert->is_collapsed() ) ) continue;
		if( (collapsedOnly && statusDegree(vert, MatStatus::Collapsed) == 0)
				||
			   (!collapsedOnly && degree(vert, true) == 0) ) {
			out.emplace_back();
			s.pos = vert->circumcircle.center_;
			s.radius = vert->circumcircle.radius_;
			out.back().push_back(s);
			out.back().push_back(s);
		}
		for( ConstEdgeIterator arc = vert->edges_.begin(); arc !=  vert->edges_.end(); ++arc ) {
			if ( visited.count(&(*arc)) > 0 ) continue;
			out.emplace_back();
			s.pos = arc->from()->circumcircle.center_;
			s.radius = arc->from()->circumcircle.radius_;
			out.back().push_back(s);
			while( true ) {
				visited.insert(&(*arc));
				visited.insert(&(*arc->twin()));
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
				for( ConstEdgeIterator ne = neighbor->edges_.begin(); ne != neighbor->edges_.end(); ++ne ) {
					if ( visited.count(&(*ne)) > 0 ) continue;
					at_least_one = true;
					arc = ne;
					break;
				}
				if( ! at_least_one ) break;
			}
		}
	}
}

void
MATGraph::
splitApexes()
{
	struct Insertion {
		EdgeIterator edge;
		Vec2d mid;
		double radius;
		Insertion(const EdgeIterator e, const Vec2d & v, const double r) : edge(e), mid(v),
		radius(r) {}
	};
	list<Insertion> insertions;
	for ( MATvert & vert : verts ) {
		for( EdgeIterator edge_it = vert.edges_.begin(); edge_it != vert.edges_.end(); ++edge_it ) {
			MATedge & edge = *edge_it;
			if ( edge.to() < edge.from() ) continue; // arbitrary ordering to avoid unneeded apex twin ;)
			if ( edge.type == EdgeType::EdgeEdge ) continue; // there is no apex

			Vec2d a, b, base, x_dir;
			if ( edge.type == EdgeType::VertVert ) {
				a = edge.site().point();
				b = edge.twin()->site().point();
				base = b;
				x_dir = (b - a).rotatedCCW();
			} else {
				assert(edge.type == EdgeType::VertEdge);
				Site& point_site = edge.site().is_point() ? edge.site() : edge.twin()->site();
				Site& segment_site = edge.site().is_point() ? edge.twin()->site() : edge.site();
				a = point_site.point();
				Vec2d sa = segment_site.source();
				Vec2d sb = segment_site.destination();
				x_dir = sb - sa;
				b = Utils::project(a, sa, sb);
				base = sa;
			}
			Utils::reduceVec(x_dir); // reduce magnitude without losing precision (by divisions by powers of 2)
			double goes_left = x_dir.dot(vert.circumcircle.center_ - a);
			double goes_right = x_dir.dot(edge.to()->circumcircle.center_ - a);
			if( goes_left * goes_right < 0.0 ) { // otherwise apex is beyond the segment
				Vec2d y_dir = x_dir.rotatedCCW();
				y_dir.normalize();
				Vec2d mid = 0.5 * (a + b);
				double radius = y_dir.dot(a-base) * 0.5;
				insertions.emplace_back(edge_it, mid, std::fabs(radius));
			}
		}
	}
	for( auto & ins : insertions ) {
		insert(nullptr, ins.edge, Disk(ins.mid, ins.radius));
	}
}

void
MATGraph::
debugCheck() const
{
#if 0
//#ifndef NDEBUG
	for ( const MATvert & vert : verts ) {
		MATvert * v = const_cast<MATvert *>(&vert);
		for( EdgeIterator e = v->edges_.begin(); e != v->edges_.end(); ++e ) {
			assert(e->to());
			assert(e->twin() != e->to()->edges_.end());
			assert(e->twin()->twin() == e);
			assert(e->from() == v);
			assert(e->next(true)->prev(true) == e);
			assert(e->from()->is_destroyed() || e->to()->is_destroyed() || e->next(false)->prev(false) == e);
		}
	}
#endif
}
