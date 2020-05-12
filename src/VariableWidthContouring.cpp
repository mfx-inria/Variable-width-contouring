/* vim: set sts=0 ts=4 sw=4 noet : */

#include "VariableWidthContouring.h"
//#include "chronograph.h"

#include <cairo/cairo-pdf.h>
#include <map>
#include <stack>

#include <CGAL/squared_distance_2.h>

using namespace std;

BBox
VariableWidthContouring::
computeBBox() const {
	BBox box(graph.verts.front().pos());
	for ( MATvert & v : graph.verts ) {
		double r = v.radius();
		Vec2d c = v.pos();
		Vec2d rv(r,r);
		BBox sub(c-rv, c+rv);
		box.extendToCover(sub);
	}
	return box;
}

void
VariableWidthContouring::
gatherCollapsedSegment(const MATedge * edge, const MATedge * & limit, CollapsedAxis & axis) const {
	//cerr << "Entering axis at " << edge->from()->pos() << " --> " << edge->to()->pos() << endl;
	// First, skip the shaved vertices until the collapsed region
	const MATedge * e = edge;
	MATvert * from = e->from();
	MATvert * to = e->to();
	bool axisCreated(false);
#define MATGRAPH_ADVANCE e = e->nextNotTrimmed(); from = e->from(); to = e->to()
	while( from->is_shaved() ) {
		if( ! from->targets_.empty() ) {
			axis.addSubAxis(from->targets_[0].edgeSource->from());
			axisCreated = true;
		}
		MATGRAPH_ADVANCE;
	}
	if( ! axisCreated ) {
		if( ! from->targets_.empty() ) {
			axis.addSubAxis(from->targets_[0].edgeSource->from());
		} else { // the case of pure collapsed cycles
			axis.addSubAxis();
		}
	}
	if( from->is_normal() ) {
		limit = e;
		//cerr << "Exiting CLIPPED axis at " << e->from()->pos() << " --> " << e->to()->pos() << endl;
		return;
	}
	axis.addPoint(from->pos());
	bool addFrom = false;
	while( ! to->is_normal() ) {
		if( (! from->is_collapsed()) || (! to->is_collapsed()) ||
		  ( to->pos() == from->pos() ) ) {
			// yes, an MA edge can have zero length...
			goto advance;
		}
		if( addFrom ) {
			addFrom = false;
			axis.addPoint(from->pos());
		}
		if( e->type == EdgeType::VertEdge ) {
			QuadraticArc qArc;
			graph.buildBezierQuadraticOfArc(*e, qArc);
			Vec2d focus;
			if( e->site().is_point() )
				focus = e->site().point();
			else
				focus = e->twin()->site().point();
			std::vector<Sample> samples;
			diceQuadratic(focus, qArc, samples);
			for( const auto & s : samples ) axis.addPoint(s.pos);
		} else {
			axis.addPoint(to->pos());
		}
		if( ! to->targets_.empty() ) {
			assert( to->targets_[0].isClipping );
			axis.addSubAxis(to->targets_[0].edgeSource->from());
			addFrom = true;
		}
advance:
		MATGRAPH_ADVANCE;
		if( e == edge ) {
			limit = e;
			//cerr << "Exiting CYCLE axis at " << e->from()->pos() << " --> " << e->to()->pos() << endl;
			return;
		}
	}
	limit = e->next(true);
	//cerr << "Exiting NORMAL axis at " << limit->from()->pos() << " --> " << limit->to()->pos() << endl;
#undef MATGRAPH_ADVANCE
}

void
MAOffsetInfo::
swap(MAOffsetInfo & other) {
	lBoundary_.swap(other.lBoundary_);
	rBoundary_.swap(other.rBoundary_);
}

void
VariableWidthContouring::
setWorkingMAOffsetInfo(MAOffsetInfo * i) {
	maoi_ = i;
}

// Returns the distance from point |p| (assumed to lie on the |arc|) to the polygon features that
// define the arc.
double
VariableWidthContouring::
distanceToBoundary(const MATedge & arc, const Vec2d & p) const {
	const Site & s1 = arc.twin()->site();
	const Site & s0 = arc.site();
	if( s0.is_point() ) {
		return (s0.point() - p).length();
	} else if( s1.is_point() ) {
		return (s1.point() - p).length();
	} else {
		assert(s0.is_segment() && s1.is_segment());
		return Utils::distanceToLine(p, s0.source(), s0.segment_vector());
	}
}

/* Return a parameter along the edge (under standard parametrisation) that represents the position
 * of the result Disk.
 */
double
VariableWidthContouring::
findPointAtDistance(const MATedge& edge, double dist, Disk & disk) const
{
	assert(edge.to());
	double param = -1;
	double a_R = edge.from()->circumcircle.radius_;
	double b_R = edge.to()->circumcircle.radius_;
	if( dist <= a_R && dist <= b_R ) return param;
	if( dist >= a_R && dist >= b_R ) return param;
	const Vec2d a = edge.from()->circumcircle.center_;
	const Vec2d b = edge.to()->circumcircle.center_;
	const Vec2d ab = b - a;
	Vec2d loc;
	switch( edge.type ) {
		case EdgeType::EdgeEdge: {
			loc = a + ab * ((dist - a_R) / (b_R - a_R));
			param = (loc - a ).length();
			//cerr << "EE, param=" << param << endl;
		} break;
		case EdgeType::VertEdge: {
			 // TODO : make the code more robust to numerical errors
			const Site& point_site = edge.site().is_point() ? edge.site() : edge.twin()->site();
			const Site& segment_site = edge.site().is_point() ? edge.twin()->site() : edge.site();
			Vec2d sa = segment_site.source();
			Vec2d sb = segment_site.destination();
			Vec2d sab = sb - sa;
			Vec2d x_dir = sab;
			Utils::reduceVec(x_dir);
			const Vec2d s1 = point_site.point(); // point site
			x_dir.normalize();
			Vec2d y_dir = x_dir.rotatedCCW();
			if( y_dir.dot(s1-sa) > 0.0 ) y_dir.negate();
			const double d = y_dir.dot(sa-s1);
			const double currentOffset = (a-s1).length() - a_R;
			const double requiredRadius = dist + currentOffset;
			// x^2 + d^2 = 2*d*y, solve x for y=requiredRadius
			double x = std::sqrt((2 * requiredRadius - d) * d);
			double lo(x_dir.dot(a-s1)), hi(x_dir.dot(b-s1));
			if( lo > hi ) std::swap(lo, hi);
			if( x <= lo || x >= hi ) x = -x;
			loc = s1 + (x * x_dir) + ( d - requiredRadius ) * y_dir;
			param = x-lo;
			//cerr << "VE, param=" << param << endl;
		} break;
		case EdgeType::VertVert: default: {
			const Vec2d s1 = edge.site().point();
			const Vec2d s2 = edge.twin()->site().point();
			Vec2d s12 = s2 - s1;
			const Vec2d mid = ( s1 + s2 ) * 0.5;
			double s12Scale = Utils::reduceVec(s12);
			const double a_mid_length_2 = s12Scale * s12Scale * s12.squaredLength() * 0.25;
			const double currentOffset = (a-s1).length() - a_R;
			const double requiredRadius = dist + currentOffset;
			bool bad = requiredRadius * requiredRadius <= a_mid_length_2 - 1e-6;
			if( bad ) {
				cerr << "requiredRadius^2=" << (requiredRadius*requiredRadius) << ", a_mid_length_2=" << a_mid_length_2
					<< ", currentOffset=" << currentOffset
					<< endl << "s1=" << s1 << ", s2=" << s2 << endl
					<< ", a=" << a << ", a_R=" << a_R << endl
					<< ", b=" << b << ", b_R=" << b_R << ", dist=" << dist << ", (dist-a_R)=" <<
					(dist-a_R)
					<< endl;
			}
			assert( ! bad );
			double x = std::sqrt(max(0.0,requiredRadius * requiredRadius - a_mid_length_2));
			Vec2d x_dir = s12.normalized().rotatedCCW();
			double lo(x_dir.dot(a-s1)), hi(x_dir.dot(b-s1));
			if( lo > hi ) std::swap(lo, hi);
			if( x <= lo || x >= hi ) x = -x;
			loc = mid + (x * x_dir);
			param = x-lo;
			//cerr << "VV, param=" << param << endl;
		} break;
	}
	disk = Disk(loc, dist);
	//graph.insert(nullptr, edge, Disk(loc, dist));
	return param;
}


// Compute a description of the current shape boundary, made a alternate straight edge and tangent
// continuous circular arcs.
// It *ALSO* fills the lBoundary_ and rBoundary_ maps. So after each medial axis trimming, we must
// call this function.
void
VariableWidthContouring::
computeSmoothPaths(const Component & component, SmoothPaths & out, bool walk_on_fire, MAOffsetInfo * maoi) const {
	if( graph.verts.empty() ) {
		return;
	}
	if( nullptr == maoi ) {
		maoi = maoi_;
	}

	NewCycleFunction newCycle = [&](const MATedge * prev, const MATvert * vert, const MATedge * e) {
		out.emplace_back();
	};

	WalkFunction walkEdge = [&](const MATedge * edge, const MATvert * vert, const MATedge * nextEdge) {
		// |vert| is |edge->to()|
		const Site& curr_site = edge->site();
		const Site& next_site = nextEdge->site();
		const Disk & maximalDisk = vert->circumcircle;
		int d;
		d = degree(vert, walk_on_fire);

		//bool debug = walk_on_fire && vert->pos().x() > 16 && vert->pos().x() < 17
		//	&& vert->pos().y() > 16 && vert->pos().y() < 17;

		if( (0 == d) || (
				((curr_site.is_segment() == next_site.is_segment()) && ( curr_site != next_site ))
				|| ((curr_site.is_segment() && next_site.is_point())   && ( ! graph.is_endpoint_of_segment(next_site, curr_site)))
				|| ((next_site.is_segment() && curr_site.is_point())   && ( ! graph.is_endpoint_of_segment(curr_site, next_site)))
		  ) ) {
			double radius = maximalDisk.radius_;
			if( walk_on_fire && vert->is_normal() ) radius += component.uniformOffset;
			BoundaryCircle bc(ToTheRight, maximalDisk.center_, radius);
			//if( debug ) {
			//	cerr << "Vertex at " << vert->pos() << " with radius " << radius
			//		<< " and uniformOffset=" << component.uniformOffset << endl;
			//}
			if( 0 == d ) {
				// Since we don't know which adjacent edge will by processed, we put the boundary disk
				// on all three *outgoing* edges. No risk, since the vertex is disconnected (degree 0)
				for( int i = 0; i < 3; ++i ) {
					const MATedge * e = & vert->edge[i];
					maoi->rBoundary_[e] = Disk(maximalDisk.center_, radius);
				}
			} else {
				maoi->rBoundary_[nextEdge] = Disk(maximalDisk.center_, radius);
			}
			if ( ! out.back().empty() && out.back().back() == bc) {
				std::cerr << "Warning: trying to put the same circle in the path twice!\n";
			} else {
				out.back().push_back(bc);
			}
		}
		if( 0 != d && next_site.is_point() && (curr_site != next_site) ) {
			// |edge| is the *correct* arc on which maximalDisk.center_ lies.
			double r = distanceToBoundary(*edge, maximalDisk.center_) - maximalDisk.radius_;
			if( walk_on_fire && vert->is_normal() ) r -= component.uniformOffset;
			BoundaryCircle bc(ToTheLeft, next_site.point(), r);
			maoi->lBoundary_[nextEdge] = Disk(next_site.point(), r);
			assert(out.back().empty() || out.back().back() != bc);
			out.back().push_back(bc);
		}
		return nullptr;
	};

	ArcToInt visited;
	
	auto walkTheMedialAxis = [&](const MATedge * edge) {
		const MATedge * e = edge;
	};

	for( const MATvert * vert : component.vertexSet) {
		if( (!walk_on_fire) && vert->is_destroyed() ) continue;
		int deg = degree(vert, walk_on_fire);
		if( 0 == deg ) {
			newCycle(&vert->edge[0], vert, &vert->edge[0]);
			walkEdge(& vert->edge[0], vert, & vert->edge[0]);
		} else {
			for( int i(0); i < 3; ++i ) {
				if( ! vert->has_neighbor(i) ) continue;
				const MATedge * edge = & vert->edge[i];
				if( ( ! walk_on_fire ) && edge->to()->is_destroyed() ) continue;
				if( visited.find(edge) != visited.end() ) continue;
				newCycle(edge->prev(walk_on_fire), vert, edge);
				while( visited.find(edge) == visited.end() ) {
					visited[edge] = 1;
					const MATedge * nextEdge = edge->next(walk_on_fire);
					walkEdge(edge, edge->to(), nextEdge);
					edge = nextEdge;
				}
			}
		}
	}
}

int
VariableWidthContouring::
degree(const MATvert * vert, bool walk_on_fire) const {
	MATedge * freeway, * deadend;
	return graph.degree(vert, freeway, deadend, walk_on_fire);
}

int
VariableWidthContouring::
degree(const MATvert * vert, MATedge *& freeway, MATedge *& deadend, bool walk_on_fire) const {
	return graph.degree(vert, freeway, deadend, walk_on_fire);
}

int
VariableWidthContouring::
statusDegree(const MATvert * vert, MatStatus status, MATedge *& freeway, MATedge *& deadend) const {
	return graph.statusDegree(vert, status, freeway, deadend);
}

MATvert * closest(const Vec2d & query, const std::vector<MATvert *> & pts) {
	double dist = (pts[0]->pos()-query).squaredLength();
	MATvert * best = pts[0];
	for( int i(1); i < pts.size(); ++i ) {
		double d = (pts[i]->pos()-query).squaredLength();
		if( d < dist ) {
			best = pts[i];
			dist = d;
		}
	}
	return best;
}

void
VariableWidthContouring::
filterSmallFeatures(double min_radius, Components & components)
{
	//FIXME: this should use Clipping-shaving, rather that trimming. Otherwise we still have a
	//little bit of overlap. e.g. test-4.txt --scaling -.55 -m 0.3 -M 1.0
	for( Component & component : components ) {
		setCurrentComponent(component);
		vector<MATvert *> toAdd;
		std::vector<MATvert *> localMinimas;
		std::vector<MATvert *> ends;
		for( auto vert_ : component.vertexSet ) {
			MATvert & vert = *vert_;
			double vrad = vert.radius();
			bool local_minima = true;
			for( int i = 0 ; i < 3 ; ++i ) {
				if( ! vert.has_neighbor(i) ) continue;
				MATedge & edge = vert.edge[i];
				if( vrad > edge.to()->radius() ) local_minima = false;
				if( edge.to() < edge.from() ) continue; // arbitrary order to prevent double handling of half-edges
				Disk disk;
				if( findPointAtDistance(edge, min_radius, disk) > 0.0 ) {
					MATvert * v = & graph.insert(nullptr, edge, disk);
					ends.push_back(v);
					toAdd.push_back(v);
				}
			}
			if( local_minima && (vrad < min_radius) )
				localMinimas.push_back(&vert);
		}
		for( auto vert : component.vertexSet ) {
			if (vert->circumcircle.radius_ < min_radius)
				vert->trim();
		}
		for( MATvert * v : ends ) {
			MATvert * t = closest(v->pos(), localMinimas);
			double d = (v->pos()-t->pos()).length();
			if( d > min_radius ) continue;
			// setup fake target
			v->targets_.emplace_back(Disk(t->pos(),min_radius), min_radius);
		}
		for( auto v : toAdd ) component.vertexSet.insert(v);
		trimLeafs(true);
	}
	graph.removeDestroyed();
}
