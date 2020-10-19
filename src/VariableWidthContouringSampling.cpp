/* vim: set sts=0 ts=4 sw=4 noet tw=0 : */

#include "VariableWidthContouring.h"

using namespace std;

static bool debug = false;

Sample
getFirstSample(BoundaryCircles * maoi, const ConstEdgeIterator edge, const double minToolRadius) {

	ConstEdgeIterator e = edge;

	BitangentComputer bmake(minToolRadius);

	if( maoi->hasToTheRightBoundary(*edge) && ( ! edge->from()->fullyCleared_) ) {
		// The samples of a fully cleared circle all project to the same point (the center
		// of the circle) on the bisector, so we don't have to sample it.
		BoundaryCircle oc0, oc1, oc2;
		maoi->findBoundaryCirclesAdjacentToTTR(*edge, oc0, oc1, oc2, 0, true);
		Vec2d a,c,tgt;
		bmake(oc0, oc1, c, a);
		tgt = (a - oc1.center()).normalized();
		tgt.rotateCCW();
		return {a, tgt, 0};
	}

	while( ! maoi->hasBoundary(*e) ) e = e->prev(true);
	const BoundaryCircle bot = maoi->boundaryCircle(*e, ToTheLeft, 0);
	e = edge->next(true);
	while( ! maoi->hasBoundary(*e) ) e = e->next(true);
	const BoundaryCircle top = maoi->boundaryCircle(*e, ToTheRight, 0);

	Vec2d from(edge->from()->circumcircle.center_);

	if( edge->site().is_point() ) { // We are sampling a ToTheLeft circular arc
		Vec2d dir = (from - bot.center()).normalized();
		Vec2d pos = bot.center() + (bot.radius() * dir);
		dir.rotateCW();
		return {pos, dir, 0};
	}
	// From now on, we are sampling a straight edge
	Vec2d p0, p1, tangent;
	bmake(bot, top, p0, p1);
	if( p0 != p1 ) {
		tangent = (p1-p0).normalized();
	} else {
		tangent = (p0 - top.circle_.center_).normalized();
		if( top.passage_ == ToTheLeft )
			tangent.rotateCW();
		else
			tangent.rotateCCW();
	}
	p0 = p0 + (((from-p0) | tangent) * tangent);
	return {p0, tangent};
}

template< typename F >
void
sampleEdge(BoundaryCircles * maoi,
		const bool walk_on_fire,
		const ConstEdgeIterator edge,
		const F & out,
		SampleFunction & project,
		bool noProject,
		const double radiusOffset,
		const double minToolRadius) {

	BitangentComputer bmake(minToolRadius);

	//const MATvert * vert = edge.from();
	//bool debug = vert->pos().x() > 5.5 && vert->pos().x() < 6.0
	//				&& vert->pos().y() < -10.6 && vert->pos().y() > -11.0;

	if( maoi->hasToTheRightBoundary(*edge) && ( ! edge->from()->fullyCleared_) ) {
		// The samples of a fully cleared circle all project to the same point (the center
		// of the circle) on the bisector, so we don't have to sample it.
		BoundaryCircle oc0, oc1, oc2;
		maoi->findBoundaryCirclesAdjacentToTTR(*edge, oc0, oc1, oc2, radiusOffset, true);
		Vec2d a,b,c;
		bmake(oc0, oc1, c, a);
		bmake(oc1, oc2, b, c);
		a = a - oc1.center();
		b = b - oc1.center();
		Sampling::sampleCircularArc(a, b, oc1, out, /*clockwise*/false, project, radiusOffset);
	}

	ConstEdgeIterator e = edge;
	while( ! maoi->hasBoundary(*e) ) e = e->prev(walk_on_fire);
	const BoundaryCircle bot = maoi->boundaryCircle(*e, ToTheLeft, radiusOffset);
	e = edge->next(walk_on_fire);
	while( ! maoi->hasBoundary(*e) ) e = e->next(walk_on_fire);
	const BoundaryCircle top = maoi->boundaryCircle(*e, ToTheRight, radiusOffset);

	// |top| and |bot| are vertices on the outer boundary of the ribbon

	const Vec2d from(edge->from()->circumcircle.center_);
	e = edge->next(walk_on_fire);
	const Vec2d to(e->from()->circumcircle.center_);

	// |from| and |to| are vertices of the medial axis

	if( edge->site().is_point() ) { // We are sampling a ToTheLeft circular arc
		Sampling::sampleCircularArc(from-bot.center(), to-bot.center(), bot, out, true/*clockwise*/, project, radiusOffset);
		return;
	}
	// From now on, we are sampling a straight edge
	if( from == to ) return;
	Vec2d p0, p1, tangent;
	bmake(bot, top, p0, p1);
	if( p0 != p1 ) {
		tangent = (p1-p0).normalized();
	} else {
		tangent = (p0 - top.circle_.center_).normalized();
		if( top.passage_ == ToTheLeft )
			tangent.rotateCW();
		else
			tangent.rotateCCW();
	}
	p1 = p0 + (((to-p0) | tangent) * tangent);
	if( noProject ) {
		// We are sampling a straight line segment and NOT projecting it: a
		// single sample is enough.
		out({p1, tangent, radiusOffset});
		return;
	}
	// We are sampling a straight line segment AND projecting it on the bisector: we have to do adaptive sampling
	p0 = p0 + (((from-p0) | tangent) * tangent);
	struct BP {
		Sample border, projected; // |border| is on the border of the variable-width ribbon; |projected| is on the ribbon medial axis or bisector.
	};
	std::list<BP> sampling;
	auto insert = [&](typename std::list<BP>::iterator at, const Sample & s) {
		return sampling.insert(at, {s, project(s)});
	};
	Sample s{p0, tangent, 0};
	insert(sampling.end(), s);
	s = {0.5*(p0+p1), tangent, 0}; // insert middle for capturing details in collapsed-axis sampling
	insert(sampling.end(), s);
	s = {p1, tangent, 0};
	insert(sampling.end(), s);
	typename std::list<BP>::iterator it0, it1, last;
	it0 = sampling.begin();
	it1 = it0; ++it1;
	last = --sampling.end();
//#define FILL_SUBDIV_LIMIT 1000
#ifdef FILL_SUBDIV_LIMIT
	int N(0);
	while( (N < FILL_SUBDIV_LIMIT) && it1 != sampling.end() ) {
		++N;
#else
	while( it1 != sampling.end() ) { // simple adaptive sampling
#endif
		if( Sampling::is_good_sample(it0->projected, it1->projected) ) {
			out(it1->projected);
			it0 = it1++;
		} else { // splitting
			Vec2d delta = it1->border.pos - it0->border.pos;
#ifdef FILL_SUBDIV_LIMIT
			if( N == 1000) {
				cerr << "sampling straight edge: " << it1->border.pos << ", from=" << from << ", to=" << to
					<< ", radiusOffset=" << radiusOffset << endl;
				cerr << p0 << " " << p1 << " " << tangent << endl;
			}
#endif
#undef FILL_SUBDIV_LIMIT
			if( it1 != last ) {
				it1 = sampling.erase(it1);
				double k = 2.0/3.0;
				it1 = insert(it1, {it0->border.pos + (k * delta), it0->border.tangent, 0});
				k = 1.0/3.0;
				it1 = insert(it1, {it0->border.pos + (k * delta), it0->border.tangent, 0});
			} else {
				double k = 1.0/2.0;
				it1 = insert(it1, {it0->border.pos + (k * delta), it0->border.tangent, 0});
			}
		}
	}
}

extern bool gDebug;

void
VariableWidthContouring::
samplePrintPath(Component & component, BoundaryCircles * outerMaoi, vector<vector<Sample>> & pts) {

	BoundaryCircles * innerMaoi = maoi_; // for consistent reading with outerMaoi
	// We ASSUME the boundary circles information is already computed on both outerMaoi and innerMaoi

	const double constantOffset = component.uniformOffset / 2.0;

	BitangentComputer bmake(minToolRadius_);
	std::unordered_set<const MATedge *> visited;
	CollapsedAxis collapsedAxis;
	vector<Sample> * path;
	SmoothPaths & voids = component.voids;

	auto out = [&path](const Sample & val) {
		if( ! path->empty() ) {
			const Sample & b = path->back();
			if( b.pos == val.pos ) return;
			if( (b.pos-val.pos).l1norm() < Sampling::len_threshold/10.0 ) return;
		}
		path->push_back(val);
	};

	SampleFunction identity = [=](const Sample & s) { return s; };

	auto newCycle = [&]() {
		pts.emplace_back();
		path = & pts.back();
		voids.emplace_back();
	};

	auto lone_vertex_function = [&](const MATvert * vert) {
		// vert has full-degree 0
		BoundaryCircle oc;
		if( ! vert->edges_.empty() ) {
			oc = outerMaoi->boundaryCircle(vert->edges_.front(), ToTheRight);
		} else {
			oc.circle_ = vert->circumcircle;
			if( vert->is_normal() ) oc.radius() += component.uniformOffset;
		}
		Disk ic = vert->circumcircle;
		if( vert->is_collapsed() ) ic.radius_ = 0;
		double trackHalfWidth = 0.5*(oc.radius() - ic.radius_);
		oc.radius() = 0.5*(oc.radius() + ic.radius_);
		if( 0.0 == oc.radius() ) {
			cerr << "ERROR: we should not have an outer boundary circle of radius 0! at "
				<< vert->pos() << endl;
			return;
		}
		Sampling::sampleCircularArc({1,0}, {-1,0}, oc, out, false/*CCW*/, identity, trackHalfWidth);
		Sampling::sampleCircularArc({-1,0}, {1,0}, oc, out, false/*CCW*/, identity, trackHalfWidth);
	};

	auto limitAndProjectorForCollapsedPart = [&](MATvert * clipV, EdgeIterator edge, EdgeIterator & limit) {
		collapsedAxis.clear();
		gatherCollapsedSegment(edge, limit, collapsedAxis);
		if( nullptr == collapsedAxis.subAxes[0].clippingVertex ) {
			// Probably, the full collapsed stuff lies inside the clipping vertex...
			collapsedAxis.subAxes[0].clippingVertex = clipV;
		}
	};

	auto edge_function = [&](EdgeIterator edge, bool makeNewCycle) {
		enum ProjectionType {
			Uniform, ToDisk, Collapse
		};
		while( true ) {
			if( visited.end() != visited.find(&(*edge)) ) break; // edge already processed
			const MATvert * startVertex = edge->from();
			const MATvert *   endVertex = edge->to();
			EdgeIterator limit = edge->next(true);
			SampleFunction project = identity;
			ProjectionType projectionType;
			bool edge_equals_limit = false;
			BoundaryCircle ic;
			double constantOffsetForCurrentEdge = constantOffset;
			//cerr << startVertex->circumcircle.center_ << '(' << startVertex->status << ") -> ";
			//cerr << endVertex->circumcircle.center_ << '(' << endVertex->status << ')' << endl;

			if( startVertex->is_normal() ) { // START IS A NORMAL VERTEX
				if( endVertex->is_normal() ) {
					limit = edge->next(true);
					projectionType = Uniform;
				} else if( endVertex->is_trimmed() ) {
					// compute limit (arrival vertex post-projection of trimmed part) and the edge
					// that has the ToTheRight boundary circle:
					limit = ccw(edge);
					EdgeIterator toTheRightBearer = limit;
					if( 0 < degree(startVertex) ) { // this guarantees that the while loop terminates
						while( ! toTheRightBearer->to()->is_normal() )
							toTheRightBearer = ccw(toTheRightBearer);
					}
					ic = innerMaoi->boundaryCircle(*toTheRightBearer, ToTheRight, 0.0);
					project = [=](const Sample & s) {
						return Sampling::moveTowardBisectorWithDisk(s, constantOffsetForCurrentEdge, ic.circle_);
					};
					projectionType = ToDisk;
				} else if( endVertex->is_collapsed() || endVertex->is_shaved() ) {
					// the call below also sets the EdgeIterator |limit|
					limitAndProjectorForCollapsedPart(edge->from(), edge->nextNotTrimmed(), limit);
					constantOffsetForCurrentEdge = 0.0;
					project = [&](const Sample & s) { return collapsedAxis(s); };
					projectionType = Collapse;
				}
			} else if( startVertex->is_collapsed() ) {
				int collapsedDegree = graph.statusDegree(startVertex, MatStatus::Collapsed);
				if( 1 <= collapsedDegree ) {
					if( ! startVertex->targets_.empty() ) return; // will be handle by a Normal vertex
					if( ! endVertex->is_collapsed() ) return; // will be handled from a different edge
					// the call below also sets the EdgeIterator |limit|
					limitAndProjectorForCollapsedPart(nullptr, edge, limit);
					constantOffsetForCurrentEdge = 0.0;
					project = [&](const Sample & s) { return collapsedAxis(s); };
					projectionType = Collapse;
					Sample s = project(getFirstSample(outerMaoi, edge, minToolRadius_));
					collapsedAxis.reorder();
				} else if( 0 == collapsedDegree ) { // START IS A COLLAPSED VERTEX ADJACENT TO ONLY SHAVED VERTICES (I THINK!)
					limit = ccw(edge);
					Disk disk(startVertex->pos(), 0);
					constantOffsetForCurrentEdge = 0.0;
					project = [=](const Sample & s) {
						return Sampling::moveTowardBisectorWithDisk(s, constantOffsetForCurrentEdge, disk);
					};
					projectionType = ToDisk;
				} else {
					return;
				}
			} else {
				cerr << "ERROR: startVertex should never be Shaved or Trimmed, but it is " << startVertex->status << endl;
				graph.print();
				exit(-1);
			}

			edge_equals_limit = ( edge == limit );

			if( makeNewCycle ) {
				makeNewCycle = false;
				newCycle();
			}

			size_t firstSampleIndex = path->size();
			const EdgeIterator startEdge = edge;
			while( (edge != limit) || edge_equals_limit ) {
				if( visited.end() != visited.find(&(*edge)) ) break;
				if( innerMaoi->hasToTheRightBoundary(*edge) ) {
					BoundaryCircle bc = innerMaoi->boundaryCircle(*edge, ToTheRight);
					if( voids.back().empty() || (voids.back().back() != bc) ) {
						//cerr << "TTR at " << voids.back().back().circle_.center_ << " [R= "
						//	<< bc.circle_.radius_ << ']' << (voids.back().empty() ? "EMPTY" : "NOT EMPTY") << endl;
						voids.back().push_back(bc);
					}
				}
				if( innerMaoi->hasToTheLeftBoundary(*edge) ) {
					voids.back().push_back(innerMaoi->boundaryCircle(*edge, ToTheLeft));
					//cerr << "TTL at " << voids.back().back().circle_.center_ << endl;
				}
				edge_equals_limit = false;
				sampleEdge(outerMaoi, true, edge, out,
						project, Uniform == projectionType,
						constantOffsetForCurrentEdge, minToolRadius_);
				visited.insert(&(*edge));
				edge = edge->next(/*walk_on_fire*/true);
				collapsedAxis.jump();
			}
			if( Collapse == projectionType  ) {
				bool addClip = startEdge->from()->is_normal() && degree(startEdge->from()) > 0;
				if( addClip ) {
					BoundaryCircle bc(ToTheRight, startEdge->from()->pos(), startEdge->from()->radius());
					if( voids.back().empty() || (voids.back().back() != bc) ) {
						//cerr << "ENDING TTR at " << startEdge->from()->pos() << " [R= "
						//		<< startEdge->from()->radius() << ']' << endl;
						voids.back().emplace_back(bc);
					}
				}
				collapsedAxis.computeBoundaryCircles(voids.back(), path->begin() + firstSampleIndex, path->end());
				collapsedAxis.clear();
			}
		}
		if( (voids.back().size() > 1) && (voids.back().back() == voids.back().front()) ) {
			// It might happen that the clipping disks that start and end the
			// Collapsed Axis are the same disks AND that this disk IS the full
			// Normal shape that remains... Then the computation of the
			// boundary will find the same disk repeated and raise an
			// exception. So we have to delete one.
			// Example: test-6.tt --scaling 3 -m 0.3 -M 1, step 11 : the small
			// collapsed axis in the low middle.
			voids.back().pop_back();
		}
	};

	/* The limits of a collapse unit are either
		 - a degree one (or zero) Collapsed vertex (except if it has a clipping disk)
		 - a degree two Collapsed vertex in case of cycles
		 - a Normal vertex, holding a clipping disk.
	   We gather collapsed segment and clipping disk along the way.
	   Trimmed edges are easier to handle: they go from a Normal vertex to the same Normal vertex,
	   and we project onto a single ToTheRight disk.
	 */

	// First pass to sample all cycles that have Normal vertex or is just an isolated vertex.
	for( MATvert * vert : component.vertexSet ) {
		if( 0 < degree(vert, true) ) {
			if( ! vert->is_normal() ) continue;
			for( EdgeIterator e = vert->edges_.begin(); e != vert->edges_.end(); ++e ) {
				// we start a sampling run at a Normal vertex:
				if( visited.end() != visited.find(&(*e)) ) continue; // edge already processed
				edge_function(e, true);
			}
		} else {
			newCycle();
			lone_vertex_function(vert);
		}
	}
	// Second pass for pure collapsed cycles
	for( MATvert * vert : component.vertexSet ) {
		if( ! vert->is_collapsed() ) continue;
		if( 0 == degree(vert, true) ) continue;
		for( EdgeIterator e = vert->edges_.begin(); e != vert->edges_.end(); ++e ) {
			// we start a sampling run at a Normal or Collapsed vertex:
			if( visited.end() != visited.find(&(*e)) ) continue; // edge already processed
			edge_function(e, true);
		}
	}

	// Remove some duplicates
	for( vector<Sample> & samples : pts ) {
		if( samples.empty() ) continue;
		const Sample & b = samples.back();
		const Sample & f = samples.front();
		if( (b.pos-f.pos).l1norm() < Sampling::len_threshold/10.0 ) {
			samples.pop_back();
		}
	}
}

