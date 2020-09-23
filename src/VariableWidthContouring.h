/* vim: set sts=0 ts=4 sw=4 noet : */

#ifndef VARIABLE_WIDTH_DENSE_INFILLER_H_
#define VARIABLE_WIDTH_DENSE_INFILLER_H_

#include <cassert>
#include "utils.h"

#include <unordered_map>
#include <unordered_set>
#include <functional> // std::function

#include "MATGraph.h"
#include "Sampling.h"

typedef struct _cairo cairo_t;

// A struct that stores the circles that support part of the boundary of a round shape.
// Each arc of the medial axis can have its "right boundary" (see paper for definition),
// supported by at most one TTR (=To The Right, i.e. convex) circle and one TTL (concave) circle.
struct BoundaryCircles {
	typedef std::unordered_map<const MATedge*, Disk> ArcToDisk;

	// When walking around the medial axis, we think of an arc/edge as the relative interior of the
	// edge and the end-point, but we don't include the start-point. In this context:
	// lBoundary_ associates the relative interior  of an MATedge, to a ToTheLeft  boundary circle, should it exist.
	// rBoundary_ associates the end-point of an MATedge, to a ToTheRight boundary circle, should it exist.
	// (One MATedge could have none, one, the other or both.)
	ArcToDisk lBoundary_; // ToTheLeft boundary disks, concave path boundary
	ArcToDisk rBoundary_; // ToTheRight boundary disks, convex path boundary

	public:

	BoundaryCircles() : lBoundary_(), rBoundary_() {}

	void clear() {
		lBoundary_.clear();
		rBoundary_.clear();
	}

	void swap(BoundaryCircles & other);

	bool hasToTheLeftBoundary(const MATedge & e) const {
		return lBoundary_.find(&e) != lBoundary_.end();
	}
	const Disk & leftBoundary(const MATedge & e) const {
		return lBoundary_.at(&e);
	}
	bool hasToTheRightBoundary(const MATedge & e) const {
		return rBoundary_.find(&e) != rBoundary_.end();
	}
	const Disk & rightBoundary(const MATedge & e) const {
		return rBoundary_.at(&e);
	}
	bool hasBoundary(const MATedge & e) const {
		return hasToTheLeftBoundary(e) || hasToTheRightBoundary(e);
	}
	BoundaryCircle boundaryCircle(const MATedge & e, Passage preferredPassage, double rd = 0.0) const {
		bool hasRight = hasToTheRightBoundary(e);
		bool hasLeft = hasToTheLeftBoundary(e);
		BoundaryCircle bc;
		if( hasLeft && ( (! hasRight) || (ToTheLeft == preferredPassage) ) ) {
			const Disk & d = lBoundary_.at(&e);
			bc = BoundaryCircle(ToTheLeft, d.center_, d.radius_ + rd);
		} else {
			const Disk & d = rBoundary_.at(&e);
			bc = BoundaryCircle(ToTheRight, d.center_, d.radius_ - rd);
		}
		return bc;
	}
	// TTR = ToTheRight
	void findBoundaryCirclesAdjacentToTTR(const MATedge & edge, BoundaryCircle & pbc,
			BoundaryCircle & curbc, BoundaryCircle & nbc, double delta=0.0, bool walk_on_fire=false)
	{
		// ASSUMES |edge| has a ToTheRight boundary circle
		curbc = boundaryCircle(edge, ToTheRight, delta);
		EdgeIterator e = edge.prev(walk_on_fire);
		while( ! hasBoundary(*e) ) e = e->prev(walk_on_fire);
		pbc = boundaryCircle(*e, ToTheLeft, delta);
		if( ! hasToTheLeftBoundary(edge) ) {
			e = edge.next(walk_on_fire);
			while( ! hasBoundary(*e) ) e = e->next(walk_on_fire);
			nbc = boundaryCircle(*e, ToTheRight, delta);
		} else {
			nbc = boundaryCircle(edge, ToTheLeft, delta);
		}
	}
};

class VariableWidthContouring
{
public:
	typedef std::unordered_map<const MATedge*, int> ArcToInt;

	// A target. In this case, the |center_| member variable is useless.
	typedef std::vector<Target> Targets;

	using VertexSet  = MATGraph::VertexSet;
	using Component  = MATGraph::Component;
	using Components = MATGraph::Components;
	
	MATGraph& graph;

public:

	using NewCycleFunction = const std::function<void (ConstEdgeIterator, const MATvert *, ConstEdgeIterator)>;
	using WalkFunction     = const std::function<void (ConstEdgeIterator, const MATvert *, ConstEdgeIterator)>;
	using SampleFunction = std::function<Sample (const Sample &)>;
	
	VariableWidthContouring(MATGraph& graph, double minToolRadius, double maxToolRadius, double
			simplifyThreshold)
	: graph(graph)
	, maoi_(nullptr)
	, minToolRadius_(minToolRadius)
	, maxToolRadius_(maxToolRadius)
	, currentComponent_(nullptr)
	, simplifyThreshold_(simplifyThreshold)
	, maxCollapseRadius_(2 * maxToolRadius)
	{
	}
	
	double minCollapseRadius() { return minToolRadius_ * 4; }
	double maxCollapseRadius() { return maxToolRadius_ * 2; }

	// Returns the distance from point |p| (assumed to lie on the |arc|) to the polygon features that
	// define the arc.
	double distanceToBoundary(const MATedge & arc, const Vec2d & p) const;
	
	/*!
	 * change the graph and add a point at a given distance
	 * \return error code: 0 = success, -1 = dist is too small, 1 = dist is too large
	 */
	double findPointAtDistance(const MATedge & edge, double dist, Disk & disk) const;
	
	// Does the circumcenter of the face |f| happen to be a vertex of the input polygon?
	bool isInputVertex(MATvert * f) const;

	void gatherCollapsedSegment(EdgeIterator startCollapsedEdge, EdgeIterator & endNormalEdge,
			CollapsedAxis & axis) const;

	BBox computeBBox() const;
	void computeSmoothPaths(const Component &, SmoothPaths & out, bool walk_on_fire = false, BoundaryCircles * maoi = nullptr) const;
	void computeMedialAxis(const Paths & slice);
	SampleFunction makeProjectorOnMA(const MATedge & arc) const;

	/*
	 * Trim an arc to clear any target
	 * 
	 * \param minRadius min tool radius
	 * \param arc edge to trim
	 * \param targets targets to clear
	 * \param[out] trim main output
	 * \param trimmed whether we succeeded
	 * \return the single leaf for which the target circle was totally cleared (if any)
	 */
	MATvert * trimArcToTarget(const MATedge & arc, const Targets & targets, Disk & trim, bool & trimmed) const;

	const Target * trimArcToTargetForCollapse(const MATedge & arc, const Targets & targets, Disk & trim, bool & trimmed) const;
	
	void shrinkMA(Components &, std::vector<Disk> & cheekPrecursors, bool sharpCut);
	void shrinkMA_trim(Components &);
	bool shrinkMA_collapse(Components &, std::vector<Disk> & cheekPrecursors, bool sharpCut);
	void shrinkMA_cheeks(Components &);
	void shrinkMA_globalOffset(Components &);
	void clipAndShaveCollapsedParts(Components &, bool noShaving);
	void initializeTargets(bool walk_on_fire = false);
	bool collapse(std::vector<Disk> & cheekPrecursors, bool sharpCut);
	
	/*!
	 * \return The largest clearedDistance
	 */
	void clip(EdgeIterator start, const Target &, vector<MATvert *> &);
	void shaveLeafs(std::list<EdgeIterator> & propagations);
	double trimLeafs(bool trimVertices);
	void reduceRadii();

	void simplifyMACombinatorics(Components * components = nullptr);
	bool simplifyCollapsedMAGeometry(Components * components);
	bool simplifyAllMAGeometry(Components * components = nullptr);

protected:
	bool simplifyMAGeometry(Components * components, bool onlyCollapsed);

public:
	
	void filterSmallFeatures(double min_radius, Components &);

	/*!
	  * Mark all branches as trimmed. Keeps only cycles, which are used as print paths.
	  */
	void removeAllBranches();

	/* Return the degree of the vertex in the medial axis (represented by the BoundaryCircles passed as
	 * last parameter). There are 2 inout parameters:
	 	- freeway is the index of the unique neighbor alive when the degree is 1.
		- deadend is the index of the unique dead neighbor when the degree is 2.
		- In all other cases, freeway and/or deadend is set to -1.
	*/
	int statusDegree(const MATvert * vert, MatStatus, EdgeIterator & freeway) const;
	int degree(const MATvert * vert, EdgeIterator & freeway, bool walk_on_fire = false) const;
	int degree(const MATvert * vert, bool walk_on_fire = false) const;
	void analyse(double & minRadius, double & minInside, double & maxRadius) const;
	void setWorkingBoundaryCircles(BoundaryCircles *);
	void samplePrintPath(Component &, BoundaryCircles * prevInfo, std::vector<std::vector<Sample>> & pts);

	void setCurrentComponent(Component & c) { currentComponent_ = & c; }

	void clear() {
		maoi_->clear();
	}

protected:
	double maxToolRadius_, minToolRadius_, simplifyThreshold_;
	BoundaryCircles * maoi_;
	Component * currentComponent_;
	double maxCollapseRadius_;
};

template< typename O >
O & operator<<(O & o, const MATedge & arc) {
	o << &(*arc.from()) << ':' << arc.to();
	return o;
}

#endif // VARIABLE_WIDTH_DENSE_INFILLER_H_
