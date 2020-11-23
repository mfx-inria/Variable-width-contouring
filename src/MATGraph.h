/* vim: set sts=0 ts=4 sw=4 noet : */

#ifndef _FILL_MAT_GRAPH_H_
#define _FILL_MAT_GRAPH_H_

#include "fill_config.h"

#ifndef FILL_NO_CAIRO
#include <cairo/cairo-pdf.h>
typedef struct _cairo cairo_t;
#endif

#include "vec.h"
#include "Sampling.h"

#include <list>
#include <unordered_set>

enum class EdgeType
{
	EdgeEdge,
	VertEdge,
	VertVert
};

enum class MatStatus
{
	Normal,
	Trimmed,
	Collapsed,
	Shaved
};

template<class O>
O & operator<<(O & o, const MatStatus & s) {
	switch( s ) {
		case MatStatus::Normal : o << "Normal"; break;
		case MatStatus::Trimmed : o << "Trimmed"; break;
		case MatStatus::Collapsed : o << "Collapsed"; break;
		case MatStatus::Shaved : o << "Shaved"; break;
	}
	return o;
}

struct Site
{
	Vec2d from;
	Vec2d to; // uninitialized if not |is_point|
	bool is_point() const { return is_point_; }
	bool is_segment() const { return !is_point(); }
	Vec2d point() const { assert(is_point()); return from; }
	Vec2d source() const { assert(is_segment()); return from; }
	Vec2d destination() const { assert(is_segment()); return to; }
	Vec2d segment_vector() const { assert(is_segment()); return to - from; }
	Site() {} // uninitialized!!
	Site(Vec2d from)
	: from(from)
	, is_point_(true)
	{}
	Site(Vec2d from, Vec2d to)
	: from(from)
	, to(to)
	, is_point_(false)
	{}
	bool operator==(const Site & rhs) const {
		if (is_point() != rhs.is_point()) return false;
		if (is_point())
			return from == rhs.from;
		else
			return (from == rhs.from && to == rhs.to) || (from == rhs.to && to == rhs.from);
	}
	bool operator!=(const Site & rhs) const { return !(*this == rhs); }
protected:
	bool is_point_;
};

template<typename S> S & operator<<(S & s, const Site & site) {
	if (site.is_point())
		s << site.from;
	else
		s << site.from << " _ " << site.to;
}

class MATvert; // forward declare
class MATedge; // forward declare
using EdgeList = std::list<MATedge>;
using EdgeIterator = EdgeList::iterator;
using ConstEdgeIterator = EdgeList::const_iterator;

// Note the MATvert * |vertexSource|. It stores the MA vertex from which the target "comes from."
struct Target {
	Disk medialAxisDisk;
	double radius;

	bool isClipping;
	MATvert * clipDirection;

	union {
		MATvert * vertexSource;
		EdgeIterator edgeSource;
	};

	Target(const Disk & d) : medialAxisDisk(d), radius(0.0),
    vertexSource(nullptr), isClipping(false) {}
	Target(const Disk & d, double tr) : medialAxisDisk(d), radius(tr),
    vertexSource(nullptr), isClipping(false) {}
  Target(const Target& t) : medialAxisDisk(t.medialAxisDisk), radius(t.radius),
    isClipping(t.isClipping), clipDirection(t.clipDirection),
    vertexSource(t.vertexSource), edgeSource(t.edgeSource) {}
  Target& operator =(const Target& t)
  {
    if (this != &t) {
      medialAxisDisk = t.medialAxisDisk;
      radius = t.radius;
      isClipping = t.isClipping;
      clipDirection = t.clipDirection;
      vertexSource = t.vertexSource;
      edgeSource = t.edgeSource;
    }
    return *this;
  }
  ~Target() {}
};

class MATedge
{
public:
	MATvert * to_;
	Site site_;
	EdgeIterator twin_;
	EdgeType type;

	MATedge() {}
	Site& site() { return site_; }
	const Site& site() const { return site_; }
	MATvert * to() const { return to_; };
	MATvert * from() const { return twin()->to(); };
	EdgeIterator twin() const { return twin_; }
	EdgeIterator next(bool walk_on_fire = false) const;
	EdgeIterator nextNormalOrTwin() const;
	EdgeIterator nextNormalOrNextAny() const;
	EdgeIterator nextNotTrimmed() const;
	EdgeIterator prev(bool walk_on_fire = false) const;
	bool operator==( const MATedge & rhs ) const { return this == &rhs; }
	bool operator!=( const MATedge & rhs ) const { return !(*this == rhs); }

	MATedge & operator=(const MATedge & rhs) =delete;
	MATedge(const MATedge & rhs) =delete;

	Vec2d getGeneratorLocation() const;
	double getRadiusGradientAtStart() const;
};

EdgeIterator ccw(EdgeIterator e);
EdgeIterator cw(EdgeIterator e);

class MATvert
{
public:
	using EdgeList = std::list<MATedge>;
	using EdgeIterator = EdgeList::iterator;
	EdgeList edges_;
	Disk circumcircle;
	// algorithm specific members
	std::vector<Target> targets_;
	MatStatus status = MatStatus::Normal;

	// fullyCleared_ is True if and only if the boundary circle does NOT need to be sampled when we
	// project it onto the print path. We don't need to sample it because this piece of circular
	// boundary is guaranteed to project onto a single point (the circle center), which will be
	// sampled anyway by the previous segment. This is a particular case which causes a non-C1 point
	// in the print path (tangent discontinuity).
	// The boolean is used only when sampling print-paths (sampling.cpp)
	// The boolean is true when the disk target is the first to be cleared and its radius is <=
	// maxToolRadius (So the track *inner* contour is tangent to the local osculating circle.)
	bool fullyCleared_ = false;

	const Vec2d & pos() const { return circumcircle.center_; }
	double radius() const { return circumcircle.radius_; }

	template<typename F>
		void iter_edges(const F & f) const {
			for( ConstEdgeIterator it = edges_.begin(); it != edges_.end(); ++it ) {
				f(it);
			}
		}
	template<typename F>
		void iter_edges(const F & f) {
			for( EdgeIterator it = edges_.begin(); it != edges_.end(); ++it ) {
				f(it);
			}
		}

	bool is_destroyed() const { return status != MatStatus::Normal; }
	bool is_normal()    const { return status == MatStatus::Normal; }
	bool is_trimmed()   const { return status == MatStatus::Trimmed; }
	bool is_collapsed() const { return status == MatStatus::Collapsed; }
	bool is_shaved()    const { return status == MatStatus::Shaved; }
	void labelAsNormal()
	{
		status = MatStatus::Normal;
	}
	void trim()
	{
		status = MatStatus::Trimmed;
	}
	void collapse()
	{
		status = MatStatus::Collapsed;
	}
	void shave()
	{
		status = MatStatus::Shaved;
	}
	MATvert(Disk circumcircle)
	: circumcircle(circumcircle)
	{
	}
	MATvert & operator=(const MATvert & rhs) =delete;
	MATvert(const MATvert & rhs) =delete;
};


inline bool operator>(const MATvert & lhs, const MATvert & rhs) {
	return lhs.circumcircle.center_.x() > rhs.circumcircle.center_.x() ||
           ( lhs.circumcircle.center_.x() == rhs.circumcircle.center_.x() && lhs.circumcircle.center_.y() > rhs.circumcircle.center_.y() );
}
inline bool operator<(const MATvert & lhs, const MATvert & rhs) {
	return lhs.circumcircle.center_.x() < rhs.circumcircle.center_.x() ||
           ( lhs.circumcircle.center_.x() == rhs.circumcircle.center_.x() && lhs.circumcircle.center_.y() < rhs.circumcircle.center_.y() );
}

class MATGraph
{
public:

	using VertexSet  = std::unordered_set<MATvert *>;
	struct Component {
		VertexSet vertexSet;
		SmoothPaths innerSmoothPaths, outerSmoothPaths, simpleOuterSmoothPaths, voids;
		double uniformOffset;
		Component() : vertexSet(), innerSmoothPaths(), outerSmoothPaths(), uniformOffset(0) {}
	};
	using Components = std::list<Component>;

	void computeConnectedComponents(Components & components);

	std::list<MATvert> verts;

	MATvert& insert(Component * component, EdgeIterator edge, const Disk & circumcircle);
	void removeDestroyed();

	/* Helper function that actually removes the MATvert's from the std::list<MATvert>
	 */
	void removeVertexSet(const std::unordered_set<MATvert *> to_remove);

	bool componentHasCollapsedPart(const VertexSet &) const;

	/* Return the degree of the vertex in the medial axis (represented by the BoundaryCircles passed as
	 * last parameter). There are 2 inout parameters:
		- freeway is the index of the unique neighbor alive when the degree is 1.
		- In all other cases, freeway to -1.
	*/
	int degree(const MATvert * vert, EdgeIterator & freeway, bool walk_on_fire = false) const;
	int degree(const MATvert * vert, bool walk_on_fire = false) const;
	int statusDegree(const MATvert * vert, const MatStatus, EdgeIterator & freeway) const;
	int statusDegree(const MATvert * vert, const MatStatus) const;

	// helpers:
#ifndef FILL_NO_CAIRO
	bool isShavedVertexConnectedToNormalAxis(const MATvert * v, const MATvert * from=nullptr) const;
	void drawMedialAxis(cairo_t * context, double lw, bool drawVertices, bool onlyTrimmed, bool noShave=false) const;
	void helper_quadratic_to(cairo_t *cr, double x1, double y1, double x2, double y2) const;
#endif

	static bool is_endpoint_of_segment(const Site & point_site, const Site & segment_site);
	void buildBezierQuadraticOfArc(const MATedge & edge, QuadraticArc & qArc) const;
	void makePSBisector(const MATedge & edge, PointSegmentBisector & psb) const;

	void linearApproxOfMedialAxis(std::vector<std::vector<Sample>> &, bool collapsedOnly = false) const;

	void splitApexes();

	void debugCheck() const;
	void print() const;
};

// - - - - - - - - - - - - - - - - - - - - COLLAPSED AXIS

// A mixture of MA (linearized) segments and maximal disks centerered at Normal vertices
struct CollapsedAxis {
	using Samples = std::vector<Sample>;
	using VertexSet = MATGraph::VertexSet;
	struct SubAxis { // A sequence of points that start with a clippingVertex
		MATvert * clippingVertex; // the vertex' circumcircle is the actual clipping disk.
		std::vector<Vec2d> points; // 2 consecutive points make a MA segment.
	};
	std::vector<SubAxis> subAxes;
	AxisPos curAP;
	AxisPos farAP;

	void clear();
	void jump();
	void reorder();
	void print() const;
	void addSubAxis(MATvert * v = nullptr);
	void addPoint(const Vec2d & p);
	// Projects the sample toward the bisector with the collapsed axis
	Sample operator()(const Sample & sample);
	void computeBoundaryCircles(SmoothPath &, Samples::const_iterator, Samples::const_iterator)
		const;
	bool incr(AxisPos &) const;
};

#endif // _FILL_MAT_GRAPH_H_
