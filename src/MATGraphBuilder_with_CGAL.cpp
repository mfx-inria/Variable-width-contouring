/* vim: set sts=0 ts=4 sw=4 noet : */

// For Delaunay/Voronoi
#include <CGAL/basic.h>
#if 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
#else
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
#endif

#if 1
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>
typedef CGAL::Segment_Delaunay_graph_filtered_traits_2<Kernel> SDGT;
#else
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
typedef CGAL::Segment_Delaunay_graph_traits_2<Kernel> SDGT;
#endif

#include <CGAL/Segment_Delaunay_graph_2.h>

typedef CGAL::Segment_Delaunay_graph_storage_traits_2<SDGT> ST;
typedef CGAL::Segment_Delaunay_graph_vertex_base_2<ST> Vb;
typedef CGAL::Triangulation_face_base_2<ST> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;

using CGALgraph = CGAL::Segment_Delaunay_graph_2<SDGT, ST, Tds>;

typedef Kernel::Vector_2   Vector_2;
typedef SDGT::Site_2   Site_2;
typedef SDGT::Point_2   Point_2;

// -------------------------------

#include <unordered_map>
#include <algorithm>
#include <unordered_set>

#include "MATGraphBuilder.h"
#include "utils.h"

using Point_handle = CGALgraph::Point_handle;
typedef std::pair<Point_handle, Point_handle> Point_handle_pair;

template<>
class std::less<Point_handle> {
	public:
		bool operator()(const Point_handle & left, const Point_handle & right) const {
			return &(*left) < &(*right);
		}
};

template<>
class std::less<Point_handle_pair> {
	public:
		bool operator()(const Point_handle_pair & left, const Point_handle_pair & right) const {
			if( &(*left.first) < &(*right.first) ) return true;
			if( &(*left.first) > &(*right.first) ) return false;
			return &(*left.second) > &(*right.second);
		}
};

class MATGraphBuilder
{
public:
	std::list<MATvert> verts_;
	void initialize(MATGraph& mat, const Paths & slice);
	void initialize(MATGraph& mat, const Paths & slice, const Paths & cuttings);
	void burnToAshes(); // remove all segments marked as burning / collapsing

protected:
	std::unordered_set<CGALgraph::Face_handle> to_be_removed;

	void markOutsideFaces();

	void burnConcaveVertices();

	void buildMAT(MATGraph& mat);

	bool isInputVertex(CGALgraph::Face_handle f) const;
	Disk voronoiCircumcircle(const CGALgraph::Face_handle fh) const;
	CGALgraph* delaunay = nullptr;
	typedef std::set<Point_handle_pair> InputEdgeSet;
	InputEdgeSet interiorEdges_;

private:
	// helpers:

	bool arcIsDegenerate(const CGALgraph::Face_handle & fh, const int i) const;

	int dumpVandEdges(int base, const Path & poly, std::vector<Point_2> & pts, std::vector<std::pair<int,int>> & indices, bool closed) const;


	Vec2d toVec2d(const Point_2 & p) const { return Vec2d(CGAL::to_double(p.x()), CGAL::to_double(p.y())); }
	Vec2d toVec2d(const Vector_2 & p) const { return Vec2d(CGAL::to_double(p.x()), CGAL::to_double(p.y())); }

	Site toSite(const Site_2 & s) const {
		if (s.is_point()) return Site(toVec2d(s.point()));
		return Site(toVec2d(s.segment().source()), toVec2d(s.segment().target()));
	}

	bool same_points(const Site_2& p, const Site_2& q) const;
	bool is_endpoint_of_segment(const Site_2& p, const Site_2& s) const;
};

void MATGraphBuilder::initialize(MATGraph& mat, const Paths & slice)
{
	if (slice.empty()) return;

	std::vector<Point_2> pts;
	std::vector<std::pair<int,int>> indices;
	int base(0);
	for( const auto & path : slice ) {
		base += dumpVandEdges(base, path, pts, indices, true);
	}
	CGALgraph delaunay_;
	//std::random_shuffle(indices.begin(), indices.end());
	this->delaunay = &delaunay_;
	delaunay->insert_segments(pts, indices.begin(), indices.end());

	markOutsideFaces();
	burnConcaveVertices();

	buildMAT(mat);
	mat.splitApexes();

	this->delaunay = nullptr;
}

void MATGraphBuilder::initialize(MATGraph& mat, const Paths & slice, const Paths & cuttings)
{
	if (slice.empty()) return;

	CGALgraph delaunay_;
	this->delaunay = &delaunay_;

	std::vector<Point_2> pts;
	std::vector<std::pair<int,int>> indices;
	int base(0);

	for( const auto & path : cuttings ) {
		base += dumpVandEdges(base, path, pts, indices, false);
	}
	delaunay->insert_segments(pts, indices.begin(), indices.end());

	interiorEdges_.clear();
	for( CGALgraph::Finite_vertices_iterator fvit = delaunay->finite_vertices_begin();
	  fvit != delaunay->finite_vertices_end(); ++fvit ) {
		auto & ss = fvit->storage_site();
		if( ss.is_point() ) continue;
		Point_handle ph0 = ss.source_of_supporting_site();
		Point_handle ph1 = ss.target_of_supporting_site();
		if( std::less<Point_handle>()(ph1, ph0) ) std::swap(ph0, ph1);
		interiorEdges_.insert(std::make_pair(ph0, ph1));
	}
	pts.clear();
	indices.clear();
	base = 0;
	for( const auto & path : slice ) {
		base += dumpVandEdges(base, path, pts, indices, true);
	}
	delaunay->insert_segments(pts, indices.begin(), indices.end());


	markOutsideFaces();
	burnConcaveVertices();

	buildMAT(mat);
	mat.splitApexes();

	this->delaunay = nullptr;
}

int
MATGraphBuilder::
dumpVandEdges(int base, const Path & poly, std::vector<Point_2> & pts, std::vector<std::pair<int,int>> & indices, bool closed) const {
	int nbv(poly.size());
	if( nbv < 2 ) return 0;
	for( int i = 0; i < nbv-1; ++i ) {
		pts.emplace_back(poly[i].x(), poly[i].y());
		indices.push_back({base + i, base + i+1});
	}
	pts.emplace_back(poly[nbv-1].x(), poly[nbv-1].y());
	if( closed ) {
		indices.push_back({base + nbv-1, base});
	}
	return nbv;
}


void
MATGraphBuilder::
buildMAT(MATGraph& mat)
{
	std::unordered_map<CGALgraph::Face_handle, MATvert *> face_to_vert;
	// First, we create all the vertices of the MATGraph.
	for ( auto it = delaunay->all_faces_begin(); it != delaunay->all_faces_end(); ++it ) {
		if ( to_be_removed.count(it) > 0) continue; // don't copy vertex into MATGraph
		if ( delaunay->is_infinite(it) ) {
			cerr << "FOUND AN INFINITE DELAUNAY FACE THAT IS NOT MARKED FOR REMOVAL !" << endl;
			continue;
		}
		CGALgraph::Face_handle fh = it;
		mat.verts.emplace_back(voronoiCircumcircle(it));
		MATvert& vert = mat.verts.back();
		face_to_vert.emplace(fh, &vert);
	}
	// Then, we connect them with half-arcs / edges
	for ( auto it = delaunay->all_faces_begin(); it != delaunay->all_faces_end(); ++it ) {
		if ( to_be_removed.count(it) > 0) continue; // don't copy vertex into MATGraph
		if ( delaunay->is_infinite(it) ) continue;
		CGALgraph::Face_handle fh = it;
		MATvert * vert_ptr = face_to_vert[fh];
		MATvert & vert = *vert_ptr;
		for (int site_i = 0; site_i < 3; ++site_i)
		{
			CGALgraph::Face_handle nh = fh->neighbor(site_i);
			auto neighbor_it = face_to_vert.find(nh);
			if (neighbor_it == face_to_vert.end()) continue;
			if( nh < fh ) continue;
			// Get the neighbor
			MATvert * neighbor = neighbor_it->second;
			// Create edge from vertex to neighbor
			vert.edges_.emplace_back();
			EdgeIterator edge = --(vert.edges_.end());
			// Create edge from neighbor to vertex
			neighbor->edges_.emplace_back();
			// setup twin_ pointers
			edge->twin_ = neighbor->edges_.end();
			--edge->twin_;
			edge->twin()->twin_ = edge;
			// setup site_
			edge->site_ = toSite(fh->vertex(CGALgraph::ccw(site_i))->site());
			edge->twin()->site_ = toSite(fh->vertex(CGALgraph::cw(site_i))->site());
			// setup to_ pointers
			edge->to_ = neighbor;
			edge->twin()->to_ = vert_ptr;
			// compute arc type
			EdgeType type;
			if( edge->site().is_point() != edge->twin()->site().is_point() )
				type = EdgeType::VertEdge;
			else if( edge->site().is_point() && edge->twin()->site().is_point() )
				type = EdgeType::VertVert;
			else type = EdgeType::EdgeEdge;
			// setup arc type
			edge->twin()->type = edge->type = type;
		}
	}
	// Finally we make sure half-arcs are ordered CCW around each vertex
	for ( auto fh = delaunay->all_faces_begin(); fh != delaunay->all_faces_end(); ++fh ) {
		auto vit = face_to_vert.find(fh);
		if( vit == face_to_vert.end() ) continue;
		MATvert * v_ptr = vit->second;
		MATvert & v = *v_ptr;
		if( v.edges_.size() < 3 ) continue;
		// So we can assume that we have a vertex of degree exactly 3 since CGAL produces at most
		// that.
		CGALgraph::Face_handle f0 = fh->neighbor(0);
		CGALgraph::Face_handle f1 = fh->neighbor(CGALgraph::ccw(0));
		MATvert * n0 = face_to_vert[f0];
		MATvert * n1 = face_to_vert[f1];
		EdgeIterator e0 = v.edges_.begin();
		while( e0->to() != n0 ) e0++;
		EdgeIterator e1 = e0;
		e1++; if( e1 == v.edges_.end() ) e1 = v.edges_.begin();
		if( e1->to() == n1 ) continue; // good ordering
		e0 = v.edges_.begin(); e1 = e0; e1++;
		//cerr << "FIRST AND 2nd = " << &(*v.edges_.begin()) << " and " << &(*(++v.edges_.begin())) << endl;
		v.edges_.splice(e0, v.edges_, e1);
		//cerr << "After swap:\n";
		//cerr << "FIRST AND 2nd = " << &(*v.edges_.begin()) << " and " << &(*(++v.edges_.begin())) << endl;
	}
}

/* [markOutsideFaces] marks as visited the faces (Voronoi vertices) that
 * we do not want : they are
 *  - outside the input polygons with holes.
 * */
void
MATGraphBuilder::
markOutsideFaces() {
	std::vector<CGALgraph::Face_handle> goodFaces;
	struct TraversalData {
		CGALgraph::Face_handle face;
		bool inside;
	};
	std::stack<TraversalData> stack;
	stack.push({delaunay->infinite_face(), false});
	while( ! stack.empty() ) {
		TraversalData td = stack.top();
		stack.pop();
		if ( to_be_removed.count(td.face) > 0 ) continue;
		to_be_removed.emplace(td.face); // mark as visited
		if( td.inside ) {
			goodFaces.push_back(td.face);
		}
		for( int i = 0; i < 3; ++i ) {
			CGALgraph::Face_handle nei = td.face->neighbor(i);
			if( to_be_removed.count(nei) > 0 ) { // the neighbor was already visited
				continue;
			}
			bool degenerate = arcIsDegenerate(td.face, i);
			if( ! degenerate ) {
				stack.push({nei, td.inside});
			} else {
				bool site_is_a_cutting;
				if( interiorEdges_.empty() ) {
					site_is_a_cutting = false;
				} else {
					auto ss = td.face->vertex((i+2)%3)->storage_site();
					if( ss.is_point() )
						ss = td.face->vertex((i+1)%3)->storage_site();
					Point_handle ph0 = ss.source_of_supporting_site();
					Point_handle ph1 = ss.target_of_supporting_site();
					if( std::less<Point_handle>()(ph1, ph0) ) std::swap(ph0, ph1);
					site_is_a_cutting = interiorEdges_.count(std::make_pair(ph0, ph1)) > 0;
				}
				if( site_is_a_cutting ) {
					stack.push({nei, td.inside});
				} else {
					stack.push({nei, ! td.inside});
				}
			}
		}
	}
	// recover the good faces
	for( const auto & fh : goodFaces ) {
		to_be_removed.erase(fh);
	}
}

/* ConVEX vertice sof input polygon are endpoint of the medial axis, so we keep them, but conCAVE
 * one are useless, we mark them are burned.
 */
void
MATGraphBuilder::
burnConcaveVertices() {
	CGALgraph::Finite_faces_iterator fh(delaunay->finite_faces_begin());
	for( ; fh != delaunay->finite_faces_end(); ++fh ) {
		if( to_be_removed.count(fh) > 0 ) continue; // outside vertex already burned
		for( int i = 0; i < 3; ++i ) {
			// If arc is not medial, continue
			if( arcIsDegenerate(fh, i) ) // this outgoing arc is a 'dummy' Voronoi arc, ignore.
				fh->set_neighbor(i, delaunay->infinite_face());
		}
		bool burn = true;
		if( ! isInputVertex(fh) ) { // not an input vertex, so it is an internal Voronoï vertex
			 burn = false;
		} else {
			// fh is an input vertex, remove concave (= not an endpoint of medial axis)
			for( int i = 0; i < 3; ++i ) {
				// If arc is not medial, continue
				if( arcIsDegenerate(fh, i) ) // this outgoing arc is a 'dummy' Voronoi arc, ignore.
					continue;
				const CGALgraph::Face_handle & end = fh->neighbor(i);
				// If arc is going outside, continue, because that means the vertex is concave
				if( to_be_removed.count(end) > 0 ) continue;
				burn = false;
			}
		}
		if( burn ) to_be_removed.emplace(fh);
	}
}



bool
MATGraphBuilder::
arcIsDegenerate(const CGALgraph::Face_handle & fh, const int i) const {
	const CGALgraph::Vertex_handle pv = fh->vertex(CGALgraph::cw(i));
	if( delaunay->is_infinite(pv) ) return false;
	const Site_2 & p = pv->site();
	const CGALgraph::Vertex_handle qv = fh->vertex(CGALgraph::ccw(i));
	if( delaunay->is_infinite(qv) ) return false;
	const Site_2 & q = qv->site();
	return
		( p.is_segment() && q.is_point() && is_endpoint_of_segment(q, p) ) ||
		( p.is_point() && q.is_segment() && is_endpoint_of_segment(p, q) );
}


Disk
MATGraphBuilder::
voronoiCircumcircle(const CGALgraph::Face_handle fh) const {
	Vec2d p = toVec2d(delaunay->circumcenter(fh));
	for( int i = 0; i < 3; ++i ) {
		const auto & s = fh->vertex(i)->site();
		if( s.is_point() ) {
			return Disk(p, (p - toVec2d(s.point())).length());
		}
	}
	const Site_2 & s = fh->vertex(0)->site();
	Vec2d sVec = toVec2d(s.segment().to_vector());
	return Disk(p, Utils::distanceToLine(p, toVec2d(s.source()), sVec));
}

// is this Voronoi vertex a vertex from the input polygon?
bool
MATGraphBuilder::
isInputVertex(CGALgraph::Face_handle f) const {
	if( delaunay->is_infinite(f) ) {
		return false;
	}
	const Site_2 & s0 = f->vertex(0)->site();
	const Site_2 & s1 = f->vertex(1)->site();
	const Site_2 & s2 = f->vertex(2)->site();
	int npts = 0;
	if ( s0.is_point() ) ++npts;
	if ( s1.is_point() ) ++npts;
	if ( s2.is_point() ) ++npts;
	if( npts != 1 ) return false;
	if ( s0.is_point() ) {
		return ( is_endpoint_of_segment(s0, s1) && is_endpoint_of_segment(s0, s2) );
	} else if ( s1.is_point() ) {
		return ( is_endpoint_of_segment(s1, s0) && is_endpoint_of_segment(s1, s2) );
	} else {
		CGAL_assertion( s2.is_point() );
		return ( is_endpoint_of_segment(s2, s0) && is_endpoint_of_segment(s2, s1) );
	}
}

bool
MATGraphBuilder::
same_points(const Site_2& p, const Site_2& q) const {
return delaunay->geom_traits().equal_2_object()(p, q);
}

bool
MATGraphBuilder::
is_endpoint_of_segment(const Site_2& p, const Site_2& s) const
{
CGAL_precondition( p.is_point() && s.is_segment() );
return ( same_points(p, s.source_site()) ||
     same_points(p, s.target_site()) );
}

void buildMATGraphWithCGAL(MATGraph& mat, const Paths & inputPoly) {
	MATGraphBuilder builder;
   	builder.initialize(mat, inputPoly);
}
