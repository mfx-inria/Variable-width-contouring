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
	void splitApexes(MATGraph& mat);

	bool isInputVertex(CGALgraph::Face_handle f) const;
	Disk voronoiCircumcircle(const CGALgraph::Face_handle fh) const;
	CGALgraph* delaunay = nullptr;
	typedef std::set<Point_handle_pair> InputEdgeSet;
	InputEdgeSet interiorEdges_;

private:
	// helpers:

	bool arcIsEE(const CGALgraph::Face_handle & fh, const int i) const;

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
	splitApexes(mat);

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
	splitApexes(mat);

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
	std::unordered_map<CGALgraph::Face_handle, MATvert*> face_to_vert;
	for ( auto it = delaunay->all_faces_begin(); it != delaunay->all_faces_end(); ++it ) {
		if ( to_be_removed.count(it) > 0) continue; // don't copy edge into MATGraph

		if ( delaunay->is_infinite(it) ) continue;

		CGALgraph::Face_handle fh = it;

		mat.verts.emplace_back(voronoiCircumcircle(it));
		MATvert& vert = mat.verts.back();

		for (int site_i = 0; site_i < 3; ++site_i)
		{
			CGALgraph::Face_handle nh = fh->neighbor(site_i);
			if (to_be_removed.count(nh) > 0) {
				            vert.edge[site_i].to_ = nullptr;
				            vert.edge[site_i].twin_ = nullptr;
			} else {
				            vert.edge[site_i].site_ = toSite(fh->vertex(CGALgraph::ccw(site_i))->site());
				auto neighbor_it = face_to_vert.find(nh);
				if (neighbor_it != face_to_vert.end())
				{
					MATvert* neighbor_vert = neighbor_it->second;
					vert.edge[site_i].twin_ = &neighbor_vert->edge[nh->index(fh)];
					vert.edge[site_i].twin()->twin_ = &vert.edge[site_i];
					vert.edge[site_i].to_ = neighbor_vert;
					vert.edge[site_i].twin()->to_ = &vert;
					EdgeType type;
					if (vert.edge[site_i].site().is_point() != vert.edge[site_i].twin()->site().is_point()) type = EdgeType::VertEdge;
					else if (vert.edge[site_i].site().is_point() && vert.edge[site_i].twin()->site().is_point()) type = EdgeType::VertVert;
					else type = EdgeType::EdgeEdge;
					vert.edge[site_i].twin()->type = vert.edge[site_i].type = type;
				} else {
					vert.edge[site_i].to_ = nullptr;
				}
			}
		}
		face_to_vert.emplace(fh, &vert);
	}
}


void
MATGraphBuilder::
splitApexes(MATGraph& mat)
{
	for ( MATvert & vert : mat.verts ) {
		for ( int i = 0 ; i < 3 ; ++i ) {
			if ( ! vert.has_neighbor(i) ) continue;
			MATedge & edge = vert.edge[i];
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
				mat.insert(nullptr, edge, Disk(mid, std::fabs(radius)));
			}
		}
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


bool
MATGraphBuilder::
arcIsEE(const CGALgraph::Face_handle & fh, const int i) const {
	if( delaunay->is_infinite(fh) ) {
		return false;
	}
	const Site_2 & p = fh->vertex(CGALgraph::cw(i))->site();
	const Site_2 & q = fh->vertex(CGALgraph::ccw(i))->site();
	return ( p.is_segment() && q.is_segment() );
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
