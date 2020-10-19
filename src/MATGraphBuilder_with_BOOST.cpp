/* vim: set sts=0 ts=4 sw=4 noet tw=0 : */

#include <boost/polygon/voronoi.hpp>

using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;

#include "MATGraphBuilder.h"
#include "utils.h"

#if ICESL_PLUGIN==1
#include <clipper.hpp>
#else
#include <polyclipping/clipper.hpp>
#endif

#include <stack>

using CLPoint = ClipperLib::IntPoint;
using CLPath = ClipperLib::Path;
using CLPaths = ClipperLib::Paths;

struct Segment {
    Vec2i p0;
    Vec2i p1;
    Segment(const CLPoint & v0, const CLPoint & v1) : p0(v0.X, v0.Y), p1(v1.X, v1.Y) {}
};

namespace boost {
    namespace polygon {

        template <>
            struct geometry_concept<Vec2i> { typedef point_concept type; };

        template <>
            struct point_traits<Vec2i> {
                typedef typename Vec2i::NumberType coordinate_type;

                static inline coordinate_type get(const Vec2i& point, orientation_2d orient) {
                    return (orient == HORIZONTAL) ? point.x() : point.y();
                }
            };

        template <>
            struct geometry_concept<Segment> { typedef segment_concept type; };

        template <>
            struct segment_traits<Segment> {
                typedef typename Vec2i::NumberType coordinate_type;
                typedef Vec2i point_type;

                static inline point_type get(const Segment& segment, direction_1d dir) {
                    return dir.to_int() ? segment.p1 : segment.p0;
                }
            };

    }  // polygon
}  // boost

inline Vec2i convert2d22i(const double scale, const Vec2d & v) {
    Vec2d vv = scale * v;
    return Vec2i(static_cast<int>(round(vv.x())), static_cast<int>(round(vv.y())));
}

using VD = voronoi_diagram<double>;
using VPtr = const VD::vertex_type *;
using EPtr = const VD::edge_type *;
using CPtr = const VD::cell_type *;

void buildMATGraphWithBOOST(MATGraph& mat, const CLPaths & paths, const double mm_per_unit);

void buildMATGraphWithBOOST(MATGraph& mat, const Paths & inputPoly, const double units_per_mm) {
    // For now we ASSUME the input is in MILLImeters and we convert to INTEGER MICROmeter
    if (inputPoly.empty()) return;

    // Scale and fix the input
    CLPaths paths;

    for( const auto & path : inputPoly ) {
        const size_t N = path.size();
		//cerr << "IN PATH OF SIZE " << N << endl;
        if( N < 3 ) continue;
        paths.emplace_back();
        CLPath & p = paths.back();
        for( const auto & pt : path ) {
            Vec2i v = convert2d22i(units_per_mm, pt);
            p.emplace_back(v.x(), v.y());
        }
    }
    ClipperLib::SimplifyPolygons(paths, ClipperLib::pftEvenOdd);
    buildMATGraphWithBOOST(mat, paths, 1.0/units_per_mm);
}

void buildMATGraphWithBOOST(MATGraph& mat, const CLPaths & paths, const double mm_per_unit) {
    if( paths.empty() ) return;

    auto isPointSite = [&](CPtr cell) -> bool {
        // FIXME: is it equivalent to cell->contains_point()?
        // if so, use it.
        if( cell->source_category() == boost::polygon::SOURCE_CATEGORY_SEGMENT_START_POINT) {
            return true;
        } else if( cell->source_category() == boost::polygon::SOURCE_CATEGORY_SEGMENT_END_POINT) {
            return true;
        }
        return false;
    };

    std::vector<Segment> segments;

    auto toSite = [&](CPtr cell) {
        std::size_t index = cell->source_index();
        const Segment & s = segments[index];
        Vec2i p0 = boost::polygon::low(s);
        Vec2d p0d = mm_per_unit * Vec2d(p0.x(), p0.y());
        Vec2i p1 = boost::polygon::high(s);
        Vec2d p1d = mm_per_unit * Vec2d(p1.x(), p1.y());
        if( cell->source_category() == boost::polygon::SOURCE_CATEGORY_SEGMENT_START_POINT) {
            return Site(p0d);
        } else if( cell->source_category() == boost::polygon::SOURCE_CATEGORY_SEGMENT_END_POINT) {
            return Site(p1d);
        }
        return Site(p0d, p1d);
    };

    auto voronoiCircumcircle = [&](const VPtr v) -> Disk {
        Vec2d p = mm_per_unit * Vec2d(v->x(), v->y());
        EPtr first_edge = v->incident_edge();
        EPtr edge = first_edge;
        do {
            CPtr cell = edge->cell();
            std::size_t index = cell->source_index();
            const Segment & s = segments[index];
            Site site = toSite(cell);
            if( site.is_point() )
                return Disk(p, (p-site.point()).length());
            edge = edge->rot_next();
        } while( edge != first_edge );
        CPtr cell = edge->cell();
        std::size_t index = cell->source_index();
        const Segment & s = segments[index];
        Vec2d p0d = mm_per_unit * Vec2d(s.p0.x(), s.p0.y());
        Vec2d p1d = mm_per_unit * Vec2d(s.p1.x(), s.p1.y());
        return Disk(p, Utils::distanceToLine(p, p0d, p1d-p0d));
    };

    auto printV = [=](const VD::vertex_type * v) {
        cerr << '(' << (v->x() * mm_per_unit) << ", " << (v->y() * mm_per_unit) << ')';
    };

    auto analyze = [&](const VD::vertex_type * v, bool & is_input, int & degree, int & n_secondaries) {
		const EPtr first_edge = v->incident_edge();
		EPtr edge = first_edge; // iterate over outgoing edges
		n_secondaries = degree = 0;
		is_input = false;
		do {
			degree++;
			if( edge->is_secondary() ) n_secondaries++;
			CPtr cell = edge->cell();
			std::size_t index = cell->source_index();
			const Segment & s = segments[index];
			if( (v->x() == s.p0.x()) && (v->y() == s.p0.y()) ) is_input = true;
			if( (v->x() == s.p1.x()) && (v->y() == s.p1.y()) ) is_input = true;
            edge = edge->rot_next();
		} while( edge != first_edge );
	};

    // Go
	for( const auto & path : paths ) {
        const size_t N = path.size();
		if( N < 3 ) continue;
		// We always insert loops.
		// So any point is adjacent to an even number of edges (unless degenerate cases).
        CLPoint p0 = path[N-1];
        for( const auto & p : path ) {
            CLPoint p1 = p;
            segments.emplace_back(p0, p1);
            p0 = p1;
        }
    }

	if( segments.empty() ) return;

    VD vd;
    construct_voronoi(segments.begin(), segments.end(), &vd);
    //cerr << "BOOST Voronoi built over " << segments.size() << " input segments\n";

	// Find start vertex
    VPtr curv = nullptr;
    const int outside_flag = 1; // we currently don't know the color of the "outside" vertices. so we store -1
    const int inside_flag = 1 - outside_flag;
	int start_flag;

	for ( auto it = vd.vertices().begin(); it != vd.vertices().end(); ++it ) {
		VPtr v = &(*it);
		EPtr first_edge = v->incident_edge();
		EPtr edge = first_edge; // iterate over outgoing edges
		do {
			if( ! edge->is_infinite() ) {
				edge = edge->rot_next();
				continue;	
			}
			curv = v;
			start_flag = edge->is_secondary() ? inside_flag : outside_flag;
			break;
		} while( edge != first_edge );
		if( curv ) break;
	}

	assert(curv != nullptr);
	//cerr << "START VERTEX IS "; printV(curv); cerr << endl;

    // Flag the vertices
    curv->color(2+start_flag); // color = 2 * visited + in/out flag // visited=1, color=0
    stack<VPtr> visitee;
    while( true ) {
		int degree, n_sec; bool is_input;
		analyze(curv, is_input, degree, n_sec);
		//cerr << "VISITING "; printV(curv); cerr << (is_input?" INPUT":" INTERNAL") << ", deg=" << degree << ", n_sec=" << n_sec << endl;
        EPtr edge = curv->incident_edge();
        int flag = curv->color() % 2;
		if( is_input && (n_sec == 0) ) {
			assert(degree % 2 == 0);
			curv->color(2+inside_flag);
			int fl = -1;
			EPtr start = curv->incident_edge();
			do {
				VPtr nv = start->vertex1();
				assert( nv != nullptr );
				int scol = nv->color();
				if( scol / 2 == 1 ) { // already visited?
					fl = scol % 2;
					break;
				}
				start = start->rot_next();
			} while (start != curv->incident_edge());
			assert( fl != -1 );
			edge = start;
			do {
				VPtr nv = edge->vertex1();
				assert( nv != nullptr );
				if( nv->color() / 2 == 1 ) { // already visited?
					assert( nv->color() % 2 == fl );
					edge = edge->rot_next();
					fl = 1 - fl;
					continue;
				}
				nv->color(2+fl);
				visitee.push(nv);
				fl = 1 - fl;
				edge = edge->rot_next();
			} while (edge != start );
		} else {
			do {
				VPtr nv = edge->vertex1();
				if( nv == curv ) { cerr << "SELF-LOOP ERROR"; exit(-1); }
				if( (nullptr == nv) || (nv->color() / 2 == 1) ) { // |nv| infinite or already visited?
					edge = edge->rot_next();
					continue;
				}
				visitee.push(nv);
				if( edge->is_primary() ) {
					nv->color(2 + flag);
				} else {
					//cerr << "SECONDARY EDGE from "; printV(curv);
					//cerr << " to "; printV(nv); cerr << endl;
					nv->color(2 + 1 - flag);
				}
				edge = edge->rot_next();
			} while (edge != curv->incident_edge());
		}
        if( visitee.empty() ) break;
        curv = visitee.top();
        visitee.pop();
    }
    //cerr << "Outside flag is " << outside_flag << endl;
    if( outside_flag == -1 ) { cerr << "outside flag ERROR"; return; }

    // Build the MAT
    std::map<VPtr, MATvert *> mat_vertex;
    int i(0);
	for ( auto it = vd.vertices().begin(); it != vd.vertices().end(); ++it ) {
		VPtr v = &(*it);
        if ( v->color() % 2 != inside_flag ) continue;

        mat.verts.emplace_back(voronoiCircumcircle(v));
        MATvert& vert = mat.verts.back();
		mat_vertex[v] = & vert;
		//cerr << "V at " << vert.circumcircle.center_ << " with radius " << vert.circumcircle.radius_ << endl;
	}
	for ( auto it = vd.vertices().begin(); it != vd.vertices().end(); ++it ) {
		VPtr v = &(*it);
		if ( v->color() % 2 != inside_flag ) continue;

		MATvert * matv = mat_vertex[v];

		EPtr first_edge = v->incident_edge();
		EPtr edge = first_edge; // iterate over outgoing edges
		do {
			auto neighbor_it = mat_vertex.find(edge->vertex1());
			if( (!edge->is_primary()) || (neighbor_it == mat_vertex.end()) ) {
				edge = edge->rot_next();
				continue;
			}
			MATvert * neighbor = neighbor_it->second;
			// create edge
			matv->edges_.emplace_back();
			MATedge & mate = matv->edges_.back();
			// setup mate.to_ pointer to neighbor
			mate.to_ = neighbor;
			// setup mate.type
			CPtr cell = edge->twin()->cell(); // the cell to the right of the edge
			CPtr opp_cell = edge->cell();
			if( isPointSite(cell) != isPointSite(opp_cell) ) {
				mate.type = EdgeType::VertEdge;
			} else if( isPointSite(cell) ) {
				mate.type = EdgeType::VertVert;
			} else {
				mate.type = EdgeType::EdgeEdge;
			}
			// setup mate.site_
			mate.site_ = toSite(cell);
			// setup mate.twin_
			if( ! neighbor->edges_.empty() ) { // the neighbor has already created its edges, let's find our twin
				bool found = false;
				for( EdgeIterator twin_candidate = neighbor->edges_.begin();
				  twin_candidate != neighbor->edges_.end(); ++twin_candidate ) {
					if( twin_candidate->to_ == matv ) {
						twin_candidate->twin_ = --matv->edges_.end();
						mate.twin_ = twin_candidate;
						found = true;
					}
				}
				if( ! found ) {
					cerr << "TWIN NOT FOUND !\n";
					exit(-1);
				}
			}
			edge = edge->rot_next();
		} while( edge != first_edge );
    }
	mat.splitApexes();
}
