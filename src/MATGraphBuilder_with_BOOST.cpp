
#include <boost/polygon/voronoi.hpp>

using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;

#include "MATGraphBuilder.h"
#include "utils.h"

#include <polyclipping/clipper.hpp>

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

Vec2i convert2d22i(const Vec2d & v) {
    Vec2d vv = 1000.0 * v;
    return Vec2i(static_cast<int>(round(vv.x())), static_cast<int>(round(vv.y())));
}

void buildMATGraphWithBOOST(MATGraph& mat, const Paths & inputPoly) {
    // For now we ASSUME the input is in MILLImeters and we convert to INTEGER MICROmeter
    if (inputPoly.empty()) return;

    // Clean the input
    CLPaths paths;

    for( const auto & path : inputPoly ) {
        const size_t N = path.size();
        if( N < 3 ) continue;
        paths.emplace_back();
        CLPath & p = paths.back();
        for( const auto & pt : path ) {
            Vec2i v = convert2d22i(pt);
            p.emplace_back(v.x(), v.y());
        }
    }
    ClipperLib::SimplifyPolygons(paths, ClipperLib::pftNonZero);

    // Go
    std::vector<Segment> segments;
    for( const auto & path : paths ) {
        const size_t N = path.size();
        CLPoint p0 = path[N-1];
        for( const auto & p : path ) {
            CLPoint p1 = p;
            segments.emplace_back(p0, p1);
            p0 = p1;
        }
    }
    using VD = voronoi_diagram<double>;
    VD vd;
    construct_voronoi(segments.begin(), segments.end(), &vd);

    // Color the vertices
    int color = 0;
    int outside_flag = -1;
    using VPtr = const VD::vertex_type *;
    VPtr curv = & vd.vertices().front();
    curv->color(2); // color = 2 * visited + in/out flag
    stack<VPtr> visitee;
    while( true ) {
        const voronoi_diagram<double>::edge_type * edge = curv->incident_edge();
        int flag = curv->color() % 2;
        do {
            VPtr nv = edge->vertex1();
            if( nv == curv ) { cerr << "ERROR"; continue; }
            if( nullptr == nv ) {
                int new_outside = edge->is_primary() ? flag : (1-flag);
                //cerr << "Outside would be " << new_outside << " at " << curv->x() << ", " << curv->y()
                //    << endl;
                if( outside_flag != -1 && outside_flag != new_outside ) {
                    cerr << "YAK " << outside_flag << ' ' << new_outside << endl;
                }
                outside_flag = new_outside;
                edge = edge->rot_next();
                continue;
            }
            if( nv->color() / 2 == 1 ) {
                edge = edge->rot_next();
                continue;
            }
            visitee.push(nv);
            if( edge->is_primary() ) {
                nv->color(2 + flag);
            } else {
                nv->color(2 + 1 - flag);
            }
            edge = edge->rot_next();
        } while (edge != curv->incident_edge());
        if( visitee.empty() ) break;
        curv = visitee.top();
        visitee.pop();
    }
    //cerr << "Outside flag is " << outside_flag << endl;
    if( outside_flag == -1 ) { cerr << "outside flag ERROR"; return; }
    int inside = 1 - outside_flag;

    // Build the MAT
    /*
    for ( auto it = vd.vertices().begin(); it != vd.vertices().end(); ++it ) {
        VPtr v = &(*it);
        if ( v->color() % 2 != inside ) continue;

        Vec2d p(v->x(), v->y());
        mat.verts.emplace_back(p);
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
    */
}
