#include "VWCInfillerPlugin.h"
#include "Path.h"

#include "VariableWidthContouring.h"

using CLint = ClipperLib::cInt;
using CLPoint = ClipperLib::IntPoint;
using CLPath = ClipperLib::Path;
using CLPaths = ClipperLib::Paths;

void buildMATGraphWithBOOST(MATGraph& mat, const CLPaths & inputPoly, const double mm_per_unit);

// -----------------------------------------------

// --- Shared library interface ---
IceSLInterface::IPlugin* createPlugin()
{
    return new VWCInfillerPlugin();
}

void destroyPlugin(IceSLInterface::IPlugin* plugin)
{
    delete plugin;
}
// ----

// -----------------------------------------------

void VWCInfiller::setNumSlices(int num)
{

}

// -----------------------------------------------

void VWCInfiller::setWorldSpace(v3f world_corner_mm, v3f extent_mm)
{

}

// -----------------------------------------------

void VWCInfiller::prepareInfill(int brush)
{
    Sampling::len_threshold = 0.01; // mm
}

// -----------------------------------------------

void VWCInfiller::terminateInfill(int brush)
{

}

// -----------------------------------------------

void VWCInfiller::prepareInfillForSlice(int id, const AAB<2, int>& xy_slice_box, float height, int brush)
{
    cerr << "VWC PREPARE INFILL FOR SLICE " << id << "... ";
    /// retrive extruder used for the infill by this brush
    int extruder_id = 0;
    {
        auto s = m_EnumerableSettings->getSettingByName("infill_extruder_" + std::to_string(brush));
        if (s == nullptr) {
            throw Fatal("This plugin is only for FDM");
        }
        s->getValue(extruder_id);
    }
    /// retrieve extruder nozzle diameter
    float nd = 0.0f;
    {
        auto s = m_EnumerableSettings->getSettingByName("nozzle_diameter_mm_" + std::to_string(extruder_id));
        if (s == nullptr) {
            throw Fatal("This plugin is only for FDM");
        }
        s->getValue(nd);
    }
    nozzle_diameter_ = nd;
    minBeadWidth_ = 0.75 * nozzle_diameter_; // mm
    maxBeadWidth_ = 2.5 * nozzle_diameter_; // mm
    /// retrieve processing scale factor
    float scale = 0.0f;
    {
        auto s = m_EnumerableSettings->getSettingByName("xy_mm_per_int");
        if (s == nullptr) {
            throw Fatal("Missing setting");
        }
        s->getValue(scale);
    }
    xy_mm_per_unit_ = scale;

    cerr << "UNIT PER MM IS " << (1.0 / xy_mm_per_unit_) << "\n";
    cerr << "MIN NEAD WIDTH (MM) IS " << (minBeadWidth_) << "\n";
    cerr << "MAX NEAD WIDTH (MM) IS " << (maxBeadWidth_) << "\n";
    /*// retrieve line width
    float line_width_mm = 0.0f;
    {
        auto s = m_EnumerableSettings->getSettingByName("line_width_mm_" + std::to_string(brush));
        if (s == nullptr) {
            throw Fatal("Missing setting");
        }
        s->getValue(line_width_mm);
    }*/
    cerr << "DONE\n";
}

// ----------------------------------------------- 

bool VWCInfiller::generateInfill(int slice_id, float slice_height_mm, int brush,
        const ClipperLib::Paths & surface,
        std::vector<std::unique_ptr<IceSLInterface::IPath> > & fills,
        bool & preserve_order)
{
    cerr << "VWC GO-GO_GADGET-AU FILL for slice " << slice_id << "... ";
    for( const auto & pp : surface ) {
        cerr << "\n path: ";
        for( const auto & p : pp ) {
            cerr << '(' << p.X << ", " << p.Y << ") ";
        }
    }
    cerr << endl;
    // work data
    MATGraph mat_;
    BoundaryCircles maoi_[2];

    std::vector<Disk> cheekPrecursors;

    // INITIALIZATION

    vector<vector<Sample>> pts;

    using Components = MATGraph::Components;
    using Component  = MATGraph::Component;
    Components components;

    MATGraph mat;
    // FOR DEBUG (dow we want to keep it in prod?
    //CLPaths slice = surface;
    //ClipperLib::SimplifyPolygons(slice, ClipperLib::pftEvenOdd);
    //buildMATGraphWithBOOST(mat, slice, xy_mm_per_unit_);
    // END OF DEBUG
    buildMATGraphWithBOOST(mat, surface, xy_mm_per_unit_);
    mat.print();
    mat.debugCheck();
    if( mat.verts.empty() ) {
        cerr << " NO MEDIAL AXIS ";
        return false;
    }

    const double minToolRadius = minBeadWidth_ / 2.0;
    const double maxToolRadius = maxBeadWidth_ / 2.0;

    VariableWidthContouring sdg(mat, minToolRadius, maxToolRadius, simplificationThreshold_);
    sdg.setWorkingBoundaryCircles(& maoi_[0]);

    // FILTER AND SIMPLIFY THE INPUT

    components.clear();
    mat.computeConnectedComponents(components);
    sdg.filterSmallFeatures(minBeadWidth_, components);
    mat.debugCheck();

    sdg.simplifyMACombinatorics();

    size_t nbSamples(0); // total number of samples generated along all print paths

    double totalLength = 0.0;

    auto insertPoint = [=](const Sample & p, CLPath & path, std::vector<double> & bw) -> bool {
        CLint x = static_cast<CLint>(round(p.pos.x() / xy_mm_per_unit_));
        CLint y = static_cast<CLint>(round(p.pos.y() / xy_mm_per_unit_));
        if( (! path.empty()) && (path.back().X == x) && (path.back().Y == y) ) return false;
        path.push_back({x,y});
        bw.push_back(p.radius*2.0);
        return true;
    };

    auto closePath = [](CLPath & path, std::vector<double> & bw) {
        if( path.size() < 3 ) { path.clear(); return; }
        const CLPoint & first = path[0];
        const CLPoint & last = path.back();
        if( first != last ) {
            path.push_back(first);
            bw.push_back(bw.front());
        }
    };

    // GO!

    int step = -1;
    while( step < 5 ) {
        ++step;
        // Compute connected components. We work on each separately, so that the uniform offset can be customized to each.
        components.clear();
        mat.computeConnectedComponents(components);
        if( 0 == components.size() ) {
            break;
        }
        else {
            if( verbose_ ) cerr << ((step%10) ? '.':'|');
        }
        // Update maoi and the smooth contour
        maoi_[0].clear();
        for( Component & comp : components ) {
            // Since we don't draw anything to PDf, the SmoothPaths are
            // actually useless. But |computeSmoothPaths()| also fills the
            // BoundaryCircles data structure (maoi_[]) which *is* used to
            // sample the print paths. So we call it. The additional useless
            // comput. time is negligible.
            sdg.computeSmoothPaths(comp, comp.outerSmoothPaths); // also fills maoi_[0]
        }

        // 1. TRIMMING
        sdg.shrinkMA_trim(components); // works on maoi_[0] and updates it when an MA edge is subdivided

        // 2. COLLAPSE
        // works on maoi_[0] and updates it when an MA edge is subdivided
        bool hadSomeCollapse = sdg.shrinkMA_collapse(components, cheekPrecursors, /* sharpCut */ false);
        /*if( hadSomeCollapse ) {

        // 3. CHEEKS
        if( ! cheekPrecursors.empty() ) {
        if( verbose_ ) std::cerr << "MODELING CHEEKS..." << endl;
        sdg.shrinkMA_cheeks(components); // works on maoi_[0] and updates it when an MA edge is subdivided
        }
        }*/

        // 4. GLOBAL OFFSET
        sdg.shrinkMA_globalOffset(components); // works on maoi_[0] and updates it when an MA edge is subdivided
        sdg.simplifyMACombinatorics(&components);

        maoi_[1] = maoi_[0];  // Save the |maoi| of the outer contour

        // 5. CLIPPING & SHAVE
        bool someEdgesWereEradicated = sdg.simplifyCollapsedMAGeometry(&components);
        sdg.clipAndShaveCollapsedParts(components, /* no shaving */ false);
        if( someEdgesWereEradicated ) {
            //cerr << "(some edges eradicated)";
            maoi_[1].clear();
            for( Component & comp : components ) {
                sdg.computeSmoothPaths(comp, comp.simpleOuterSmoothPaths, /*walk_on_fire*/true, &maoi_[1]);
            }
        }
        sdg.simplifyMACombinatorics(&components);

        // Re-compute smoothPaths and maoi for inner contour
        maoi_[0].clear();
        for( Component & comp : components ) {
            sdg.computeSmoothPaths(comp, comp.innerSmoothPaths); // also fills maoi_[0]
        }
        cerr << "A ";

        pts.clear();
        for( Component & comp : components ) {
            sdg.samplePrintPath(comp, & maoi_[1], pts);
        }
        cerr << "B ";
        if( ! hadSomeCollapse ) {
            for( Component & comp : components ) {
                comp.voids.clear();
            }
        }
        cheekPrecursors.clear();

        cerr << "F ";
        // output
        CLPath clipperPath;
        std::vector<double> beadWidth;
        for( const auto & samples : pts ) {
            nbSamples += samples.size();
            if( samples.empty() ) continue;
            cerr << "G ";
            Vec2d q = samples.back().pos;
            for( const auto & p : samples ) {
                if( p.radius > maxToolRadius * 1.01 )
                    cerr << "ERROR : output extrusion radius is too large: "
                        << p.radius << " >> " << maxToolRadius << " at " << p.pos << endl;
                if( p.radius < minToolRadius / 1.01 )
                    cerr << "ERROR : output extrusion radius is too small: "
                        << p.radius << " << " << minToolRadius << " at " << p.pos << endl;
                //outpath << p.pos.x() << ' ' << p.pos.y() << ' ' << p.radius << ' ' << p.tangent.x() << ' ' << p.tangent.y() << endl;
                totalLength += (q-p.pos).length();
                q = p.pos;
                if( insertPoint(p, clipperPath, beadWidth) ) {
                }
            }
            closePath(clipperPath, beadWidth);
            fills.push_back(std::unique_ptr<IceSLInterface::IPath>(new IceSLInterface::Path()));
            fills.back()->createPath(clipperPath, true);
            // customize flow
            int flow_idx = fills.back()->addPerVertexAttribute("flow_multiplier");
            for( int i = 0; i < beadWidth.size(); ++i ) {
                fills.back()->setPerVertexAttributeValue(flow_idx, i, beadWidth[i]);
            }
            clipperPath.clear();
            beadWidth.clear();
            cerr << "In step " << step << ", now up to " << nbSamples << " samples\n";
#if 0 // Check angle formed by 3 consecutive samples
                const int N = samples.size();
                for( int i1 = 0; i1 < N; ++i1 ) {
                    int i0 = (i1+N-1)%N;
                    int i2 = (i1+1)%N;
                    const Vec2d & p0 = samples[i0].pos;
                    const Vec2d & p1 = samples[i1].pos;
                    const Vec2d & p2 = samples[i2].pos;
                    const Vec2d a = p1 - p0;
                    const Vec2d b = p2 - p1;
                    const double dx = a.dot(b);
                    const double dy = det(a,b);
                    double angle = std::abs(atan2(dy, dx));
                    if( angle > 1.1*3.14159/2.0 ) {
                        cout << setprecision(10) << std::fixed <<
                            "VERY BIG ANGLE " << angle << ' '<< dx << ' '<<dy<< " FOUND IN " << gOutputFile << " at pos " << p1 << std::endl
                            << i0 << ": " << p0 << endl
                            << i1 << ": " << p1 << endl
                            << i2 << ": " << p2 << endl;
                    }
                }
#endif
        }

        pts.clear();

        mat.removeDestroyed();
    }
    if( verbose_ ) cerr << endl << nbSamples << " samples.\n";
    cerr << "DONE\n";
    return true;
}

// -----------------------------------------------

bool VWCInfillerPlugin::addInfillerSettings(IceSLInterface::EnumerableSettingsInterface& enumerable_settings)
{/*
    // allocates an array with one value per brush
    m_LineWidth_mm.resize(enumerable_settings.getNumberOfBrushes(),0.4f);
    // generates setting per-brush for this infiller
    for (int i = 0; i < (int)m_LineWidth_mm.size(); i++) {
    std::unique_ptr<IceSLInterface::EnumerableSettingsInterface::SettingInterface> s = 
    enumerable_settings.addSettingFromPlugin(
    &m_LineWidth_mm[i],
    &m_LineWidth_min_mm,
    &m_LineWidth_max_mm,
    "line_width_mm_" + std::to_string(i),
    "Line width",
    "Brush_" + std::to_string(i),
    "Specifies the line width to be used during contouring. The flow is adjusted accordingly.",
    1000, // rank to order multiple settings, lower appears before
    [this, &enumerable_settings, i]() -> bool { // show setting only when the infiller is selected
    std::string infiller;
    enumerable_settings.getSettingByName("infill_type_" + std::to_string(i))->getValue(infiller);
    return infiller == this->name();
    }
    );
    if (s == nullptr) {
    return false;
    }
    }
  */
    return true;
}

// -----------------------------------------------
