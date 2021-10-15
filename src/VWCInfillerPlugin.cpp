#include "VWCInfillerPlugin.h"
#include "Path.h"
#include "PathTypeBitField.h"

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
    std::unique_ptr<IceSLInterface::EnumerableSettingsInterface::SettingInterface> s;
    s = m_EnumerableSettings->getSettingByName("vwc_overtakes_icesl_"+std::to_string(brush));
    if (s == nullptr) { throw Fatal("Missing setting vwc_overtakes_icesl_"); }
    s->getValue(overtakeIcesl_);
    if( ! overtakeIcesl_ ) return;
    s = m_EnumerableSettings->getSettingByName("path_priority_use_default");
    if (s == nullptr) { overtakeIcesl_ = false; throw Fatal("Missing setting path_priority_use_default"); }
    s->getValue(path_priority_use_default);
    s = m_EnumerableSettings->getSettingByName("path_priority_perimeter");
    if (s == nullptr) { overtakeIcesl_ = false; throw Fatal("Missing setting path_priority_perimeter"); }
    s->getValue(path_priority_perimeter);
    s = m_EnumerableSettings->getSettingByName("path_priority_shell");
    if (s == nullptr) { overtakeIcesl_ = false; throw Fatal("Missing setting path_priority_shell"); }
    s->getValue(path_priority_shell);
    s = m_EnumerableSettings->getSettingByName("path_priority_infill");
    if (s == nullptr) { overtakeIcesl_ = false; throw Fatal("Missing setting path_priority_infill"); }
    s->getValue(path_priority_infill);
    s = m_EnumerableSettings->getSettingByName("path_priority_bridge");
    if (s == nullptr) { overtakeIcesl_ = false; throw Fatal("Missing setting path_priority_bridge"); }
    s->getValue(path_priority_bridge);
    s = m_EnumerableSettings->getSettingByName("cover_thickness_mm_"+std::to_string(brush));
    if (s == nullptr) { overtakeIcesl_ = false; throw Fatal("Missing setting cover_thickness_mm_X"); }
    s->getValue(cover_thickness_mm);
    s = m_EnumerableSettings->getSettingByName("print_perimeter_"+std::to_string(brush));
    if (s == nullptr) { overtakeIcesl_ = false; throw Fatal("Missing setting print_perimeter_X"); }
    s->getValue(print_perimeter);
    s = m_EnumerableSettings->getSettingByName("num_shells_"+std::to_string(brush));
    if (s == nullptr) { overtakeIcesl_ = false; throw Fatal("Missing setting num_shells_X"); }
    s->getValue(num_shells);
    m_EnumerableSettings->getSettingByName("path_priority_use_default")->setValue(false);
    m_EnumerableSettings->getSettingByName("path_priority_perimeter")->setValue(2);
    m_EnumerableSettings->getSettingByName("path_priority_shell")->setValue(3);
    m_EnumerableSettings->getSettingByName("path_priority_infill")->setValue(4);
    m_EnumerableSettings->getSettingByName("path_priority_bridge")->setValue(1);
    m_EnumerableSettings->getSettingByName("cover_thickness_mm_"+std::to_string(brush))->setValue(0.0f);
    m_EnumerableSettings->getSettingByName("print_perimeter_"+std::to_string(brush))->setValue(false);
    m_EnumerableSettings->getSettingByName("num_shells_"+std::to_string(brush))->setValue(0);
}

// -----------------------------------------------

void VWCInfiller::terminateInfill(int brush)
{
    if( ! overtakeIcesl_ ) return;
    m_EnumerableSettings->getSettingByName("path_priority_use_default")->setValue(path_priority_use_default);
    m_EnumerableSettings->getSettingByName("path_priority_perimeter")->setValue(path_priority_perimeter);
    m_EnumerableSettings->getSettingByName("path_priority_shell")->setValue(path_priority_shell);
    m_EnumerableSettings->getSettingByName("path_priority_infill")->setValue(path_priority_infill);
    m_EnumerableSettings->getSettingByName("path_priority_bridge")->setValue(path_priority_shell);
    m_EnumerableSettings->getSettingByName("cover_thickness_mm_"+std::to_string(brush))->setValue(cover_thickness_mm);
    m_EnumerableSettings->getSettingByName("print_perimeter_"+std::to_string(brush))->setValue(print_perimeter);
    m_EnumerableSettings->getSettingByName("num_shells_"+std::to_string(brush))->setValue(num_shells);
}

// -----------------------------------------------

void VWCInfiller::prepareInfillForSlice(int id, const AAB<2, int>& xy_slice_box, float, double, int brush)
{
    // retrieve extruder used for the infill by this brush
    {
        auto s = m_EnumerableSettings->getSettingByName("infill_extruder_" + std::to_string(brush));
        if (s == nullptr) {
            throw Fatal("This plugin is only for FDM");
        }
        s->getValue(extruder_id_);
    }
    float v;
    // retrieve extruder nozzle diameter
    {
        auto s = m_EnumerableSettings->getSettingByName("nozzle_diameter_mm_" + std::to_string(extruder_id_));
        if (s == nullptr) {
            throw Fatal("This plugin is only for FDM");
        }
        s->getValue(v); nozzle_diameter_ = v;
    }
    {
        auto s = m_EnumerableSettings->getSettingByName("vwc_num_beads_" + std::to_string(extruder_id_));
        if (s == nullptr) {
            throw Fatal("This plugin is only for FDM");
        }
        s->getValue(num_beads_);
    }
    {
        auto s = m_EnumerableSettings->getSettingByName("vwc_bead_min_width_mm_" + std::to_string(extruder_id_));
        if (s == nullptr) {
            throw Fatal("Fatality kill (VWC min)");
        }
        s->getValue(v); minBeadWidth_ = v;
    }
    {
        auto s = m_EnumerableSettings->getSettingByName("vwc_bead_max_width_mm_" + std::to_string(extruder_id_));
        if (s == nullptr) {
            throw Fatal("Fatality kill (VWC max)");
        }
        s->getValue(v); maxBeadWidth_ = v;
    }
    if( maxBeadWidth_ < 2.0 * minBeadWidth_ ) {
        throw Fatal("Fatality kill (VWC 2*min>max)");
    }
    // retrieve processing scale factor
    {
        auto s = m_EnumerableSettings->getSettingByName("xy_mm_per_int");
        if (s == nullptr) {
            throw Fatal("Missing setting");
        }
        s->getValue(v); xy_mm_per_unit_ = v;
    }
}

// ----------------------------------------------- 

void printPaths(const CLPaths & paths) {
    for( auto & path : paths ) {
        cerr << path.size() << '\n';
        for( const auto & p : path ) {
            cerr << p.X << ' ' << p.Y << '\n';
        }
    }
}

//static void cleanInput(const CLPaths & surface, CLPaths & clean, bool verbose = false) {
//    CLPaths tempPaths = surface;
//    if( verbose ) {
//        cerr << "----------------------------------\n";
//        for( auto & path : surface ) {
//            cerr << "INPUT PATH :\n" << path.size() << '\n';
//            for( const auto & p : path ) {
//                cerr << p.X << ' ' << p.Y << '\n';
//            }
//        }
//        cerr << " - - - - - - - - - - -> \n";
//    }
//    for( const CLPath & temp : tempPaths ) {
//        if( temp.size() <= 2 ) continue;
//        clean.emplace_back();
//        CLPath & path = clean.back();
//        path.push_back(temp[0]);
//        int i(0), j(1), k(2), numV(temp.size());
//        while( j > i ) {
//            if( temp[i] == temp[j] ) { // duplicate vertices
//                cerr << "KILLED a Duplicate vertex: bad input file" << endl;
//                j = k;
//                k = (k+1) % numV;
//                continue;
//            }
//            ClipperLib::IntPoint v0 = subPoints(temp[j], temp[i]);
//            ClipperLib::IntPoint v1 = subPoints(temp[k], temp[j]);
//            double lenProd = cautiousLength(v0) * cautiousLength(v1);
//            double ddot, ddet;
//#ifndef use_int32
//            if (true) {//UseFullInt64Range) {
//                ddot = dot128(v0, v1);
//                ddet = det128(v0, v1);
//            }
//            else
//#endif
//            {
//                ddot = dot(v0, v1);
//                ddet = det(v0, v1);
//            }
//            if( (ddot < (-0.999999)*lenProd) && (std::abs(ddet) < 1e-6*lenProd) ) {
//                cerr << "KILLED a U-turn: bad input file" << endl;
//                j = k;
//                k = (k+1)%numV;
//                continue;
//            }
//            if( (ddot > 0.999999*lenProd) && (std::abs(ddet) < 1e-6*lenProd) ) {
//                cerr << "KILLED a colinearity at " << temp[j] << ": bad input file" << endl;
//                j = k;
//                k = (k+1)%numV;
//                continue;
//            }
//            path.push_back(temp[j]);
//            i = j;
//            j = k;
//            k = (k+1)%numV;
//        }
//        if( verbose ) {
//            cerr << path.size() << '\n';
//            for( const auto & p : path ) {
//                cerr << p.X << ' ' << p.Y << '\n';
//            }
//        }
//    }
//        if( verbose ) {
//            cerr << "=========================\n";
//        }
//    ClipperLib::CleanPolygons(clean);
//}

// ----------------------------------------------- 

bool VWCInfiller::generateInfill(int slice_id, float layer_height_mm, double layer_thickness_mm, int brush,
        const ClipperLib::Paths & surface,
        std::vector<std::unique_ptr<IceSLInterface::IPath> > & fills,
        bool & preserve_order, ClipperLib::Paths & fallback_surface)
{
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
#if 1
    CLPaths input = surface;
    bool debug = false;//slice_id == 54;
    if( debug ) {
        cerr << "SLICE " << slice_id << ", input paths:\n";
        //printPaths(surface);
    }
    ClipperLib::CleanPolygons(input);
    if( debug ) {
        cerr << "Cleaned paths (mm_per_unit = "<<xy_mm_per_unit_<<"):\n";
        printPaths(input);
    }
    buildMATGraphWithBOOST(mat, input, xy_mm_per_unit_);
#else
    buildMATGraphWithBOOST(mat, surface, xy_mm_per_unit_);
#endif
    mat.debugCheck();
    if( mat.verts.empty() ) {
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

    auto insertPoint = [=](const Sample & p, CLPath & path, std::vector<double> * bw = nullptr) -> bool {
        CLint x = static_cast<CLint>(round(p.pos.x() / xy_mm_per_unit_));
        CLint y = static_cast<CLint>(round(p.pos.y() / xy_mm_per_unit_));
        if( (! path.empty()) && (path.back().X == x) && (path.back().Y == y) ) return false;
        path.push_back({x,y});
        if( ! bw ) return true;
        bw->push_back(p.radius*2.0/nozzle_diameter_);
        return true;
    };

    auto closePath = [](CLPath & path, std::vector<double> * bw = nullptr) {
        if( path.size() < 3 ) { path.clear(); if( bw ) { bw->clear(); } return false; }
        const CLPoint & first = path[0];
        const CLPoint & last = path.back();
        if( first != last ) {
            path.push_back(first);
            if( bw ) {
                bw->push_back(bw->front());
            }
        }
        return true;
    };

    // GO!

    int step = 0;
    const int numSteps = num_beads_ == 0 ? 5000 : num_beads_;
    while( step < numSteps ) {
        ++step;
        // Compute connected components. We work on each separately, so that the uniform offset can be customized to each.
        components.clear();
        mat.computeConnectedComponents(components);
        if( 0 == components.size() ) {
            break;
        }
        else {
            if( debug ) {
                //cerr << "VWC STEP " << step << endl;
            }
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

        pts.clear();
        for( Component & comp : components ) {
            sdg.samplePrintPath(comp, & maoi_[1], pts);
        }
        if( ! hadSomeCollapse ) {
            for( Component & comp : components ) {
                comp.voids.clear();
            }
        }
        cheekPrecursors.clear();

        // output
        CLPath clipperPath;
        std::vector<double> beadWidth;
        const auto pathType = overtakeIcesl_ ? (step == 1 ? PathType::e_Perimeter : PathType::e_Shell) :
            PathType::e_Infill;
        for( const auto & samples : pts ) {
            nbSamples += samples.size();
            if( samples.empty() ) continue;
            Vec2d q = samples.back().pos;
            for( const auto & p : samples ) {
                if( p.radius > maxToolRadius * 1.01 )
                    cerr << "[Variable-width contouring] ERROR at slice " << slice_id << " step " << step << ": output extrusion radius is too large: "
                        << p.radius << " >> " << maxToolRadius << " at " << p.pos << endl;
                if( p.radius < minToolRadius / 1.01 )
                    cerr << "[Variable-width contouring] ERROR at slice " << slice_id << " step " << step << " : output extrusion radius is too small: "
                        << p.radius << " << " << minToolRadius << " at " << p.pos << endl;
                //outpath << p.pos.x() << ' ' << p.pos.y() << ' ' << p.radius << ' ' << p.tangent.x() << ' ' << p.tangent.y() << endl;
                totalLength += (q-p.pos).length();
                q = p.pos;
                if( insertPoint(p, clipperPath, & beadWidth) ) {
                }
            }
            if( closePath(clipperPath, & beadWidth) ) {
                fills.push_back(std::unique_ptr<IceSLInterface::IPath>(new IceSLInterface::Path()));
                fills.back()->createPath(clipperPath, true);
                fills.back()->setPathType(pathType);
                fills.back()->setTool(extruder_id_);
                // customize flow
                int flow_idx = fills.back()->addPerVertexAttribute("flow_multiplier");
                for( int i = 0; i < static_cast<int>(beadWidth.size()); ++i ) {
                    fills.back()->setPerVertexAttributeValue(flow_idx, i, 1.5 * beadWidth[i] / nozzle_diameter_);
                }
                //int cst_flow_idx = fills.back()->addPathAttribute("flow_multiplier");
                //fills.back()->setPathAttributeValue(cst_flow_idx, minBeadWidth_ / nozzle_diameter_);
            }
            clipperPath.clear();
            beadWidth.clear();
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

        mat.removeDestroyed();
    }

    // Fill the fallback_surface for remaining infill
    for( Component & comp : components ) {
        Sampling::sampleSmoothPaths(comp.innerSmoothPaths, pts, minToolRadius);
        for( const auto & samples : pts ) {
            if( samples.empty() ) continue;
            fallback_surface.emplace_back();
            ClipperLib::Path & clipperPath = fallback_surface.back();
            for( const auto & p : samples ) {
                if( insertPoint(p, clipperPath) ) {
                }
            }
            if( ! closePath(clipperPath) ) {
                fallback_surface.pop_back();
            }
        }
    }

    //if( verbose_ ) cerr << endl << nbSamples << " samples.\n";
    //cerr << "Total length = " << totalLength << endl;
    return true;
}

// -----------------------------------------------

bool VWCInfillerPlugin::addPluginSettings(IceSLInterface::EnumerableSettingsInterface & enumerable_settings)
{
    using IE = IceSLInterface::EnumerableSettingsInterface;
    const int numBrushes = enumerable_settings.getNumberOfBrushes();
    minBeadWidth_.resize(numBrushes, 0.3f);
    maxBeadWidth_.resize(numBrushes, 0.7f);
    numBeads_.resize(numBrushes, 4);
    std::unique_ptr<IE::SettingInterface> s;
    for( int i = 0; i < numBrushes; ++i ) {
        auto Zarathustra = [this, &enumerable_settings, i]() -> bool { // show setting only when the infiller is selected
            std::string infiller;
            enumerable_settings.getSettingByName("infill_type_" + std::to_string(i))->getValue(infiller);
            return infiller == this->name();
        };
        s = enumerable_settings.addSettingFromPlugin(
                & numBeads_[i],
                & minBeads_,
                & maxBeads_,
                "vwc_num_beads_" + std::to_string(i),
                "Number of contours",
                "Brush_" + std::to_string(i),
                "Specifies the number of beads used during contouring.Zero beads means as many as possible.",
                912, // rank to order multiple settings, lower appears before
                Zarathustra,
                IE::e_NoUnit,
                IE::e_NoFlag
                );
        if (s == nullptr) {
            return false;
        }
        s = enumerable_settings.addSettingFromPlugin(
                & minBeadWidth_[i],
                & minBeadWidth_min_,
                & minBeadWidth_max_,
                "vwc_bead_min_width_mm_" + std::to_string(i),
                "Minimal perimeter width (mm)",
                "Brush_" + std::to_string(i),
                "Specifies the minimal bead width in millimeters used during contouring.",
                913, // rank to order multiple settings, lower appears before
                Zarathustra,
                IE::e_MM,
                IE::e_NoFlag
            );
        if (s == nullptr) {
            return false;
        }
        s = enumerable_settings.addSettingFromPlugin(
                & maxBeadWidth_[i],
                & maxBeadWidth_min_,
                & maxBeadWidth_max_,
                "vwc_bead_max_width_mm_" + std::to_string(i),
                "Maximal perimeter width (mm)",
                "Brush_" + std::to_string(i),
                "Specifies the Maximal bead width in millimeters used during contouring.",
                914, // rank to order multiple settings, lower appears before
                Zarathustra,
                IE::e_MM,
                IE::e_NoFlag
            );
        if (s == nullptr) {
            return false;
        }
        s = enumerable_settings.addSettingFromPlugin(
                & overtakeIcesl_,
                & overtakeIcesl_min_,
                & overtakeIcesl_max_,
                "vwc_overtakes_icesl_" + std::to_string(i),
                "Produce perimeter and shells",
                "Brush_" + std::to_string(i),
                "The VWC plugin overrides IceSL's perimeter and shells and produces its own.",
                914, // rank to order multiple settings, lower appears before
                Zarathustra,
                IE::e_NoUnit,
                IE::e_NoFlag
                );
        if (s == nullptr) {
            return false;
        }
    }
    return true;
}

// -----------------------------------------------
