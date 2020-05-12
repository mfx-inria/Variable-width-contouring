/* vim: set sts=0 ts=4 sw=4 noet tw=0 : */

#include "fill_config.h"

#ifndef FILL_NO_CAIRO
#include <cairo/cairo-pdf.h>
#endif
#include <tclap/CmdLine.h>

#include <limits>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "vec.h"
#include "chronograph.h"
#include "VariableWidthContouring.h"

#include "MATGraph.h"
#ifndef FILL_NO_CAIRO
#include "CairoSurface.h"
#endif
#include "MATGraphBuilder.h"
#include "Sampling.h"

using namespace std;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - DEBUG

bool gDebug = false;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - COMMAND LINE

static TCLAP::CmdLine gCmdLine("fill, fills you with warmth and comfort.", ' ', "1.0-alpha");
#ifndef FILL_NO_CAIRO
// drawing options
static TCLAP::SwitchArg tightArg("t", "tight", "Don't draw tool bounds and use tight bounding box", false);
static TCLAP::SwitchArg drawGridArg("", "grid", "Draw grid", false);
static TCLAP::SwitchArg noDrawBisectorArg("", "no-bisector", "Do not draw tracks bisector", false);
static TCLAP::SwitchArg drawSubStepsArg("", "sub", "Draw substeps", false);
static TCLAP::SwitchArg noPDFArg("", "no-pdf", "Don't output PDF file", false);
static TCLAP::SwitchArg drawInputVerticesArg("", "draw-input-vertices", "Draw input vertices", false);
static TCLAP::SwitchArg drawMAVerticesArg("", "draw-ma-vertices", "Draw medial-axis vertices", false);
static TCLAP::ValueArg<int> resetArg("r", "reset", "Reset at given step", false /* required? */, -1, "non-negative number");
static TCLAP::ValueArg<double> lineWidthArg("l", "line-width-multiplier", "Multiply line width", false /* required? */, 1, "non-negative number");
static TCLAP::ValueArg<string> toolPosArg("", "tool-pos", "Tool position: n t b l r", false /* required? */, "br", "[n,t,b,l,r,bl,br,tl,tr]");
#endif
// process options
static TCLAP::ValueArg<int> stepsArg("s", "steps", "Number of steps", false /* required? */, 9999, "non-negative number");
static TCLAP::ValueArg<double> scalingArg("", "scaling", "Input scaling", false /* required? */,	1.0, "positive number");
static TCLAP::ValueArg<double> simplifyArg("", "simplify", "Simplify the medial axis", false /* required? */,	1.05, "float > 1");
static TCLAP::SwitchArg simplifyAllArg("", "simplify-all", "Simplify the full medial axis at the beginning", false);
static TCLAP::ValueArg<double> minOffsetArg("m", "min-offset", "Minimal offset value", false /* required? */,	0, "positive number");
static TCLAP::ValueArg<double> maxOffsetArg("M", "max-offset", "Maximal offset value", false /* required? */,	0, "positive number");
static TCLAP::ValueArg<double> maxCollapseFactorArg("C", "max-collapse-factor", "Multiplier applied to the minimal collapse radius to compute the maximal collapse radius", false /* required? */,	0, "must be > 1");
static TCLAP::ValueArg<double> samplingPeriodArg("", "sampling-period", "Maximal output segment length in mm (not strictly enforced)", false /* required? */, 0.005, "float >= 0.005 [unit = mm]");
static TCLAP::SwitchArg innerUnderfillArg("v", "inner-underfill", "Draw the inner underfill", false);
static TCLAP::SwitchArg outerUnderfillArg("", "outer-underfill", "Draw the outer underfill", false);
static TCLAP::SwitchArg sharpCutArg("k", "sharp-cut", "collapse cuts sharp at 4gamma", false);
static TCLAP::SwitchArg noShaveArg("", "no-shave", "Disable shaving", false);
static TCLAP::SwitchArg verboseArg("", "verbose", "Print info and progress", false);
static TCLAP::SwitchArg doSampleArg("", "do-sample-paths", "Do sample the print paths even if without output file", false);
// I/O arguments
static TCLAP::ValueArg<string> polygonInputFileArg("p", "polygon", "Input file for polygon", false /* required? */, "-", "path to file");
#ifndef FILL_NO_CAIRO
#define RIBBON_FILE_REQUIRED false
#else
#define RIBBON_FILE_REQUIRED true
#endif
static TCLAP::ValueArg<string> ribbonOutputFileArg("o", "output", "Output file for ribbons", RIBBON_FILE_REQUIRED, "", "path to file");

static double gMinOffset, gMaxOffset, gScale, gSimplifyThreshold;
static string gInputFile, gOutputFile;
static bool computeBisectors, gSubSteps, gTight, gGrid, gVerbose, gReallyReallyDoSample;
string gToolPos;
int gCurrentStep, gTotalSteps;
bool gInnerUnderfill, gOuterUnderfill, gUnderfill;

bool readCommandLine(int argc, char **argv) {
	try {
#ifndef FILL_NO_CAIRO
		gCmdLine.add(tightArg);
		gCmdLine.add(drawGridArg);
		gCmdLine.add(noDrawBisectorArg);
		gCmdLine.add(drawSubStepsArg);
		gCmdLine.add(noPDFArg);
		gCmdLine.add(drawInputVerticesArg);
		gCmdLine.add(drawMAVerticesArg);
		gCmdLine.add(resetArg);
		gCmdLine.add(lineWidthArg);
		gCmdLine.add(toolPosArg);
#endif
		gCmdLine.add(stepsArg);
		gCmdLine.add(scalingArg);
		gCmdLine.add(simplifyArg);
		gCmdLine.add(simplifyAllArg);
		gCmdLine.add(minOffsetArg);
		gCmdLine.add(maxOffsetArg);
		gCmdLine.add(maxCollapseFactorArg);
		gCmdLine.add(samplingPeriodArg);
		gCmdLine.add(innerUnderfillArg);
		gCmdLine.add(outerUnderfillArg);
		gCmdLine.add(sharpCutArg);
		gCmdLine.add(noShaveArg);
		gCmdLine.add(verboseArg);
		gCmdLine.add(doSampleArg);
		gCmdLine.add(polygonInputFileArg);
		gCmdLine.add(ribbonOutputFileArg);
		gCmdLine.parse(argc, argv);

#ifndef FILL_NO_CAIRO
		gTight = tightArg.getValue();
		gToolPos = toolPosArg.getValue();
		for( char c : gToolPos ) {
			if( (c == 'c') || (c == 'n') || (c == 'b') || (c == 't') || (c == 'l') || (c == 'r') ) continue;
			exit(-1);
		}
		if( gToolPos.empty() ) gToolPos = "br";
		gSubSteps = drawSubStepsArg.getValue() && (!noPDFArg.getValue());
		gGrid = drawGridArg.getValue();
		computeBisectors = ! noDrawBisectorArg.getValue();
#endif
		gVerbose = verboseArg.getValue();
		gReallyReallyDoSample = doSampleArg.getValue();
		double sp = samplingPeriodArg.getValue();
		Sampling::len_threshold = std::max(sp, 0.005); // in mm.
		gSimplifyThreshold = std::max(1.0, simplifyArg.getValue());
		if( minOffsetArg.getValue() < 0.0f )
			throw TCLAP::ArgException("bad value", "-m --min-offset");
		gMinOffset = minOffsetArg.getValue();
		if( maxOffsetArg.getValue() < 0.0f )
			throw TCLAP::ArgException("bad value", "-M --max-offset");
		gMaxOffset = maxOffsetArg.getValue();
		gScale = scalingArg.getValue();
		gInnerUnderfill = innerUnderfillArg.getValue();
		gOuterUnderfill = outerUnderfillArg.getValue();
		if( gInnerUnderfill && gOuterUnderfill )
			throw TCLAP::ArgException("conflict", "--inner-underfill --outer-underfill");
		gUnderfill = gOuterUnderfill || gInnerUnderfill;
		if( gScale <= 0.0 )
			throw TCLAP::ArgException("bad value", "-s --scaling");
		gInputFile = polygonInputFileArg.getValue();
		gOutputFile = ribbonOutputFileArg.getValue();
		if( gUnderfill ) {
			gTight = true;
			computeBisectors = false;
			gSubSteps = false;
			gToolPos = "n";
		}
		return false;
	}
	catch (const TCLAP::ArgException & e) {
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
	} catch (...) { // catch any exceptions
		cerr << "Error: unknown exception caught" << endl;
	}
	return true;
}


void readPolygons(istream & in, Paths & polys, BBox & bbox) {
	polys.clear();
	double b = numeric_limits<double>::max();
	bbox.min_ = Vec2d(+b,+b);
	bbox.max_ = Vec2d(-b,-b);
	double minO, maxO;
	in >> minO >> maxO;
	if( gMinOffset <= 0.0 ) gMinOffset = minO;
	if( gMaxOffset <= 0.0 ) gMaxOffset = maxO;
	while(in) {
		size_t numV;
		in >> numV;
		if( ! in ) break;
		polys.emplace_back();
		auto & path = polys.back();
		Path temp;
		for( size_t v(0); v < numV; ++v ) {
			double x, y;
			in >> x >> y;
			x *= gScale;
			y *= gScale;
			temp.emplace_back(x, y);
		}
		if( temp.size() <= 2 ) {
			polys.pop_back();
			continue;
		}
		int i(0), j(1), k(2);
		path.push_back(temp[0]);
		bbox.extendToCover(BBox(temp[0]));
		while( j > i ) {
			if( temp[i] == temp[j] ) { // duplicate vertices
				if( gVerbose ) cerr << "KILLED a duplicate vertex: bad input file" << endl;
				j = k;
				k = (k+1)%numV;
				continue;
			}
			Vec2d v0 = temp[j] - temp[i];
			Vec2d v1 = temp[k] - temp[j];
			if( ((v0 | v1) < (-0.999999)*v0.length()*v1.length()) &&
					(std::abs(det(v0, v1)) < 1e-6*v0.length()*v1.length()) ) {
				if( gVerbose ) cerr << "KILLED a U-turn: bad input file" << endl;
				j = k;
				k = (k+1)%numV;
				continue;
			}
			if( ((v0 | v1) > 0.999999*v0.length()*v1.length()) &&
					(std::abs(det(v0, v1)) < 1e-6*v0.length()*v1.length()) ) {
				if( gVerbose ) cerr << "KILLED a colinearity at " << temp[j] << ": bad input file" << endl;
				j = k;
				k = (k+1)%numV;
				continue;
			}
			path.push_back(temp[j]);
			bbox.extendToCover(BBox(temp[j]));
			i = j;
			j = k;
			k = (k+1)%numV;
		}
	}
	if (polys.empty())
    {
        std::cerr << "ERROR: Couldn't read polygons! (or empty polygons file)\n";
    }
}

double polygonArea(const Paths & paths, const Vec2d & center) {
	double area(0);
	for( auto path : paths ) {
		int n = path.size();
		int i = n-1;
		for( int j = 0; j < n; ++j ) {
			area += det(path[i]-center, path[j]-center);
			i = j;
		}
	}
	return area/2.0;
}

int main(int argc, char **argv) {

	if( readCommandLine(argc, argv) ) exit(EXIT_FAILURE);
	gTotalSteps = stepsArg.getValue();

	Paths inputPoly;
	BBox bbox({0,0});

	string outFileName = "example.pdf";

	if( gInputFile == "-" )
		readPolygons(cin, inputPoly, bbox);
	else {
		ifstream ifs(gInputFile.c_str());
        if (! ifs.good() )
        {
            cerr << "Cannot open \"" << gInputFile << "\" !\n";
            exit(-1);
        }
		readPolygons(ifs, inputPoly, bbox);
		outFileName = gInputFile + ".pdf";
	}

	bool output = gOutputFile != "";
	ofstream out;
	if( output ) {
		if( gVerbose ) cerr << "Output ribbons to " << gOutputFile << endl;
		out = ofstream(gOutputFile.c_str());
		out << std::fixed << std::setprecision(7);
	}

	Chronograph TotalTimer("Beads computation");
	TotalTimer.start();

	double minToolRadius = gMinOffset / 2.0;
	double maxToolRadius = gMaxOffset / 2.0;
	if( gVerbose ) cerr << "Tool radius ∈ [" << minToolRadius << ", " << maxToolRadius << ']';
	if( gVerbose ) cerr << ". Bounding box = " << bbox.min_ << " -> " << bbox.max_;
	if( gVerbose ) cerr << " (" << bbox.size().x() << " x " << bbox.size().y() << " mm).\n";
	if( gVerbose ) cerr << "(Loose) sampling period set to " << Sampling::len_threshold << " mm." << endl;
	if( gUnderfill ) {
		cout << setprecision(10) << std::fixed <<
			(bbox.size().y()*bbox.size().x())
			<< ' '
			<< polygonArea(inputPoly, bbox.center());
	}

#ifndef FILL_NO_CAIRO
	CairoSurface * pdf, * contours;
	if( ! noPDFArg.getValue() ) {
		// The output PDF
		pdf = new CairoSurface(bbox, outFileName.c_str(), minToolRadius, maxToolRadius, gTight);
		pdf->setDrawMAVertices(drawMAVerticesArg.getValue());
		pdf->setDrawGrid(gGrid);
		pdf->setLineWidthMultiplier(lineWidthArg.getValue());
		pdf->setUnderfill(gUnderfill);
		// The surface that accumulates all the tracks/beads/ribbons and the background input polygon
		contours = new CairoSurface(bbox, "", minToolRadius, maxToolRadius, gTight, true);
		contours->setDrawVertices(drawInputVerticesArg.getValue());
		contours->setLineWidthMultiplier(lineWidthArg.getValue());
		contours->setUnderfill(gUnderfill);
		contours->drawPolygons(inputPoly);
	}
#endif

	std::vector<Disk> cheekPrecursors;

	// INITIALIZATION

	MAOffsetInfo maoi[2];
	vector<vector<Sample>> pts;

	using Components = MATGraph::Components;
	using Component  = MATGraph::Component;
	Components components;

	MATGraph mat;
	Chronograph MATTimer("MAT building");
	MATTimer.start();
	buildMATGraphWithCGAL(mat, inputPoly);
	buildMATGraphWithBOOST(mat, inputPoly);
	MATTimer.stop();
	mat.debugCheck();

	VariableWidthContouring sdg(mat, minToolRadius, maxToolRadius);
	sdg.setMaxCollapseFactor(maxCollapseFactorArg.getValue());
	sdg.setSimplifyThreshold(gSimplifyThreshold);
	sdg.setWorkingMAOffsetInfo(& maoi[0]);

#ifndef FILL_NO_CAIRO
	// DRAW FIRST PAGE : input & its medial axis
	if( ! noPDFArg.getValue() ) {
		CairoSurface::fuse(bbox, *contours, *pdf, mat, components, pts,  minToolRadius,
				maxToolRadius, computeBisectors);
		contours->step();
	}
#endif

	// FILTER AND SIMPLIFY THE INPUT

	components.clear();
	mat.computeConnectedComponents(components);
	sdg.filterSmallFeatures(minToolRadius * 2, components);
	mat.debugCheck();

	sdg.simplifyMACombinatorics();
	if( simplifyAllArg.getValue() )
		sdg.simplifyAllMAGeometry();

	size_t nbSamples(0); // total number of samples generated along all print paths

	bool notYetReset = true;

	double totalLength = 0.0;

	// GO!

	for( int i(1); i <= gTotalSteps+1; ++i ) {

		gCurrentStep = i;

		// Compute connected components. We work on each separately, so that the uniform offset can be customized to each.
		components.clear();
		mat.computeConnectedComponents(components);
		if( 0 == components.size() ) {
			gTotalSteps = i-1;
		}
		else {
			//cerr << "Shrink step " << i << " on " << components.size() << " components\n";
			if( gVerbose ) cerr << ((i%10) ? '.':'|');
		}
		// Update maoi and the smooth contour
		maoi[0].clear();
		for( Component & comp : components ) {
			sdg.computeSmoothPaths(comp, comp.outerSmoothPaths); // also fills maoi[0]
		}

#ifndef FILL_NO_CAIRO
		if( gOuterUnderfill ) {
			pdf->polygonPath(inputPoly);
			for( Component & comp : components ) {
				pdf->drawSmoothPaths(comp.outerSmoothPaths);
			}
			cairo_set_source_rgb(pdf->context_, 0., 0., 0.);
			cairo_set_fill_rule(pdf->context_, CAIRO_FILL_RULE_EVEN_ODD);
			cairo_fill(pdf->context_);
			break;
		}

		if( i == resetArg.getValue() && notYetReset ) { // FUN HACK
			notYetReset = false;
			bbox = sdg.computeBBox();
			//cerr << "New bounding box = " << bbox.min_ << " -> " << bbox.max_ << endl;
			if( ! noPDFArg.getValue() ) {
				delete pdf;
				delete contours;
				pdf = new CairoSurface(bbox, outFileName.c_str(), minToolRadius, maxToolRadius, gTight);
				contours = new CairoSurface(bbox, "", minToolRadius, maxToolRadius, gTight, true);
				contours->setDrawVertices(drawInputVerticesArg.getValue());
				contours->setLineWidthMultiplier(lineWidthArg.getValue());
				contours->setUnderfill(gUnderfill);
				pdf->setDrawMAVertices(drawMAVerticesArg.getValue());
				pdf->setDrawGrid(gGrid);
				pdf->setUnderfill(gUnderfill);
				pdf->setLineWidthMultiplier(lineWidthArg.getValue());
				cairo_t * c = contours->context_;
				cairo_save(c);
				cairo_set_source_rgb(c, 0.8, 0.8, 0.8);
				for( Component & comp : components ) {
					contours->drawSmoothPaths(comp.outerSmoothPaths);
				}
				cairo_set_fill_rule(c, CAIRO_FILL_RULE_EVEN_ODD);
				cairo_fill(c);
				cairo_restore(c);
				cairo_t * cr = cairo_create(pdf->surface_);
				cairo_set_source_surface(cr, contours->surface_, 0, 0);
				cairo_paint(cr);
				cairo_destroy(cr);
				pdf->drawMA(mat);
				cairo_show_page(pdf->context_);
			}
			//steps = steps - (i-1);
			i = 1;
		}
		if( (! noPDFArg.getValue()) && gSubSteps && (i > 1) ) pdf->step();
#endif

		auto drawSubStep = [&](bool drawInner = false, bool noShave = true, bool noMA = false) {
#ifndef FILL_NO_CAIRO
			if( ! gSubSteps ) return;
			if( noPDFArg.getValue() ) return;
			if( gUnderfill ) return;
			if( i == gTotalSteps+1 ) return;
			for( Component & comp : components ) {
				pdf->drawContour(comp);
			}
			if( drawInner ) {
				if( noShave ) maoi[1].clear();
				for( Component & comp : components ) {
					pdf->step();
					if( noShave ) {
						SmoothPaths paths;
						sdg.computeSmoothPaths(comp, paths, /*walk_on_fire*/false, &maoi[1]);
						pdf->drawContour(paths);
					} else { // we want to draw Shaved edges; we assume innerSmoothPaths has been computed
						// because we can't touch maoi[1] anymore
						pdf->drawContour(comp.innerSmoothPaths);
					}
					pdf->step();
				}
			}
			if( ! noMA ) {
				pdf->drawMA(mat, noShave);
			}
			//pdf->step();
			pdf->drawPaths(pts); // draw bisector
			//pdf->step();
			pdf->drawDisks(cheekPrecursors);
			if( gGrid ) {
				pdf->drawGrid(1.0, 0.1, true);
				pdf->drawGrid(10.0, 0.3);
			}
			cairo_set_source_rgb(pdf->context_, 0, 0, 0);
			pdf->drawTools(bbox);
			cairo_show_page(pdf->context_);
#endif
		};

		drawSubStep();
		//sdg.shrinkMA(components, cheekPrecursors); // works on maoi[0] and updates it when an MA edge is subdivided
		// 1. TRIMMING
		sdg.shrinkMA_trim(components); // works on maoi[0] and updates it when an MA edge is subdivided
		drawSubStep(true);

		// 2. COLLAPSE
		// works on maoi[0] and updates it when an MA edge is subdivided
		bool hadSomeCollapse = sdg.shrinkMA_collapse(components, cheekPrecursors, sharpCutArg.getValue());
		if( hadSomeCollapse ) {
			drawSubStep(true);

			// 3. CHEEKS
			if( ! cheekPrecursors.empty() ) {
				if( gVerbose ) std::cerr << "MODELING CHEEKS..." << endl;
				sdg.shrinkMA_cheeks(components); // works on maoi[0] and updates it when an MA edge is subdivided
				drawSubStep(true);
			}
		}

		// 4. GLOBAL OFFSET
		sdg.shrinkMA_globalOffset(components); // works on maoi[0] and updates it when an MA edge is subdivided
		sdg.simplifyMACombinatorics(&components);
		drawSubStep(true);

		maoi[1] = maoi[0];  // Save the |maoi| of the outer contour

		// 5. CLIPPING & SHAVE
		bool someEdgesWereEradicated = sdg.simplifyCollapsedMAGeometry(&components);
		sdg.clipAndShaveCollapsedParts(components, noShaveArg.getValue());
		if( someEdgesWereEradicated ) {
			//cerr << "(some edges eradicated)";
			maoi[1].clear();
			for( Component & comp : components ) {
				sdg.computeSmoothPaths(comp, comp.simpleOuterSmoothPaths, /*walk_on_fire*/true, &maoi[1]);
			}
		}
		sdg.simplifyMACombinatorics(&components);

		// Re-compute smoothPaths and maoi for inner contour
		maoi[0].clear();
		for( Component & comp : components ) {
			sdg.computeSmoothPaths(comp, comp.innerSmoothPaths); // also fills maoi[0]
		}
		if( hadSomeCollapse ) {
			drawSubStep(true, false);
		}

#ifndef FILL_NO_CAIRO
		if( gReallyReallyDoSample || output || ( ! noPDFArg.getValue()) ) {
#else
		{
#endif
			pts.clear();
			for( Component & comp : components ) {
				sdg.samplePrintPath(comp, & maoi[1], pts);
			}
		}
		if( ! hadSomeCollapse ) {
			for( Component & comp : components ) {
				comp.voids.clear();
			}
		}
		drawSubStep(true, false);
#if 0 // ONLY FOR GENERATING THE BASIC COLLAPSE FIGURE!
		drawSubStep(true, false, true);
#endif
		cheekPrecursors.clear();

		// output
		for( const auto & samples : pts ) {
			nbSamples += samples.size();
			if( (gReallyReallyDoSample || output) && ! samples.empty() ) {
				out << "closed " << samples.size() << " " << (i-1) << endl;
				Vec2d q = samples.back().pos;
				for( const auto & p : samples ) {
					if( p.radius > maxToolRadius * 1.01 )
						cerr << "ERROR : output extrusion radius is too large: "
							<< p.radius << " >> " << maxToolRadius << " at " << p.pos << endl;
					if( p.radius < minToolRadius / 1.01 )
						cerr << "ERROR : output extrusion radius is too small: "
							<< p.radius << " << " << minToolRadius << " at " << p.pos << endl;
					out << p.pos.x() << ' ' << p.pos.y()
						<< ' ' << p.radius
						<< ' ' << p.tangent.x() << ' ' << p.tangent.y()
						<< endl;
					totalLength += (q-p.pos).length();
					q = p.pos;
				}
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
		}

#ifndef FILL_NO_CAIRO
		// Draw new page in PDF output
		if( ! noPDFArg.getValue() ) {
			CairoSurface::fuse(bbox, *contours, *pdf, mat, components, pts, minToolRadius,
					maxToolRadius, computeBisectors);
		}
#endif
		pts.clear();

		mat.removeDestroyed();
	}
	if( gVerbose ) cerr << endl << nbSamples << " samples.\n";

	if( output ) {
		out.close();
	}

#ifndef FILL_NO_CAIRO
	if( ! noPDFArg.getValue() ) {
		delete pdf;
		delete contours;
	}
#endif

	TotalTimer.stop();
	if( gUnderfill ) {
		cout << setprecision(10) << std::fixed << ' ' << totalLength << endl;
	} else {
		TotalTimer.print(cout);
		cout << " (" << (100.0-100.0*MATTimer.elapsed_time()/TotalTimer.elapsed_time()) << " % without MAT)" << endl;
		//MATTimer.print(cout);
		//cout << " (" << (100.0*MATTimer.elapsed_time()/TotalTimer.elapsed_time()) << " %)" << endl;

	}

	return EXIT_SUCCESS;
}
