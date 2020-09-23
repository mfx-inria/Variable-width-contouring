/* vim: set sts=0 ts=4 sw=4 noet tw=0 : */

#include "fill_config.h"

#ifndef FILL_NO_CAIRO
#include <cairo/cairo-pdf.h>
#endif
#include <tclap/CmdLine.h>

#include <fstream>
#include <iostream>
//#include <iomanip>

//#include "vec.h"
#include "chronograph.h"
#include "VariableWidthContouring.h"

#ifndef FILL_NO_CAIRO
#include "CairoSurface.h"
#endif
#include "MATGraphBuilder.h"

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
static TCLAP::SwitchArg innerUnderfillArg("v", "inner-underfill", "Draw the inner underfill", false);
static TCLAP::SwitchArg outerUnderfillArg("", "outer-underfill", "Draw the outer underfill", false);
#endif
// process options
#if (FILL_MA_BOOST && FILL_MA_CGAL)
static TCLAP::SwitchArg boostArg("b", "boost", "Use BOOST instead of CGAL", false);
#endif
static TCLAP::ValueArg<int> stepsArg("s", "steps", "Number of steps", false /* required? */, 9999, "non-negative number");
static TCLAP::ValueArg<double> scalingArg("", "scaling", "Input scaling", false /* required? */,	1.0, "positive number");
static TCLAP::ValueArg<double> simplifyArg("", "simplify", "Simplify the medial axis", false /* required? */,	1.05, "float > 1");
static TCLAP::SwitchArg simplifyAllArg("", "simplify-all", "Simplify the full medial axis at the beginning", false);
static TCLAP::ValueArg<double> minOffsetArg("m", "min-offset", "Minimal offset value", false /* required? */,	0, "positive number");
static TCLAP::ValueArg<double> maxOffsetArg("M", "max-offset", "Maximal offset value", false /* required? */,	0, "positive number");
static TCLAP::ValueArg<double> samplingPeriodArg("", "sampling-period", "Maximal output segment length in mm (not strictly enforced)", false /* required? */, 0.005, "float >= 0.005 [unit = mm]");
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

// Processing parameters
static double gMinOffset, gMaxOffset, gScale, gSimplifyThreshold;
static bool gUseBoost;
// I/O parameters
static string gInputFile, gOutputFile;
// Drawing/benchmarking parameters
static bool gDrawPrintPaths, gSubSteps, gTight, gGrid, gVerbose, gReallyReallyDoSample;
string gToolPos;
int gTotalSteps;
#ifndef FILL_NO_CAIRO
bool gInnerUnderfill, gOuterUnderfill, gUnderfill;
#endif

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
		gCmdLine.add(innerUnderfillArg);
		gCmdLine.add(outerUnderfillArg);
#endif
#if (FILL_MA_BOOST && FILL_MA_CGAL)
		gCmdLine.add(boostArg);
#endif
		gCmdLine.add(stepsArg);
		gCmdLine.add(scalingArg);
		gCmdLine.add(simplifyArg);
		gCmdLine.add(simplifyAllArg);
		gCmdLine.add(minOffsetArg);
		gCmdLine.add(maxOffsetArg);
		gCmdLine.add(samplingPeriodArg);
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
			cerr << "The tool-pos arguments is invalid.\n";
			exit(-1);
		}
		if( gToolPos.empty() ) gToolPos = "br";
		gSubSteps = drawSubStepsArg.getValue() && (!noPDFArg.getValue());
		gGrid = drawGridArg.getValue();
		gDrawPrintPaths = ! noDrawBisectorArg.getValue();
		gInnerUnderfill = innerUnderfillArg.getValue();
		gOuterUnderfill = outerUnderfillArg.getValue();
		gUnderfill = gOuterUnderfill || gInnerUnderfill;
		if( gInnerUnderfill && gOuterUnderfill )
			throw TCLAP::ArgException("conflict", "--inner-underfill --outer-underfill");
		if( gUnderfill ) {
			gTight = true;
			gDrawPrintPaths = false;
			gSubSteps = false;
			gToolPos = "n";
		}
#endif
#if (FILL_MA_BOOST && FILL_MA_CGAL)
		gUseBoost = boostArg.getValue();
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
		if( gScale <= 0.0 )
			throw TCLAP::ArgException("bad value", "-s --scaling");
		gInputFile = polygonInputFileArg.getValue();
		gOutputFile = ribbonOutputFileArg.getValue();
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

#ifndef FILL_NO_CAIRO
	if( gUnderfill ) {
		cout << setprecision(10) << std::fixed <<
			(bbox.size().y()*bbox.size().x())
			<< ' '
			<< polygonArea(inputPoly, bbox.center());
	}
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

	std::vector<Disk> cheekPrecursors; // not used anymore. not in the paper. Left here as we might come back to it later (low prob.)

	// INITIALIZATION

	BoundaryCircles maoi[2]; // stores boundary circles (supports of the circular
	// arcs along the shape boundary) for two consecutive round shapes. (The bead
	//is roughly the  bisecotr between the two shapes.
	vector<vector<Sample>> pts; // stores the sampled print paths.

	using Components = MATGraph::Components;
	using Component  = MATGraph::Component;
	Components components;
	// stores connected components of the MA (one set<MAvert*> per component),
	// together with per-component data (some smooth paths and, importantly,
	// the "uniform offset".

	MATGraph mat; // The medial axis. THE MEDIAL AXIS!!
	Chronograph MATTimer("MAT building");
	MATTimer.start();
#if (FILL_MA_BOOST && FILL_MA_CGAL)
	if( gUseBoost )
		buildMATGraphWithBOOST(mat, inputPoly);
	else
		buildMATGraphWithCGAL(mat, inputPoly);
#elif FILL_MA_BOOST
		buildMATGraphWithBOOST(mat, inputPoly);
#else
		buildMATGraphWithCGAL(mat, inputPoly);
#endif
	//mat.print();
	MATTimer.stop();
	mat.debugCheck();

	// instantiate the class that implements the algorithm
	VariableWidthContouring vwc(mat, minToolRadius, maxToolRadius, gSimplifyThreshold);
	vwc.setWorkingBoundaryCircles(& maoi[0]);

#ifndef FILL_NO_CAIRO
	// DRAW FIRST PAGE : input & its medial axis
	if( ! noPDFArg.getValue() ) {
		CairoSurface::fuse(bbox, *contours, *pdf, mat, components, pts,  minToolRadius,
				maxToolRadius, gDrawPrintPaths, 0);
		contours->step(); // For bead color alternation
	}
#endif

	// FILTER AND SIMPLIFY THE INPUT

	components.clear();
	mat.computeConnectedComponents(components);
	vwc.filterSmallFeatures(minToolRadius * 2, components); // makes the shape 2gamma-fat
	mat.debugCheck();

	// Processing of the MA : (filtering, trimming, collapsing) may introduce
	// unecessary vertices on the MA graph: degree-2 vertices where nothing is
	// different on either side. These vertices can be safely removed with simplifyMACombinatorics()
	vwc.simplifyMACombinatorics();
	// MA simplification as described in the paper, but on the full MA.
	// Useful for comparison and fun, but not useful in practice.
	if( simplifyAllArg.getValue() ) vwc.simplifyAllMAGeometry();

	size_t nbSamples(0); // total number of samples generated along all print paths

	// To simulate "round" input shapes (as opposed to polygonal ones), we can
	// "reset" the processing by forgetting the first few iterations. This flag
	// is used to indicate that no reset has happened yet. See |resetArg|.
	bool notYetReset = true;

	double totalLength = 0.0; // total length of all the output beads.

	// GO!

	for( int stepNum(1); stepNum <= gTotalSteps+1; ++stepNum ) { // we do one more steps in order to output a nice final page result, without MA in the PDF.

		// Compute connected components. We work on each separately, so that the uniform offset can be customized to each.
		components.clear();
		mat.computeConnectedComponents(components);
		if( 0 == components.size() ) { // there is nothing left to do.
			gTotalSteps = stepNum-1;
		}
		else { // the usual case
			if( gVerbose ) cerr << ((stepNum%10) ? '.':'|');
		}
		// Compute the boundary circles in maoi[0] and the smooth contours
		maoi[0].clear();
		for( Component & comp : components ) {
			vwc.computeSmoothPaths(comp, comp.outerSmoothPaths); // also fills maoi[0]
		}

#ifndef FILL_NO_CAIRO
		if( gOuterUnderfill ) {
			// Outer underfill is the underfill due to the difference netween
			// the polygonal input shape and its 2-gamma-fat version over which
			// we actually work. Here, we want to estimate this outer
			// underfill, so we draw a black and white PDF page, and exit. We
			// alter use image processing to count pixels and estimate the
			// area.
			pdf->polygonPath(inputPoly);
			for( Component & comp : components ) {
				pdf->drawSmoothPaths(comp.outerSmoothPaths);
			}
			cairo_set_source_rgb(pdf->context_, 0., 0., 0.);
			cairo_set_fill_rule(pdf->context_, CAIRO_FILL_RULE_EVEN_ODD);
			cairo_fill(pdf->context_);
			break;
		}

		// If we "reset" the processing, forgetting all that was done before
		// now, then we have a bunch of data to manipulate...
		if( stepNum == resetArg.getValue() && notYetReset ) { // FUN HACK
			notYetReset = false; // To avoid an infinite loop of resets
			bbox = vwc.computeBBox();
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
			//steps = steps - (stepNum-1);
			stepNum = 1;
		}
		if( (! noPDFArg.getValue()) && gSubSteps && (stepNum > 1) ) pdf->step();
#endif

		// A function to draw the current state of affairs in a single page, so
		// we can see the various steps of creating a single bead in several
		// pages of the pDF.
		auto drawSubStep = [&](bool drawInner = false, bool noShave = true, bool noMA = false) {
#ifndef FILL_NO_CAIRO
			if( ! gSubSteps ) return;
			if( noPDFArg.getValue() ) return;
			if( gUnderfill ) return;
			if( stepNum == gTotalSteps+1 ) return;
			for( Component & comp : components ) {
				pdf->drawContour(comp);
			}
			if( drawInner ) {
				if( noShave ) maoi[1].clear();
				for( Component & comp : components ) {
					pdf->step();
					if( noShave ) {
						SmoothPaths paths;
						vwc.computeSmoothPaths(comp, paths, /*walk_on_fire*/false, &maoi[1]);
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
		// shrinkMA() does all the steps at once. We don't use it anymore,
		// preferring calling function for each substep, so that we can draw
		// each substep with drawSubStep(). So it may be outdated.
		//vwc.shrinkMA(components, cheekPrecursors); // works on maoi[0] and updates it when an MA edge is subdivided

		// 1. TRIMMING
		vwc.shrinkMA_trim(components); // works on maoi[0] and updates it when an MA edge is subdivided
		drawSubStep(true);

		// 2. COLLAPSE
		// works on maoi[0] and updates it when an MA edge is subdivided
		bool hadSomeCollapse = vwc.shrinkMA_collapse(components, cheekPrecursors, sharpCutArg.getValue());
		if( hadSomeCollapse ) {
			drawSubStep(true);

			// 3. CHEEKS // NOT USE ANYMORE |cheekPrecursors| is always empty (could be changed in VariableWidthContouringLabeling.cpp)
			if( ! cheekPrecursors.empty() ) {
				if( gVerbose ) std::cerr << "MODELING CHEEKS..." << endl;
				vwc.shrinkMA_cheeks(components); // works on maoi[0] and updates it when an MA edge is subdivided
				drawSubStep(true);
			}
		}

		// 4. GLOBAL PARALLEL OFFSET
		// Remember that each connected component might have a different
		// uniform parallel offset. See Section 5.4 of the paper.
		vwc.shrinkMA_globalOffset(components); // works on maoi[0] and updates it when an MA edge is subdivided
		vwc.simplifyMACombinatorics(&components); // Keep it clean and tidy
		drawSubStep(true);

		maoi[1] = maoi[0];  // Save the |maoi| of the outer contour

		// 5. CLIPPING & SHAVE
		// Simplify the input shape around the collapsed axis
		bool someEdgesWereEradicated = vwc.simplifyCollapsedMAGeometry(&components);
		// Then shave the collapsed axis. Clipping refers to the points where
		// the MA goes outside the inner shape and inside the outer shape,
		// passing through the boundary of a maximal disk centered at the root
		// of a maximal trimmed sub-tree.
		vwc.clipAndShaveCollapsedParts(components, noShaveArg.getValue());
		if( someEdgesWereEradicated ) {
			// If some collapsed arc were eradicated (simplified), then the
			// INPUT shape is changed, so we must compute the smooth contour
			// (smooth path) of the simplified input shape:
			maoi[1].clear();
			for( Component & comp : components ) {
				vwc.computeSmoothPaths(comp, comp.simpleOuterSmoothPaths, /*walk_on_fire*/true, &maoi[1]);
			}
		}
		vwc.simplifyMACombinatorics(&components); // Keep it clean and tidy

		// Re-compute smoothPaths and maoi for inner contour. The outer contour
		// / boundaryDisks is now in maoi[1]
		maoi[0].clear();
		for( Component & comp : components ) {
			vwc.computeSmoothPaths(comp, comp.innerSmoothPaths); // also fills maoi[0]
		}
		if( hadSomeCollapse ) {
			drawSubStep(true, false);
		}

#ifndef FILL_NO_CAIRO
		if( gReallyReallyDoSample || output || ( ! noPDFArg.getValue()) ) {
#else
		{
#endif
			// If required, we now sample the print paths.
			pts.clear();
			for( Component & comp : components ) {
				vwc.samplePrintPath(comp, & maoi[1], pts);
			}
		}
		if( ! hadSomeCollapse ) {
			// Voids are the outline of the underfill, for drawing them
			// red/pink in the PDF. (!hadSomeCollapse) means that comp.voids
			// still stores voids from previous steps... To avoid messing our
			// brain, we clear them.
			for( Component & comp : components ) {
				comp.voids.clear();
			}
		}
		drawSubStep(true, false);
#if 0 // ONLY FOR GENERATING THE BASIC COLLAPSE FIGURE!
		drawSubStep(true, false, true);
#endif
		cheekPrecursors.clear();

		// Now we output the sampled print paths. Nothing fancy. Along the way,
		// we check that each sample is within the prescribed bounds, with some
		// tolerance of ≈ 1 %.
		for( const auto & samples : pts ) {
			nbSamples += samples.size();
			if( (gReallyReallyDoSample || output) && ! samples.empty() ) {
				out << "closed " << samples.size() << " " << (stepNum-1) << endl;
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
#if 0
				// Check angle formed by 3 consecutive samples. (That code was
				// there to detect and fix a bug.) Might still be useful in the
				// future.)
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
		// Draw new page in PDF output, showing the resulting bead and labeled medial-axis.
		if( ! noPDFArg.getValue() ) {
			CairoSurface::fuse(bbox, *contours, *pdf, mat, components, pts, minToolRadius,
					maxToolRadius, gDrawPrintPaths, stepNum);
		}
#endif
		// cleared the sampled print-paths, making room for the next iteration.
		pts.clear();

		// finally we cleanup the mess : remove all part of the MA that is non-Normal.
		mat.removeDestroyed();
	}// END OF THE MAIN LOOP

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
#ifndef FILL_NO_CAIRO
	if( gUnderfill ) {
#else
	if( false ) {
#endif
		cout << setprecision(10) << std::fixed << ' ' << totalLength << endl;
	} else {
		TotalTimer.print(cout);
		cout << " (" << (100.0-100.0*MATTimer.elapsed_time()/TotalTimer.elapsed_time()) << " % without MAT)" << endl;
		//MATTimer.print(cout);
		//cout << " (" << (100.0*MATTimer.elapsed_time()/TotalTimer.elapsed_time()) << " %)" << endl;

	}

	return EXIT_SUCCESS;
}
