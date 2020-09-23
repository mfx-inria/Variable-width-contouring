/* vim: set sts=0 ts=4 sw=4 noet tw=0 : */

#include "VariableWidthContouringForIceSL.h"

#include <iostream>

void buildMATGraphWithBOOST(MATGraph& mat, const CLPaths & inputPoly);

using namespace std;

VariableWidthContouringForIceSL::
VariableWidthContouringForIceSL() {
	debug_ = false;
	verbose_ = true;

	minBeadWidth_ = 0.3; // mm
	maxBeadWidth_ = 1.0; // mm
	simplificationThreshold_ = 1.05; // unit-less

	Sampling::len_threshold = 0.01; // mm

	if( verbose_ ) {
		cerr << "Variable-width Contouring bead width ∈ [" << minBeadWidth_ << ", " << maxBeadWidth_ << ']'
			<< ". (Loose) sampling period set to " << Sampling::len_threshold << " mm." << endl;
	}
}

void
VariableWidthContouringForIceSL::
go(const CLPaths & input, CLPaths & output, int numberOfBeads) {

	std::vector<Disk> cheekPrecursors;

	// INITIALIZATION

	vector<vector<Sample>> pts;

	using Components = MATGraph::Components;
	using Component  = MATGraph::Component;
	Components components;

	MATGraph mat;
	buildMATGraphWithBOOST(mat, input);
	//mat.print();
	mat.debugCheck();

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

	// GO!

	for( int step(1); step <= numberOfBeads; ++step ) {

		// Compute connected components. We work on each separately, so that the uniform offset can be customized to each.
		components.clear();
		mat.computeConnectedComponents(components);
		if( 0 == components.size() ) {
			break;
		}
		else {
			//cerr << "Shrink step " << step << " on " << components.size() << " components\n";
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
		for( const auto & samples : pts ) {
			nbSamples += samples.size();
			if( ! samples.empty() ) {
				output.emplace_back();
				auto & outpath = output.back();
				break;
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

		pts.clear();

		mat.removeDestroyed();
	}
	if( verbose_ ) cerr << endl << nbSamples << " samples.\n";
}
