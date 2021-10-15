/* vim: set sts=0 ts=4 sw=4 noet : */

#ifndef _CAVITIES_SAMPLING_H_
#define _CAVITIES_SAMPLING_H_

#include "vec.h"

#include <list>

#include "utils.h"

using namespace std;

/*
 */
struct AxisPos : Vec2i {
	AxisPos() : Vec2i(0, -1) {}
	int axis() const { return x(); }
	int & axis() { return x(); }
	int pos() const { return y(); }
	int & pos() { return y(); }
};

// This struct stores a sample point along a curve. The curve is either the boundary of the shape
// (input polygon, or offset versions of it), OR a print-path.
// When the curve is a print-path, the radius indicates the radius of deposition/extrusion; this is
// the local width of that ribbon (radius of maximal disk in the ribbon, centered at |pos|
// The tangent is useful for adaptive sampling.
struct Sample {
	Vec2d pos, tangent;
	double radius;
	AxisPos axisPos; // for computing the track contour along collapsed axes
};

using SampleFunction = std::function<Sample (const Sample &)>;

// insert linear approximation of the quadratic curve, except for the first point.
void diceQuadratic(const Vec2d & focus, const QuadraticArc & q, std::vector<Sample> & path);

class Sampling {

public:

	static constexpr double cos_5_degree = 0.99619469809;
	static constexpr double cos_4_degree = 0.99756405026;
	static constexpr double cos_3_degree = 0.99862953475;
	static constexpr double cos_2_degree = 0.99939082702;
	static constexpr double cos_1_degree = 0.99984769516;
	static constexpr double ang_threshold = cos_3_degree;
	static double len_threshold;

public:

	double min_radius, max_radius;
	Sampling(double min_radius, double max_radius)
	: min_radius(min_radius)
	, max_radius(max_radius)
	{}

	static Sample moveTowardBisectorWithDisk(const Sample & sample, double delta, const Disk & disk,
			bool * bad = nullptr);
	static Sample moveTowardBisectorWithSegment(const Sample & sample, const Vec2d &
			s0, const Vec2d & s1, bool & bad);
	static bool is_good_sample(const Sample & s0, const Sample & s1);
	static Vec2d splitArcInHalf(const Vec2d & v0, const Vec2d & v1, bool clockwise);

	static void sampleSmoothPaths(const SmoothPaths & smoothPaths, std::vector<std::vector<Sample>> & samples, double minToolRadius);

	// Used by sampleSmoothPath
	template< typename F >
	static void sampleCircularArc(const Vec2d & v0, const Vec2d & v1, const BoundaryCircle & circ,
			const F & out, bool clockwise) {
		struct VP {
			Vec2d v;
			Sample projected; // |border| is on the border of the variable-width ribbon; |projected| is on the ribbon medial axis or bisector.
		};
		auto makeVP = [&](const Vec2d & v) {
			Vec2d vn = v.normalized();
			Vec2d tangent = vn;
			if( clockwise ) tangent.rotateCW(); else tangent.rotateCCW();
			Sample border{circ.center() + (circ.radius() * vn), tangent, 0.0};
			return VP{vn, border};
		};
		if( circ.radius() <= 0.0 ) {
			Sample s{circ.center(), Vec2d(0,0), 0.0};
			out(s);
			return;
		}
		const double clockFactor = clockwise ? -1.0 : 1.0;
		list<VP> sampling;
		sampling.push_back(makeVP(v0));
		sampling.push_back(makeVP(v1));
		typename list<VP>::iterator it1 = sampling.begin(), it0 = it1++, last;
		last = --sampling.end();
		while( it1 != sampling.end() ) {
			if( is_good_sample(it0->projected, it1->projected) ) {
				out(it1->projected);
				it0 = it1++;
			} else { // splitting
				auto it2 = it1; ++it2;
				if( (it1 != last) && ((it0->v | it2->v) > 0.0) && (clockFactor*det(it0->v,it2->v) > 0.0) ) {
					it1 = sampling.erase(it1); // it1 is now equal to it2
					Vec2d d = it1->v - it0->v;
					it1 = sampling.insert(it1, makeVP(it0->v + (0.66666 * d)));
					it1 = sampling.insert(it1, makeVP(it0->v + (0.33333 * d)));
				} else { // split in 2
					it1 = sampling.insert(it1, makeVP(splitArcInHalf(it0->v, it1->v, clockwise)));
				}
			}
		}
	}

	template< typename F >
	static void sampleCircularArc(const Vec2d & v0, const Vec2d & v1, const BoundaryCircle & circ,
		const F & out, bool clockwise, SampleFunction & project, double delta) {
		struct VP {
			Vec2d v;
			Sample projected; // |border| is on the border of the variable-width ribbon; |projected| is on the ribbon medial axis or bisector.
		};
		auto makeVP = [&](const Vec2d & v) {
			Vec2d vn = v.normalized();
			Vec2d tangent = vn;
			if( clockwise ) tangent.rotateCW(); else tangent.rotateCCW();
			Sample proj, border{circ.center() + (circ.radius() * vn), tangent, delta}; // {position, tangent, radius}
			proj =  project(border);
			return VP{vn, proj};
		};
		const double clockFactor = clockwise ? -1.0 : 1.0;
		list<VP> sampling;
		sampling.push_back(makeVP(v0));
		sampling.push_back(makeVP(v1));
		typename list<VP>::iterator it1 = sampling.begin(), it0 = it1++, last;
		last = --sampling.end();
#define FILL_SUBDIV_LIMIT 300
		/* All (but one) input files in the 300-database and test_geometry* never go beyond N == 250.
		 * One file, input/test_geometry__svg_txt/doubleOutSpike.txt, goes to N == 268.
		 * */
#ifdef FILL_SUBDIV_LIMIT
		int N(0);
		while( (N < FILL_SUBDIV_LIMIT) && it1 != sampling.end() ) {
			++N;
#else
		while( it1 != sampling.end() ) {
#endif
			if( is_good_sample(it0->projected, it1->projected) ) {
				out(it1->projected);
				it0 = it1++;
			} else { // splitting
				auto it2 = it1; ++it2;
				if( (it1 != last) && ((it0->v | it2->v) > 0.0) && (clockFactor*det(it0->v,it2->v) > 0.0) ) {
					it1 = sampling.erase(it1); // it1 is now equal to it2
					Vec2d d = it1->v - it0->v;
					it1 = sampling.insert(it1, makeVP(it0->v + (0.66666 * d)));
					it1 = sampling.insert(it1, makeVP(it0->v + (0.33333 * d)));
				} else { // split in 2
					it1 = sampling.insert(it1, makeVP(splitArcInHalf(it0->v, it1->v, clockwise)));
				}
			}
		}
	}
}; // end of class Sampling

#undef FILL_SUBDIV_LIMIT

#endif // _CAVITIES_SAMPLING_H_
