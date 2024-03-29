/* vim: set sts=0 ts=4 sw=4 noet tw=0 : */

#include "Sampling.h"

#include <stack>
#include <algorithm>

using namespace std;

double Sampling::len_threshold = 0.005;

// insert linear approximation of the quadratic curve, except for the first point.
void diceQuadratic(const Vec2d & focus, const QuadraticArc & q, std::vector<Sample> & path) {
	static std::stack<QuadraticArc> stack;
	bool work = true;
	QuadraticArc curQ = q;
	Sample s;
	s.tangent = Vec2d(0,0);
	while( work ) {
		double acLen = (curQ.c - curQ.a).length();
		Vec2d abn = (curQ.b - curQ.a).normalized();
		Vec2d bcn = (curQ.c - curQ.b).normalized();
		double cosine = abn | bcn;
		if( (acLen <= Sampling::len_threshold) || (cosine >= Sampling::ang_threshold) ) { // almost a line
			s.pos = curQ.c;
			s.radius = (curQ.c-focus).length();
			path.push_back(s);
			work = ! stack.empty();
			if( work ) {
				curQ = stack.top();
				stack.pop();
			}
		} else {
			QuadraticArc q0, q1;
			Utils::quadraticBezierSplitAtT(curQ, 0.5, q0, q1);
			stack.push(q1);
			curQ = q0;
		}
	}
}

inline
double distance(const Vec2d a, const Vec2d b) {
	return (a-b).length();
}

bool
Sampling::
is_good_sample(const Sample & s0, const Sample & s1) {
	double len = distance(s0.pos, s1.pos);
	bool too_small = len <= len_threshold;
	if( too_small ) return true;
	double cosine = s0.tangent | s1.tangent;
	//if( s0.tangent == s1.tangent ) return len <= 0.1;// <- TEMP HACK FOR BORDER SAMPLING
	double radius_ratio = ( s0.radius <= 0.0 || s1.radius <= 0.0 ) ? 0.0 :
		s0.radius > s1.radius ? (s0.radius/s1.radius) : (s1.radius/s0.radius);
	bool good_ratio = radius_ratio < 1.1; // by construction, 1.0 <= radius_ratio
	bool good_angle = cosine > ang_threshold;
	return good_angle && good_ratio;
}

Sample
Sampling::
moveTowardBisectorWithSegment(const Sample & sample, const Vec2d & P, const Vec2d & Q, bool & bad) {
	assert(sample.pos.eval(Utils::is_finite));
	assert(std::isfinite(sample.radius));

	const Vec2d U = sample.tangent.rotatedCCW();
	const Vec2d segment_tangent = (Q - P).normalized();
	const Vec2d V = segment_tangent.rotatedCW();
	const Vec2d PO = sample.pos - P;
	const double numerator = PO | V;
	if( numerator < 0.0 ) { bad = true; return {P, P, 0}; }
	const double denom = std::max(0.0, 1.0 - (U | V));
	double t;
	bool bad_seg = denom < 1e-6;
	if( ! bad_seg ) {
		t = numerator / denom;
		bad_seg = t < -1e-6;
	}
	if( bad_seg ) t = 0.0;
	bad = false;
	const Vec2d bisector_pos = sample.pos + (t*U);
	double val = ( bisector_pos - P ) | segment_tangent;
	if( val < 0.0 ) { // maximal disk is not tangent in interior point of PQ
		return moveTowardBisectorWithDisk(sample, 0.0, Disk(P, 0), & bad);
	} else if( val > (Q - P).length() ) {
		return moveTowardBisectorWithDisk(sample, 0.0, Disk(Q, 0), & bad);
	} else if( bad_seg ) {
		bad = true;
		return {P, P, 0}; // dummy values
	}
	const Vec2d bisector_tangent = (sample.tangent + segment_tangent).normalized();

	assert(bisector_pos.eval(Utils::is_finite));
	assert(bisector_tangent.eval(Utils::is_finite));

	Sample s = {bisector_pos, bisector_tangent, t};
	return s;
}

Sample
Sampling::
moveTowardBisectorWithDisk(const Sample & sample, double delta, const Disk & disk, bool * bad) {
	assert(sample.pos.eval(Utils::is_finite));
	assert(std::isfinite(sample.radius));
	assert(std::isfinite(delta));

	const Vec2d dir = sample.tangent.rotatedCCW();
	const Vec2d O = sample.pos - delta * dir; // cancel any previous inward parallel offset of the outer boundary
	const Vec2d OC = disk.center_ - O;
	const double R = disk.radius_;
	const double denom = R + (dir | OC);
	const double t = 0.5 * (OC.squaredLength() - R*R) / denom;
	if( bad != nullptr ) {
		*bad = (t + 1e-6 <= delta) || (fabs(denom) < 1e-6);
		if( *bad ) return {dir, dir, 0}; // dummy value
	}
	const Vec2d bisector_pos = sample.pos + ((t-delta)*dir);
	Vec2d circle_tangent = (bisector_pos - disk.center_).normalized().rotatedCCW();
	if( (circle_tangent | sample.tangent) < 0.0 )
		circle_tangent.negate();
	const Vec2d bisector_tangent = (sample.tangent + circle_tangent).normalized();

	assert(bisector_pos.eval(Utils::is_finite) && bisector_tangent.eval(Utils::is_finite));
	assert(t + 1e-6 > delta);

	Sample s = {bisector_pos, bisector_tangent, t};
	return s;
}

Vec2d
Sampling::
splitArcInHalf(const Vec2d & v0, const Vec2d & v1, bool clockwise) {
	Vec2d mid = 0.5*(v0 + v1);
	double lmid = mid.length();
	double d = det(v0, v1);
	if( lmid > 1e-3 ) {
		if( ((!clockwise)&&(d<0.0)) || (clockwise&&(d>0.0)) )
			mid.negate();
	} else {
		mid = clockwise ? (v1-v0) : (v0-v1);
		mid.rotateCCW();
	}
	return mid;
}

void
Sampling::
sampleSmoothPaths(
		const SmoothPaths & smoothPaths,
		vector<vector<Sample>> & samples,
		double minToolRadius) {
	samples.clear();
	for( const SmoothPath & sp : smoothPaths ) {
		size_t N = sp.size();
		if( N == 0 ) continue;
		samples.emplace_back();
		vector<Sample> & path = samples.back();
		auto inserter = [&](const Sample & s) { path.push_back(s); };
		if( N == 1 ) {
			const Vec2d c = sp[0].center();
			const double r =  sp[0].radius();
			sampleCircularArc(c+r*Vec2d(0,1), c-r*Vec2d(0,1), sp[0], inserter, /*clockwise=*/false);
			sampleCircularArc(c-r*Vec2d(0,1), c+r*Vec2d(0,1), sp[0], inserter, /*clockwise=*/false);
			continue;
		}
		Vec2d P, Q;
		BitangentComputer bmake(minToolRadius);
		bmake(sp[N-1], sp[0], P, Q);
		for( size_t i(0); i < N; ++i ) {
			size_t j = (i+1) % N;
			Vec2d p0, p1;
			bmake(sp[i], sp[j], p0, p1);
			// draw MATedge from Q to p0
			const Vec2d & c = sp[i].center();
			if( sp[i].radius() >= 1e-5 ) {
				if( sp[i].passage_ == ToTheRight ) {
					sampleCircularArc(Q-c, p0-c, sp[i], inserter, false);
				} else {
					sampleCircularArc(Q-c, p0-c, sp[i], inserter, true);
				}
			}
			Q = p1;
			Sample s{Q, (Q-p0).normalized(), 0.0};
			inserter(s);
		}
	}
}
