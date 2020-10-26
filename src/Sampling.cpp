/* vim: set sts=0 ts=4 sw=4 noet : */

#include "Sampling.h"

#include <stack>
#include <algorithm>


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
moveTowardBisectorWithSegment(const Sample & sample, double delta, const Vec2d & P, const Vec2d & Q, bool & bad) {
	assert(sample.pos.eval(Utils::is_finite));
	assert(std::isfinite(sample.radius));
	assert(std::isfinite(delta));
	Vec2d U = sample.tangent;
	U.rotateCCW();
	Vec2d segment_tangent = (Q - P).normalized();
	Vec2d V = segment_tangent;
	V.rotateCW();
	Vec2d PO = sample.pos - P;
	const double denom = std::max(0.0, 1.0 - (U | V));
	double t;
	bool bad_seg = denom < 1e-6;
	if( ! bad_seg ) {
		t = ((PO | V) - delta) / denom;
		bad_seg = t < -1e-6;
	}
	if( bad_seg ) t = 0.0;
	bad = false;
	Vec2d bisector_pos = sample.pos + (t*U);
	double val = ( bisector_pos - P ) | segment_tangent;
	if( val < 0.0 ) { // maximal disk is not tangent in interior point of PQ
		return moveTowardBisectorWithDisk(sample, delta, Disk(P, 0), & bad);
	} else if( val > (Q - P).length() ) {
		Disk disk(Q, 0.0);
		return moveTowardBisectorWithDisk(sample, delta, Disk(Q, 0), & bad);
	} else if( bad_seg ) {
		bad = true;
		return {P, P, 0}; // dummy values
	}
	Vec2d bisector_tangent = (sample.tangent + segment_tangent).normalized();
	assert(bisector_pos.eval(Utils::is_finite));
	assert(bisector_tangent.eval(Utils::is_finite));
	Sample s = {bisector_pos, bisector_tangent, delta + distance(sample.pos, bisector_pos)};
	return s;
}

Sample
Sampling::
moveTowardBisectorWithDisk(const Sample & sample, double delta, const Disk & disk, bool * bad) {
	assert(sample.pos.eval(Utils::is_finite));
	assert(std::isfinite(sample.radius));
	assert(std::isfinite(delta));
	Vec2d OC = disk.center_ - sample.pos;
	double r = disk.radius_ + delta;
	Vec2d dir = sample.tangent;
	dir.rotateCCW();
	double denom = r + (dir | OC);
	double t = 0.5 * (OC.squaredLength() - r*r) / denom;
	if( bad != nullptr ) {
		*bad = (t <= 0.0) || (fabs(denom) < 1e-6);
		if( *bad ) return {dir, dir, 0}; // dummy value
	}
	Vec2d bisector_pos = sample.pos + (t*dir);
	Vec2d circle_tangent = (bisector_pos - disk.center_).normalized();
	circle_tangent.rotateCCW();
	if( (circle_tangent | sample.tangent) < 0.0 )
		circle_tangent.negate();
	Vec2d bisector_tangent = (sample.tangent + circle_tangent).normalized();
	assert(bisector_pos.eval(Utils::is_finite) && bisector_tangent.eval(Utils::is_finite));
	Sample s = {bisector_pos, bisector_tangent, delta + distance(sample.pos, bisector_pos)};
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
