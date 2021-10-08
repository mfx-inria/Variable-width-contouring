/* vim: set sts=0 ts=4 sw=4 noet : */

#ifndef _CAVITIES_UTILS_H_
#define _CAVITIES_UTILS_H_

#include "vec.h"

#include <vector>
#include <iostream> // std::cout

using Path = std::vector<Vec2d>;
using Paths = std::vector<Path>;

// For sorting radially using a reference center and direction
struct RadialCompare {
	Vec2d center, reference, dir;
	RadialCompare(const Vec2d & c, const Vec2d & r) : center(c), reference(r) {
		dir = (reference - center).normalized();
	}
	/* Returns true if and only if we meet 'a' strictly before 'b' when the ray from 'center' to
	 * 'reference' is rotating counter-clockwise from its start position (towards 'reference').
	 * In angle terms, if we see 'a' and 'b' as angles in the interval [0, 2*pi) (0 included, 2*pi
	 * excluded), the function returns true if and only if 'a' < 'b'.
	 */
	bool operator()(const Vec2d & a, const Vec2d & b);
};

// For conic piece of medial axis, in the familiar quadratic Bezier curve form.

// A plain quadratic Bezier curve.
struct QuadraticArc {
	Vec2d a, b, c; // 'a' and 'c' are endpoint, 'b' is tangents' crossing
	void reverse() {
		std::swap(a, c);
	}
	// Returns the length of the quadratic Bezier arc
	double analyticLength() const;
};

/* This struct stores a geometric description of the parabolic bisector between a point |point| and
   a segment.
   |x0| and |x1| are the signed distance from the |point| to the endpoints of the parabolic bisector
   arc, taken along |u| (x_i = dot(P_i - point, u)).
   So, we must have x0 <= x1.
   |u| has the same direction as the bisector arc.
   |v| points towards the half-plane that contains |point|, so that (u,v) can be CW or CCW, depending
   on the situation.
   |leftToRight| is true if when visualizing the point above the segment, then the bisector arc goes
   left-to-right. (The 3 points (arc-origin, arc-end, point) form a right turn.)
   */
struct PointSegmentBisector {
	Vec2d point;
	Vec2d u; // unit vector parallel to the segment
	Vec2d v; // unit vector orthogonal, from the segment to the point
	double h; // distance between segment and point
	double x0, x1;
	bool leftToRight;
	void print() const {
		std::cout << "PSB( point: " << point
			<< ", u: " << u << ", v: " << v
			<< ", h:" << h << ", x0: " << x0 << ", x1: " << x1 << ", leftToRight: " << leftToRight
			<< ')';
	}
};

template<typename T>
struct MinMax {
	typedef MinMax<T> Self;
	T min_;
	T max_;
	MinMax(const T & m) : min_(m), max_(m) {}
	MinMax(const T & min, const T & max) : min_(min), max_(max) {}
	void extendToCover(const Self & other) {
		min_ = min(min_, other.min_);
		max_ = max(max_, other.max_);
	}
};

typedef MinMax<double> Interval;

/* BBox means Bounding Box
 */
struct BBox : public MinMax<Vec2d>
{
	void extendToCover(const MinMax<Vec2d> & other) {
		min_ = min_.compMin(other.min_);
		max_ = max_.compMax(other.max_);
	}
	BBox(const Vec2d & m) : MinMax<Vec2d>(m) {};
	BBox(const Vec2d & min, const Vec2d & max) : MinMax<Vec2d>(min, max) {};
	Vec2d center() const;
	Vec2d size() const;
	double diagonal() const;
};

class Utils {
public:
	static const std::function<bool (double)> is_finite;
	// compute length of Bezier arc by recursive subdivision
	static double quadraticBezierNumLength(const QuadraticArc & q);
	static void test();

	// de Casteljeau split
	static void quadraticBezierSplitAtT(
			const QuadraticArc & q, double t,
			QuadraticArc & left, QuadraticArc & right);

	// bisection to find parameter corresponding to desired arc-length
	static void quadraticBezierSplitAtArcLength(
			const QuadraticArc & q, double len,
			QuadraticArc & left, QuadraticArc & right);

	// smallest of {r0, r1} that lies in the interval [min, max]
	static void smallestInInterval(double min, double max, int & numSols, double & r0, double & r1);

	// Solves polynomial a*X^2 + b*X + c = 0, puts solutions in r0 and r1
	static void solveQuadratic(double a, double b, double c, int & numSols, double & r0, double & r1);

	// Returns the orthogonal projection of point 'x' onto the line supporting 'a' and 'b'
	static Vec2d project(const Vec2d & x, const Vec2d & a, const Vec2d b);

	/* scale the input vector (passed by reference) so that its Linf norm is in [1/2, 1).
	   Returns the required scaling factor.
	   This is useful to compute the euclidean length of the vector without overflow.
	   Mostly useful with vectors that store integer coordinates, where the computation of the squared
	   length, prior to applying sqrt() is likely to overflow.
	 */
	static double reduceVec(Vec2d & v) {
		int exponent;
		double x = v.x();
		double y = v.y();
		if( std::abs(x) > std::abs(y) ) {
			x = std::frexp(x, &exponent);
			y = std::ldexp(y, -exponent);
		} else {
			y = std::frexp(y, &exponent);
			x = std::ldexp(x, -exponent);
		}
		v = Vec2d(x,y);
		return ldexp(1.0, exponent);
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - Segment Delaunay Graph

	// distance from point |query| to the line with direction |lineDir| passing through
	// |pointOnLine|.
	static double distanceToLine(const Vec2d & query, const Vec2d & pointOnLine, const Vec2d & lineDir);
};

#endif // _CAVITIES_UTILS_H_
