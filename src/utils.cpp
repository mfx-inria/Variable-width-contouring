/* vim: set sts=0 ts=4 sw=4 noet : */

#include "utils.h"

#include <cmath>
#include <stack>
#include <iostream>
#include <algorithm>

using namespace std;

#ifdef WIN32
double drand48()
{
  return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}
#endif

bool
RadialCompare::
operator()(const Vec2d & a, const Vec2d & b) {
	double ax = dir | (a-center);
	double bx = dir | (b-center);
	double ay = det(dir, a-center);
	double by = det(dir, b-center);
	if( ay > 0.0 ) {
		return (by < 0.0) || (bx < ax);
	} else if( ay < 0.0 ) {
		return (by < 0.0) && (ax < bx);
	} else { // ay == 0
		if( by > 0.0 ) {
			return ax > 0.0;
		} else if( by < 0.0 ) {
			return true;
		} else { // by == 0
			return (ax > 0.0) && (bx < 0.0);
		}
	}
}

double
QuadraticArc::
analyticLength() const {
	// The vectors B and A come from the polynomial representation of the quadratic arc :
	// P(t) = C + 2*B*t + A*t^2
	// with
	// 		C = a
	// 		B = b - a
	// 		A = a - 2*b + c
	// The derivative is P'(t) = 2*B + 2*A*t whose length is 2 * sqrt( Z + Y*t + X*t^2 )
	// with
	// 		X = A.A
	// 		Y = 2*B.A
	// 		Z = B.B

	// Notation Li = `length of L to the power i`

	const Vec2d A((a - b) + (c - b));
	const double A1 = A.length();
	if( A1 < 1e-6 ) { return (c-a).length(); }
	const double A2 = A.dot(A);

	const Vec2d  B  = b - a;
	const double B1 = B.length();
	const double B2 = B.dot(B);

	const Vec2d  F  = c - b;
	const double F1 = F.length();
	const double AB = A.dot(B);
	const double AF = A.dot(F);

	// TODO: Maybe replace /A2 by /A1/A1? (computing A2 might overflow)
	return
		(F1 * AF / A2) - (B1 * AB / A2)
		+
		(B2 - AB * AB / A2) / A1
		*
		(std::log(AF + A1*F1) - std::log(AB + A1*B1));
}


Vec2d
BBox::
center() const {
	return 0.5*Vec2d(max_.x()+min_.x(), max_.y()+min_.y());
}

Vec2d
BBox::
size() const {
	return Vec2d(max_.x()-min_.x(), max_.y()-min_.y());
}

double
BBox::
diagonal() const {
	return size().length();
}

const std::function<bool (double)> Utils::is_finite = [](double d) { return std::isfinite(d); };

double
Utils::
quadraticBezierNumLength(const QuadraticArc & iq) {
	QuadraticArc q = iq;
	double accum = 0.0;
	std::stack<QuadraticArc> jobs;
	int n(0);
	while( true ) {
		Vec2d ab = q.b - q.a;
		Vec2d bc = q.c - q.b;
		if( (ab.dot(bc) > 0.0) && (std::abs(det(ab, bc)) * 1e2 < ab.length()*bc.length()) ) {
			accum += (q.c-q.a).length();
			++n;
			if( jobs.empty() ) {
				break;
			}
			q = jobs.top();
			jobs.pop();
		} else {
			QuadraticArc p, r;
			quadraticBezierSplitAtT(q, 0.5, p, r);
			jobs.push(p);
			q = r;
		}
	}
	cerr << '(' << n << ')';
	return accum;
}

void
Utils::
// Compares the accuracies of two methods for computing the length of a quadratic Bezier arc.
test() {
	for( int i(0); i < 100; ++i ) {
		QuadraticArc q;
		q.a = 10.0 * Vec2d(drand48(), drand48());
		q.b = 10.0 * Vec2d(drand48(), drand48());
		q.c = 10.0 * Vec2d(drand48(), drand48());
		double anaL = q.analyticLength();
		double numL = quadraticBezierNumLength(q);
		cerr << numL << ' ' << anaL << endl;
	}
}

// de Casteljeau split
void Utils::quadraticBezierSplitAtT(
		const QuadraticArc & q, double t,
		QuadraticArc & left, QuadraticArc & right) {
	// left or right CAN be ref to q
	left.a = q.a;
	right.c = q.c;
	Vec2d b0   = lerp(q.a, q.b, t);
	Vec2d b1   = lerp(q.b, q.c, t);
	Vec2d mid  = lerp(b0, b1, t);
	left.b = b0;
	left.c = mid;
	right.a = mid;
	right.b = b1;
}

// bisection to find parameter corresponding to desired arc-length
void Utils::quadraticBezierSplitAtArcLength(
		const QuadraticArc & q, double len,
		QuadraticArc & left, QuadraticArc & right) {
	double lb(0.0), rb(1.0);
	double l = 0.0;
	while( ::fabs(l-len) > 0.0005 ) { // 0.0005 milli-meter
		double midparam = 0.5*(lb+rb);
		quadraticBezierSplitAtT(q, midparam, left, right);
		l = left.analyticLength();
		if( l > len ) {
			rb = midparam;
		} else {
			lb = midparam;
		}
	}
}

// smallest of {r0, r1} that lies in the interval [min, max]
void
Utils::
smallestInInterval(double min, double max, int & numSols, double & r0, double & r1) {
	// ASSUMES r0 <= r1 and numSols in {0, 1, 2}
	if( numSols == 0 ) return; // there is no valid solution
	// Here, numSols >= 1, so r0 is valid. Is it in the interval?
	if( (min <= r0) && (r0 <= max) ) { numSols = 1; return; }
	// r0 is not in the interval
	if( numSols == 1 ) { numSols = 0; return; } // r1 is not valid, invalidate all and return
	// Here, numSols == 2, so r1 is valid. Is it in the interval?
	if( (min <= r1) && (r1 <= max) ) { numSols = 1; r0 = r1; return; }
	// r1 was not in the interval, invalidate all and return
	numSols = 0;
}

void
Utils::
solveQuadratic(double a, double b, double c, int & numSols, double & r0, double & r1) {
	numSols = 0;
	const double eps = 1e-6;
	double disc;
	if( std::abs(a) < eps ) {
		if( std::abs(b) < eps ) return;
		r0 = - c / b;
		numSols = 1;
		return;
	}
	b /= 2.0;
	disc = b * b - a * c;
	if( disc < -eps*eps ) return;
	double sqrtOfDisc = std::sqrt(std::max(0.0, disc));
	if( sqrtOfDisc < eps ) {
		r0 = - b / a;
		numSols = 1;
		return;
	}
	double q;
	if( b < 0.0 ) q = - (b - sqrtOfDisc); // q > 0
	else	  q = - (b + sqrtOfDisc); // q < 0
	// q is the numerator with highest absolute value, which is better numerically for the division
	// by q below:
	r0 = q / a;
	r1 = c / q; // Because r0*r1 == c/a
	if( r0 > r1 ) std::swap(r0, r1);
	numSols = 2;
}

// returns projection of point x on line (a--b)
Vec2d
Utils::
project(const Vec2d & x, const Vec2d & a, const Vec2d b) {
	Vec2d ab = b - a;
	reduceVec(ab);
	return a + ab * ab.dot(x - a) / ab.squaredLength();
}

// squared distance from point |query| to the line with direction |lineDir| passing through
// |pointOnLine|.
double
Utils::
distanceToLine(const Vec2d & query, const Vec2d & pointOnLine, const Vec2d & lineDir) {
	Vec2d ortho = lineDir.rotatedCCW().normalized();
	return std::abs(ortho.dot(query - pointOnLine));
}

// - - - - - - - - - - - - - - BitangentComputer

bool
BitangentComputer::
rightRight(const BoundaryCircle & n1, const BoundaryCircle & n2, Vec2d & p1, Vec2d & p2) const {
	Vec2d D = n2.center() - n1.center();
	double L = D.length();
	double dr = n1.radius() - n2.radius();
	if( L*L + bitangent_eps <= dr*dr ) { return false; }
	D = D / L;
	dr /= L;
	const Vec2d Dcw = D.rotatedCW();
	const Vec2d v(dr*D + ::sqrt(abs(1.0-dr*dr))*Dcw);
	p1 = n1.center() + n1.radius() * v;
	p2 = n2.center() + n2.radius() * v;
	assert(std::isfinite(p1.x()) && std::isfinite(p1.y()) && std::isfinite(p2.x()) && std::isfinite(p2.y()));
	return true;
}

bool
BitangentComputer::
leftLeft(const BoundaryCircle & n1, const BoundaryCircle & n2, Vec2d & p1, Vec2d & p2) const {
	Vec2d D = n2.center() - n1.center();
	double L = D.length();
	double dr = n1.radius() - n2.radius();
	if( L*L + bitangent_eps <= dr*dr ) { return false; }
	D = D / L;
	dr /= L;
	const Vec2d Dccw = D.rotatedCCW();
	const Vec2d v(dr*D + ::sqrt(abs(1.0-dr*dr))*Dccw);
	p1 = n1.center() + n1.radius() * v;
	p2 = n2.center() + n2.radius() * v;
	assert(std::isfinite(p1.x()) && std::isfinite(p1.y()) && std::isfinite(p2.x()) && std::isfinite(p2.y()));
	return true;
}

bool
BitangentComputer::
rightLeft(const BoundaryCircle & n1, const BoundaryCircle & n2, Vec2d & p1, Vec2d & p2) const {
	Vec2d D = n2.center() - n1.center();
	double L = D.length();
	double sr = n1.radius() + n2.radius();
	if( L*L + bitangent_eps <= sr*sr ) { return false; }
	D = D / L;
	sr /= L;
	const Vec2d Dcw = D.rotatedCW();
	const Vec2d v(sr*D + ::sqrt(abs(1.0-sr*sr))*Dcw);
	p1 = n1.center() + n1.radius() * v;
	p2 = n2.center() - n2.radius() * v;
	assert(std::isfinite(p1.x()) && std::isfinite(p1.y()) && std::isfinite(p2.x()) && std::isfinite(p2.y()));
	return true;
}

bool
BitangentComputer::
leftRight(const BoundaryCircle & n1, const BoundaryCircle & n2, Vec2d & p1, Vec2d & p2) const {
	Vec2d D = n2.center() - n1.center();
	double L = D.length();
	double sr = n1.radius() + n2.radius();
	if( L*L + bitangent_eps <= sr*sr ) { return false; }
	D = D / L;
	sr /= L;
	const Vec2d Dccw = D.rotatedCCW();
	const Vec2d v(sr*D + ::sqrt(abs(1.0-sr*sr))*Dccw);
	p1 = n1.center() + n1.radius() * v;
	p2 = n2.center() - n2.radius() * v;
	assert(std::isfinite(p1.x()) && std::isfinite(p1.y()) && std::isfinite(p2.x()) && std::isfinite(p2.y()));
	return true;
}

bool
BitangentComputer::
touching(const BoundaryCircle & n1, const BoundaryCircle & n2, Vec2d & p) const {
	Vec2d D = n2.center() - n1.center();
	double L = D.length();
	double sr = n1.radius() + n2.radius();
	double squared_sr = sr * sr;
	if( ( L*L > squared_sr + bitangent_eps ) ||
		( L*L < squared_sr - bitangent_eps ) ) { return false; }
	if( L > 0.0 ) {
		D = D / L; // normalize D
	}
	p = n1.center() + n1.radius() * D;
	assert(std::isfinite(p.x()) && std::isfinite(p.y()));
	return true;
}

void
BitangentComputer::
operator()(const BoundaryCircle & left, const BoundaryCircle & right, Vec2d & p0, Vec2d & p1) const {
	Vec2d p;
	bool touch(false);
	if( touching(left, right, p) ) {
		p0 = p;
		p1 = p;
		touch = true;
	}
	if( left.passage_ == ToTheLeft ) {
		if( right.passage_ == ToTheLeft ) {
			if( ! leftLeft(left, right, p0, p1) ) { cerr << "Bad Left Left " << left << right << endl; }
		} else {
			if( touch ) return;
			if( ! leftRight(left, right, p0, p1) ) { cerr << "Bad Left Right " << left << right << endl; }
		}
	} else {
		if( right.passage_ == ToTheLeft ) {
			if( touch ) return;
			if( ! rightLeft(left, right, p0, p1) ) { cerr << "Bad Right Left " << left << right << endl; }
		} else {
			if( ! rightRight(left, right, p0, p1) ) { cerr << "Bad Right Right " << left << right << endl; }
		}
	}
}
