/* vim: set sts=0 ts=4 sw=4 noet tw=0 : */

#include "VariableWidthContouring.h"

#include <set>
#include <stack>
#include <unordered_set>
#include <algorithm>

using namespace std;

bool
VariableWidthContouring::
simplifyMAGeometry(Components * components, bool onlyCollapsed) {
	if( simplifyThreshold_ <= 1.0 ) return false;
	bool retVal = false;
	VertexSet to_remove;
	auto do_the_job = [&](MATvert * vert) {
		if( onlyCollapsed && ( ! vert->is_collapsed()) ) return;
		EdgeIterator freeway;
		if( 1 != degree(vert, freeway, /*walk_on_fire*/true) ) return;
		if( onlyCollapsed ) {
			const int collapsedDeg = statusDegree(vert, MatStatus::Collapsed, freeway);
			if( 1 != collapsedDeg ) return;
		}
		// Now we are sure that our vertex has degree-1, and is collapsed if required
		MATvert * neighbor = freeway->to();
		const int neighborCollapsedDeg = onlyCollapsed ? graph.statusDegree(neighbor, MatStatus::Collapsed) : degree(neighbor, true);
		if( 3 > neighborCollapsedDeg ) return; // we don't want to get a lonely collapsed vertex.
		if( vert->radius() + (vert->pos() - neighbor->pos()).length() < simplifyThreshold_ * neighbor->radius() ) {
			to_remove.insert(vert);
			retVal = true;
		}
	};
	if( nullptr == components ) {
		for ( MATvert & vert_ : graph.verts ) {
			do_the_job(& vert_);
		}
		graph.removeVertexSet(to_remove);
	} else {
		for( Component & component : *components ) {
			for( MATvert * vert : component.vertexSet ) {
				do_the_job(vert);
			}
			for( auto v : to_remove ) {
				component.vertexSet.erase(v);
			}
			graph.removeVertexSet(to_remove);
			to_remove.clear();
		}
	}
	return retVal;
}
bool
VariableWidthContouring::
simplifyCollapsedMAGeometry(Components * components) {
	return simplifyMAGeometry(components, true);
}
bool
VariableWidthContouring::
simplifyAllMAGeometry(Components * components) {
	return simplifyMAGeometry(components, false);
}

void
VariableWidthContouring::
simplifyMACombinatorics(Components * components) {
	bool writeEOL(false);
	for ( auto vertIt = graph.verts.begin(); vertIt != graph.verts.end(); ) {
		MATvert * vert = &(*vertIt);
		EdgeIterator freeway;
		if( 2 != degree(vert, freeway, true) ) {
			goto nextVert;
		} else { // the vertex has degree 2 (including all labels), we check if
			// we can delete if if it is redundant.
			EdgeIterator out0 = vert->edges_.begin();
			EdgeIterator out1 = out0; out1++;
			if( out0->to()->status != vert->status ) goto nextVert;
			if( out1->to()->status != vert->status ) goto nextVert;
			if( maoi_->hasBoundary(*out0) ) goto nextVert;
			if( maoi_->hasBoundary(*out1) ) goto nextVert;
			if( out0->site() != out1->twin()->site() ) goto nextVert;
			if( out1->site() != out0->twin()->site() ) goto nextVert;

			bool clipping(false);
			for( int j(0); j < static_cast<int>(vert->targets_.size()); ++j ) {
				Target & target = vert->targets_[j];
				if( target.isClipping ) {
					clipping = true;
					break;
				}
			}
			if( clipping ) goto nextVert;

			double r = vert->circumcircle.radius_;
			double r0 = out0->to()->circumcircle.radius_;
			double r1 = out1->to()->circumcircle.radius_;
			if( r0 > r1 ) swap(r0, r1);
			if( ( r < r0 ) || ( r1 < r ) ) goto nextVert; // local minima
			out0->twin()->twin_ = out1->twin();
			out0->twin()->to_   = out1->to();
			out1->twin()->twin_ = out0->twin();
			out1->twin()->to_   = out0->to();
			if( nullptr != components ) {
				for( Component & comp : *components ) {
					auto it = comp.vertexSet.find(vert);
					if( comp.vertexSet.end() != it ) {
						comp.vertexSet.erase(it);
					}
				}
			}
			vertIt = graph.verts.erase(vertIt);
			//writeEOL = true;
			continue;
		}
nextVert:
		++vertIt;
	}
	if( writeEOL ) {
		cerr << endl;
	}
}

function< Sample(const Sample &) >
VariableWidthContouring::
makeProjectorOnMA(const MATedge & arc) const {

	MATvert * vert = arc.from();

	MATvert * neighbor = arc.to();
	//MATvert * clearedTarget;

	if( arc.type == EdgeType::VertVert ) { // point/point bisector : closed form solution is too complex, should use dichotomic search.
		Vec2d bot = arc.site().point();
		Vec2d top = arc.twin()->site().point();
		Vec2d m(bot + 0.5 * (top-bot));
		Vec2d n(top-bot);
		n.normalize();
		Vec2d tgt = n;
		tgt.rotateCCW();
		return [=](const Sample & s) {
			Vec2d u = s.tangent;
			u.rotateCCW();
			double t = ((m - s.pos) | n) / (u|n);
			Sample ret;
			ret.pos = s.pos + (t * u);
			ret.tangent = ((tgt | s.tangent) > 0.0 ) ? tgt : tgt.negated();
			return ret;
		};
	} else if( arc.type == EdgeType::VertEdge ) { // segment/point bisector (parabola)
		PointSegmentBisector psb;
		graph.makePSBisector(arc, psb);
		Vec2d p(psb.point);
		Vec2d u(psb.u);
		Vec2d v(psb.v);
		double h(psb.h);
		double y = distanceToBoundary(arc, vert->circumcircle.center_);
		if( arc.twin()->site().is_point() )
			return [=](const Sample & s) {
				double ux = (s.pos - p)|u;
				double uy = (s.pos - p)|v;
				double A = ux;
				double B = -2.0*uy*h;
				double C = -ux*h*h;
				double r0, r1;
				int numSols;
				Utils::solveQuadratic(A, B, C, numSols, r0, r1);
				Sample ret;
				if( 0 == numSols ) {
					cerr << "OOOPS"; exit(-1);
				} else if( 1 == numSols ) {
					ret.pos = p - (0.5*h)*v;
					ret.tangent = u;
				} else {
					double x = r0*ux > 0.0 ? r0 : r1;
					double y = (h+x*x/h)/2.0-h;
					ret.pos = p + (x*u) + (y*v);
					ret.tangent = u + ((x/h)*v);
					ret.tangent.normalize();
				}
				return ret;
			};
		else
			return [=](const Sample & s) {
				double x = (s.pos-p)|u;
				double y = (h+x*x/h)/2.0-h;
				Sample ret;
				ret.pos = p + (x*u) + (y*v);
				ret.tangent = u + ((x/h)*v);
				ret.tangent.normalize();
				return ret;
			};
	} else { // segment/segment bisector
		Disk trim = vert->circumcircle;
		Vec2d m(trim.center_);
		Vec2d n(neighbor->circumcircle.center_);
		n = n - m;
		n.normalize();
		n.rotateCCW();
		return [=](const Sample & s) {
			Vec2d u = s.tangent;
			u.rotateCCW();
			double t = ((m - s.pos) | n) / (u|n);
			Sample ret;
			ret.pos = s.pos + (t * u);
			ret.tangent = n;
			ret.tangent.rotateCW();
			return ret;
		};
	}
}

const Target *
VariableWidthContouring::
trimArcToTargetForCollapse(const MATedge & arc, const Targets & targets, Disk & trim, bool & trimmed) const {
	MATvert * vert = arc.from();

	MATvert * neighbor = arc.to();
	const Target * clearedTarget = nullptr;

	trim = vert->circumcircle;
	trimmed = false;
	bool clearFound(false);
	double firstClear = 0.0;

	if ( targets.empty() ) {
		return nullptr;
	}

	auto updateClearedTarget = [&](const Target & target, double r0) {
		bool update = (! clearFound) || (r0 < firstClear);
		clearFound = true;
		if( update ) {
			firstClear = r0;
			clearedTarget = & target;
		}
	};

	if( arc.type == EdgeType::VertVert ) {
		Vec2d bot = arc.site().point();
		Vec2d top = arc.twin()->site().point(); // FOR COLLAPSE
		double h = (bot - top).length() / 2.0;
		Vec2d up = (top - bot)/(2.0*h);
		Vec2d right(up.y(), -up.x());
		double y = distanceToBoundary(arc, vert->circumcircle.center_);
		double xmin = (vert->circumcircle.center_-bot) * right;
		double xmax = (neighbor->circumcircle.center_-bot) * right;
		for( int ti(0); ti < static_cast<int>(targets.size()); ++ti ) {
			const Target & target = targets[ti];
			const double targetDist = 2.0 * target.radius - target.medialAxisDisk.radius_;
			if( targetDist <= 0.0 ) continue;
			double X0 = (target.medialAxisDisk.center_ - bot) * right;
			double Y0 = (target.medialAxisDisk.center_ - bot) * up - h;
			const double UU = targetDist * targetDist - Y0 * Y0;
			if( UU <= 0.0 ) continue;
			double U = ::sqrt(UU);
			double r0 = X0 - U, r1 = X0 + U;
			int numSols(2);
			Utils::smallestInInterval(xmin, xmax, numSols, r0, r1);
			if( numSols > 0 ) { // FIXME: we should check that the
				//	solutions are indeed solutions. because we
				//	squared 2 numbers in order to obtain
				//	the quadratic polynomial...
				updateClearedTarget(target, r0);
			}
		}
		if( clearFound ) {
			//cerr << " clear Found at " << firstClear;
			trim.center_ += (firstClear-xmin) * right;
			double y2 = distanceToBoundary(arc, trim.center_);
			trim.radius_ += y2 - y;
			trimmed = true;
		}
	} else if( arc.type == EdgeType::EdgeEdge ) { // segment/segment bisector
		Vec2d right = arc.site().segment_vector();
		Vec2d  left =  arc.twin()->site().segment_vector(); // FOR COLLAPSE
		Vec2d dest   = neighbor->circumcircle.center_;
		Vec2d dir = dest - trim.center_;
		//if( right * dir < 0.0 ) right = - right;
		//if( left * dir < 0.0 ) left = - left;
		double dirLen = dir.length();
		dir = dir / dirLen;
		for( int ti(0); ti < static_cast<int>(targets.size()); ++ti ) {
			const Target & target = targets[ti];
			const double targetDist = 2.0 * target.radius - target.medialAxisDisk.radius_;
			if( targetDist <= 0.0 ) continue;
			double X0 = (target.medialAxisDisk.center_ - trim.center_) * dir;
			double Y0 = det(target.medialAxisDisk.center_ - trim.center_, dir);
			const double UU = targetDist * targetDist - Y0 * Y0; // FIXME: use a hypot-like function
			if( UU <= 0.0 ) continue;
			double U = ::sqrt(UU);
			double r0 = X0 - U, r1 = X0 + U;
			int numSols(2);
			Utils::smallestInInterval(0.0, dirLen, numSols, r0, r1);
			if( numSols > 0 ) {
				updateClearedTarget(target, r0);
			}
		}
		if( clearFound ) {
			trim.center_ += firstClear * dir;
			trim.radius_ = (firstClear * neighbor->circumcircle.radius_ + (dirLen-firstClear)*trim.radius_) / dirLen;
			trimmed = true;
		}
	} else { // segment/point bisector (parabola)
		assert(arc.type == EdgeType::VertEdge);
		PointSegmentBisector psb;
		graph.makePSBisector(arc, psb); // FOR COLLAPSE
		double h = psb.h;
		double y = distanceToBoundary(arc, vert->circumcircle.center_);
		for( int ti(0); ti < static_cast<int>(targets.size()); ++ti ) {
			const Target & target = targets[ti];
			double X0 = (target.medialAxisDisk.center_ - psb.point) * psb.u;
			double Y0 = (target.medialAxisDisk.center_ - psb.point) * psb.v + h;
			// Equation to solve is polynomial of degree 4.
			// I think that in this problem the value is actually monotonous (property of Voronoi diagrams/MA)
			// so we go for a direct dichotomic search.
			auto val = [&](double qx) {
				double qy = h/2.0 + qx*qx/2.0/h;
				Vec2d diff = Vec2d(qx, qy) - Vec2d(X0, Y0);
				return diff.length() - 2.0*target.radius + target.medialAxisDisk.radius_;
			};
			if( val(psb.x0) >= 0.0 ) {
				cerr << "BIZARRE value " << val(psb.x0) << " at " << arc.from()->circumcircle.center_
					<< " MA disk=" << target.medialAxisDisk.center_ << ", rad="
					<< target.medialAxisDisk.radius_
					<< endl;
				continue;
			}
			if( val(psb.x1) < 0.0 ) continue;
			double lo = psb.x0, hi = psb.x1;
			double mid = 0.5*(hi+lo);
			while( hi-lo > 1e-5 ) {
				if( val(mid) < 0.0 )
					lo = mid;
				else
					hi = mid;
				mid = 0.5*(hi+lo);
			}
			updateClearedTarget(target, mid);
		}
		if( clearFound ) {
			double newY = h/2+firstClear*firstClear/2.0/h;
			trim.center_ = psb.point + firstClear * psb.u + (newY-h) * psb.v;
			trim.radius_ += newY - y;
			trimmed = true;
		}
	}
	return clearedTarget;
}

static bool debug = false;

MATvert *
VariableWidthContouring::
trimArcToTarget(const MATedge & arc, const Targets & targets, Disk & trim, bool & trimmed) const {
	double d = 2.0 * minToolRadius_;

	MATvert * vert = arc.from();

	MATvert * neighbor = arc.to();
	MATvert * clearedTarget = nullptr;

	trim = vert->circumcircle;
	trimmed = false;
	bool clearFound(false);
	double firstClear = 0.0;

	if ( targets.empty() ) {
		return clearedTarget;
	}

	//bool debug = vert->pos().x() > 55 && vert->pos().x() < 56 && vert->pos().y() > 998 && vert->pos().y() < 1000;

	auto updateClearedTarget = [&](const Target & target, double r0) {
		bool update = (! clearFound) || (r0 < firstClear);
		clearFound = true;
		if( update ) {
			firstClear = r0;
			if( target.medialAxisDisk.radius_ <= maxToolRadius_ ) // we have cleared the actual osculating circle
				clearedTarget = target.vertexSource;
			else
				clearedTarget = nullptr;
		}
	};

	if( arc.type == EdgeType::VertVert ) {
		Vec2d bot = arc.site().point();
		Vec2d top = arc.twin()->site().point();
		double h = (bot - top).length() / 2.0;
		Vec2d up = (top - bot)/(2.0*h);
		Vec2d right(up.y(), -up.x());
		double y = distanceToBoundary(arc, vert->circumcircle.center_);
		double xmin = (vert->circumcircle.center_-bot) * right;
		double xmax = (neighbor->circumcircle.center_-bot) * right;
		for( int ti(0); ti < static_cast<int>(targets.size()); ++ti ) {
			const Target & target = targets[ti];
			double QX = (target.medialAxisDisk.center_ - bot) * right;
			double QY = (target.medialAxisDisk.center_ - bot) * up - h;
			double Delta = y - vert->radius();
			double U = /* W(q) = 2 tgt.radius - d */ 2.0 * target.radius - d  - /* \dsi(q) */ target.medialAxisDisk.radius_ - Delta;
			double K = QX*QX + QY*QY - h*h - U*U;
			double A = QX*QX - U*U;
			double B = -QX*K;
			double C = K*K/4.0 - U*U*h*h;
			double r0(-1), r1(-2);
			int numFakeSols, numSols;
			Utils::solveQuadratic(A, B, C, numFakeSols, r0, r1);
			if( debug ) {
				cerr << "for tgt " << target.medialAxisDisk.center_
					<< ": numFakeSols=" << numFakeSols << ", r0=" << r0 << ", r1=" << r1
					<< "\n A=" << A << ", B=" << B << ", C=" << C << ", discr=" << ((B/2)*(B/2)-A*C)
					<< endl;
			}
			//cerr << ' ' << numFakeSols << " sols, " << r0 << " and " << r1 << ". ";
			numSols = numFakeSols;
			if( numFakeSols > 0 ) {
				double X = r0;
				double right = std::hypot(X, h)+U;
				double check = (X-QX)*(X-QX) + QY*QY - right*right;
				if( abs(check) > 1e-4 ) { // r0 is NOT a solution
					numSols--; r0 = r1;
				}
				if( numFakeSols > 1 ) {
					double X = r1;
					double right = std::hypot(X, h)+U;
					double check = (X-QX)*(X-QX) + QY*QY - right*right;
					if( abs(check) > 1e-4 ) { // r1 is NOT a solution
						numSols--;
					}
				}
			}
			Utils::smallestInInterval(xmin, xmax, numSols, r0, r1);
			if( numSols > 0 ) {
				updateClearedTarget(target, r0);
			}
		}
		if( clearFound ) {
			//cerr << " clear Found at " << firstClear;
			trim.center_ += (firstClear-xmin) * right;
			double y2 = distanceToBoundary(arc, trim.center_);
			trim.radius_ += y2 - y;
			trimmed = true;
		}
	} else if( arc.type == EdgeType::VertEdge ) { // segment/point bisector (parabola)
		PointSegmentBisector psb;
		graph.makePSBisector(arc, psb);
		double h = psb.h;
		double y = distanceToBoundary(arc, vert->circumcircle.center_);
		for( int ti(0); ti < static_cast<int>(targets.size()); ++ti ) {
			const Target & target = targets[ti];
			double X0 = (target.medialAxisDisk.center_ - psb.point) * psb.u;
			double Y0 = (target.medialAxisDisk.center_ - psb.point) * psb.v + h;
			double U = 2.0 * target.radius - target.medialAxisDisk.radius_ - y - d + vert->circumcircle.radius_;
			double A = 1.0 - ( U + Y0 ) / h;
			double B = -2.0 * X0;
			double C = -h * U - U * U + X0 * X0 - h * Y0 + Y0 * Y0;
			double r0, r1;
			int numSols;
			Utils::solveQuadratic(A, B, C, numSols, r0, r1);
			Utils::smallestInInterval(psb.x0, psb.x1, numSols, r0, r1);
			if( numSols > 0 ) {
				updateClearedTarget(target, r0);
			}
		}
		if( clearFound ) {
			double newY = h/2+firstClear*firstClear/2.0/h;
			trim.center_ = psb.point + firstClear * psb.u + (newY-h) * psb.v;
			trim.radius_ += newY - y;
			trimmed = true;
		}
	} else { // segment/segment bisector
		assert(arc.type == EdgeType::EdgeEdge);
		Vec2d right = arc.site().segment_vector();
		Vec2d  left =  arc.twin()->site().segment_vector();
		Vec2d dest   = neighbor->pos();
		Vec2d dir = dest - trim.center_;
		if( right * dir < 0.0 ) right = - right;
		if( left * dir < 0.0 ) left = - left;
		double cos2a = (right * left) / right.length() / left.length();
		double squaredSina = max(0.0, (1.0-cos2a)*0.5);
		double sina = sqrt(squaredSina);
		if( det(right, left) < 0.0 ) {
			sina = - sina;
		}
		double dirLen = dir.length();
		dir = dir / dirLen;
		for( int ti(0); ti < static_cast<int>(targets.size()); ++ti ) {
			const Target & target = targets[ti];
			double X0 = (target.medialAxisDisk.center_ - trim.center_) * dir;
			double Y0 = det(target.medialAxisDisk.center_ - trim.center_, dir);
			double U = trim.radius_ + 2.0 * target.radius - target.medialAxisDisk.radius_ - d;
			double A = 1.0 - squaredSina; // cos^2
			double B = - 2.0 * (X0 + U * sina);
			double C = X0*X0 + Y0*Y0 - U*U;
			if( debug ) {
				cerr << "A=" << A << ", B=" << B << ", C=" << C << ", U=" << U
					<< ", dirLen=" << dirLen << ", dir=" << dir
					<< endl;
			}
			double r0, r1;
			int numSols;
			Utils::solveQuadratic(A, B, C, numSols, r0, r1);
			if( debug ) {
				cerr << "numSols=" << numSols << ", r0=" << r0 << ", r1=" << r1 << endl;
			}
			Utils::smallestInInterval(0.0, dirLen, numSols, r0, r1);
			if( numSols > 0 ) {
				if( debug ) {
					cerr << "numSols=" << numSols << ", r0=" << r0 << ", r1=" << r1 << endl;
				}
				updateClearedTarget(target, r0);
			}
		}
		if( clearFound ) {
			trim.center_ += firstClear * dir;
			trim.radius_ += firstClear * sina;
			trimmed = true;
		}
	}
	return clearedTarget;
}

void
VariableWidthContouring::
shrinkMA(Components & components, vector<Disk> & cheekPrecursors, bool sharpCut) {
	initializeTargets(/*walk_on_fire*/false);
	for( Component & component: components ) {
		setCurrentComponent(component);
		// trim branches, starting at each leaf.
		component.uniformOffset = trimLeafs(/*trim vertices*/true);
		collapse(cheekPrecursors, sharpCut); // collapse untrimmed parts which wouldn't fit another minToolDiameter
		// trim the collapse targets
		double clearance_update = trimLeafs(/*collapse vertices*/false);
		component.uniformOffset = min(component.uniformOffset, clearance_update);
		/* Finally, clear the |targets_| in the remaining Voronoi vertices and reduce the radius of non-leaf Voronoi
		 * vertices. */
		reduceRadii();
	}
}

void
VariableWidthContouring::
shrinkMA_trim(Components & components) {
	initializeTargets(/*walk_on_fire*/false);
	// Separate loop for each sub-step so that we can DRAW the result at each
	// sub-step for all components
	for( Component & component: components ) {
		setCurrentComponent(component);
		// trim branches, starting at each leaf.
		component.uniformOffset = trimLeafs(/*trim vertices*/true);
	}
}
bool
VariableWidthContouring::
shrinkMA_collapse(Components & components, vector<Disk> & cheekPrecursors, bool sharpCut) {
	bool hadSomeCollapse = false;
	for( Component & component: components ) {
		setCurrentComponent(component);
		bool r = collapse(cheekPrecursors, sharpCut); // collapse untrimmed parts which wouldn't fit another minToolDiameter
		if( r && (component.vertexSet.size() > 1) ) {// && (component.uniformOffset > 2.0 * minToolRadius_ + 1e-4) ) {
			//cerr << "In the presence of collapse, the component's uniformOffset must be minimal ("
			//	<< (2.0*minToolRadius_) << ") but it is "
			//	<< component.uniformOffset << ". Changing it..." << endl;
			component.uniformOffset = 2 * minToolRadius_;
		}
		hadSomeCollapse = hadSomeCollapse || r;
	}
	return hadSomeCollapse;
}
void
VariableWidthContouring::
shrinkMA_cheeks(Components & components) {
	for( Component & component: components ) {
		setCurrentComponent(component);
		// trim the collapse targets
		double clearance_update = trimLeafs(/*collapse vertices*/false);
		component.uniformOffset = min(component.uniformOffset, clearance_update);
	}
}
void
VariableWidthContouring::
shrinkMA_globalOffset(Components & components) {
	for( Component & component: components ) {
		setCurrentComponent(component);
		/* Finally, clear the |targets_| in the remaining Voronoi vertices and reduce the radius of non-leaf Voronoi
		 * vertices. */
		reduceRadii();
	}
}

void
VariableWidthContouring::
initializeTargets(bool walk_on_fire) {
	for( MATvert & vert_ : graph.verts) {
		MATvert * vert = & vert_;
		//Vec2d pos = vert->circumcircle.center_;
		//bool debug = (pos.x() > -8.0) && (pos.x() < -7.0) && (pos.y() > 12.0) && (pos.y() < 13.0);
		//if( debug ) cerr << "Init targets " << (walk_on_fire?"on fire":"on normal") << " at " << pos;
		vert->targets_.clear();
		const int deg = degree(vert, walk_on_fire);
		bool addTarget = (1 == deg) && (!walk_on_fire || vert->is_collapsed());
		//if( debug ) cerr << ". degree=" << deg << ". 1st test: " << addTarget;
		if( (! addTarget) && (! walk_on_fire) && (2 <= deg) ) {
			for( EdgeIterator edge = vert->edges_.begin(); edge != vert->edges_.end(); ++edge ) {
				if( edge->to()->is_destroyed() ) continue;
				if( ! maoi_->hasToTheRightBoundary(*edge) ) continue;
				// We have a protruding ToTheRight boundary circle at a normal, degree>=2 vertex. Must add it as a target
				// FIXME: should we also do it when walk[ing]_on_fire?
				addTarget = true;
			}
			//if( debug ) cerr << ". 2nd test: " << addTarget;
		}
		if( addTarget ) {
			if( walk_on_fire ) { // reset the fullyCleared_ flag
				vert->fullyCleared_ = false;
			}
			if( walk_on_fire && debug ) cerr << "Adding Collapsed target at " << vert->pos() << endl;
			double targetRadius = min(vert->circumcircle.radius_, maxToolRadius_);
			vert->targets_.emplace_back(vert->circumcircle, targetRadius);
			vert->targets_.back().vertexSource = vert;
		}
		//if( debug ) cerr << endl;
	}
}

double
VariableWidthContouring::
trimLeafs(bool trimVertices) {

	VertexSet & component_verts = currentComponent_->vertexSet;

	list<MATvert *> leafs;

	for ( MATvert * vert : component_verts ) {
		if( vert->is_destroyed() ) continue; // vert was collapsed
		if( vert->targets_.empty() ) continue;
		const int deg = degree(vert);
		if( 1 == deg ) { // init leaf
			leafs.emplace_back(vert);
		}
	}

	double global_clearance = 2 * maxToolRadius_;

	while( ! leafs.empty() ) {
		MATvert * vert = leafs.back();
		leafs.pop_back();

		//bool debug = vert->pos().x() > 55 && vert->pos().x() < 56 && vert->pos().y() > 998 && vert->pos().y() < 1000;

		if ( vert->targets_.empty() ) continue;
		EdgeIterator freeway;
		int deg = degree(vert, freeway); // degree == 1 in general so freeway != vert->edges_.end()
		if( 1 < deg ) {
			cerr << "ERROR: having a leaf with degree " << deg << " > 1. BUG." << endl;
			exit(-1);
		}
		if( 0 == deg ) { // It means the neighbor leaf was shrunk and died.
			continue;
		}

		//debug = (vert->pos().x() > 114) && (vert->pos().x() < 116);
		//debug = debug && (vert->pos().y() > 951) && (vert->pos().y() < 953);

		bool trimmed(false);
		Disk trim;
		MATvert * firstCleared;
		if( debug ) cerr << "Trimming... have " << vert->targets_.size() << " targets." << endl;
		firstCleared = trimArcToTarget(*freeway, vert->targets_, trim, trimmed);
		if( trimmed ) { // we succesfully found a position along the branch |freeway| where we can trim the branch
			MATvert & inserted_vert = graph.insert(currentComponent_, freeway, trim);
			if( firstCleared != nullptr ) {
				firstCleared->fullyCleared_ = true;
				if( debug ) cerr << "[TrimFullyCleared at "<<firstCleared->pos()<<"]";
				global_clearance = 2 * minToolRadius_;
			} else {
				vert->targets_.swap(inserted_vert.targets_); // works because inserted_vert.targets_ is empty.
			}
			//vert->targets_.clear();
			if( trimVertices ) {
				vert->trim();
				if( debug ) cerr << "<cut-trimmed at " << inserted_vert.pos() << ",r=" << inserted_vert.radius() << '>' << endl;
			} else {
				vert->collapse();
			}
		} else {
			/* trimming failed: push the next MA vertex onto |leafs| in order to continue the
			 * trimming further up the branch. If the degree of |neighbor| is > 1, we do nothing: the targets_ of
			 * |neighbor| will be picked up by another branch arriving there, in the spirit of the "grassfire
			 * transform" */
			MATvert * neighbor = freeway->to();

			// handle targets
			Targets & neighbor_targets = neighbor->targets_;
			neighbor_targets.insert(neighbor_targets.end(), vert->targets_.begin(), vert->targets_.end());
			vert->targets_.clear();
			if( trimVertices ) {
				vert->trim();
				if( debug ) cerr << "<prop-trimmed " << vert->pos() << '>' << endl;
			} else {
				vert->collapse();
			}

			if( ( ! neighbor->is_destroyed()) && (degree(neighbor) == 1) ) {
				leafs.push_back(neighbor);
			}
		}
	}
	// Compute new global clearance and clear the targets of normal-degree 1 vertices.
	// Targets at normal-degree 2 vertices stay there.
	for ( MATvert * vert : component_verts ) {
		for ( const Target & target : vert->targets_ ) {
			double localTrackWidth = (vert->circumcircle.center_ - target.medialAxisDisk.center_).length()
				+ target.medialAxisDisk.radius_ - vert->circumcircle.radius_;
			double targetWidth = 2.0 * target.radius;
			double target_clearance = targetWidth - localTrackWidth;
			global_clearance = min(global_clearance, target_clearance);
		}
		//if( 1 == degree(vert) ) vert->targets_.clear();
	}
	return global_clearance;
}

void
VariableWidthContouring::
reduceRadii() {
	VertexSet & vertices = currentComponent_->vertexSet;
	double global_clearance = currentComponent_->uniformOffset;
	// global clearance should never be less than 2*minToolRadius_, but when we hit a fully cleared
	// target, we set it to zero, so...
	double reduction = max(global_clearance, 2 * minToolRadius_); // most optimal step distance

	//cerr << "REDUCTION = " << reduction;

	for ( MATvert * vert : vertices ) {
		if( vert->is_destroyed() ) continue;
		if( vert->circumcircle.radius_ - reduction < 2 * minToolRadius_ // reducing maximally will cause impossibly small feature
		) {
			reduction = min(reduction, vert->circumcircle.radius_ - 2 * minToolRadius_);
			reduction = max(reduction, 2 * minToolRadius_); // don't make reduction zero if the step size is the feature size
		}
	}

	//cerr << " then " << reduction << endl;

	currentComponent_->uniformOffset = reduction;
	for ( MATvert * vert : vertices ) {
		vert->targets_.clear();
		if( vert->is_destroyed() ) continue;
		vert->circumcircle.radius_ -= reduction;
		if (vert->circumcircle.radius_ < 1.99 * minToolRadius_) {
			cerr << "SHOULD NEVER HAPPEN" << endl;
			vert->circumcircle.radius_ = 0; // avoid generating impossibly small features
		}
	}
}

extern int gCurrentStep;

bool
VariableWidthContouring::
collapse(vector<Disk> & cheekPrecursors, bool sharpCut) {
	const double goodGradient = minToolRadius_ / maxToolRadius_;
	bool hadSomeCollapse = false;

	unordered_set< MATvert * > & component_verts = currentComponent_->vertexSet;

	vector<MATvert *> to_add_in_component;
	stack<EdgeIterator> to_check;

	auto checkTargets = [&](MATvert * vert) {
		if( ! vert->is_trimmed() ) return;
		for( auto target : vert->targets_ ) {
			target.vertexSource->fullyCleared_ = false;
		}
	};

	auto markCollapsedAndPropagate = [&](EdgeIterator edge) { // Propagate collapse
		MATvert * to = edge->to();
		assert(to->circumcircle.radius_ <= maxCollapseRadius()*(1.0+1.0e-6));
		to->collapse();
		for( EdgeIterator e = to->edges_.begin(); e != to->edges_.end(); ++e ) {
			if( e->to()->is_collapsed() ) continue; // already collapsed
			if( e == edge->twin() ) continue; // don't go back!
			to_check.push(e);
		}
	};

	//bool debug = gCurrentStep == 4;

	for( MATvert * vert : component_verts ) {
		//if( debug ) cerr << "[Collapse] Looking at " << vert->pos() << ": " << endl;

		if (vert->is_destroyed()) continue; // vert is trimmed
		if (vert->radius() >= minCollapseRadius()) continue; // don't collapse

		vert->collapse();
		hadSomeCollapse = true;

		//double x = vert->pos().x();
		//double y = vert->pos().y();
		const bool debug = false;//(x>53)&&(x<60);//&&(y>59)&&(y<60);

		if( debug ) cerr << "[Collapse] Initiated at " << vert->pos() << endl;

		assert(to_check.empty());
		for( EdgeIterator e = vert->edges_.begin(); e != vert->edges_.end(); ++e ) {
			if ( e->to()->is_collapsed() ) continue; // already collapsed
			to_check.push(e);
		}

		while( ! to_check.empty() ) {
			EdgeIterator edge = to_check.top(); to_check.pop();
			if( edge->to()->is_collapsed() ) continue; // avoid going around in circles

			// FIXME: maybe a fixme... I comment the line because AFAIK, no target should be cleared at this point.
			// Useless ones have been removed already at the end of trimming and the one remaining are useful.
			// Also, if we clear anything, that should be edge->to()->targets_, not from().
			//// In preparation for putting cheeks as new targets
			////edge->from()->targets_.clear();
			// Instead, let's check that we never did harm:
			// End of fixme.

			if (debug ) {
				cerr << "At " << edge->from()->pos() << " (r=" << edge->from()->radius() << "). ";
				cerr << "Checking neighbor at " << edge->to()->pos() << " (r=" << edge->to()->radius() << "): ";
			}

			double r0 = edge->from()->radius();
			double r1 = edge->to()->radius();
			bool fromIsVerySmall = r0 <= minCollapseRadius() - 1e-4;
			bool toIsSmall = r1 < minCollapseRadius() + 1e-4;

			//if( edge->from()->radius() <= minCollapseRadius() ) { // too small, propagate
			if( fromIsVerySmall || toIsSmall ) { // propagate
				if( debug ) {
					cerr << "smalltown, propagate\n";
				}
				checkTargets(edge->to());
				Disk disk;
				if( findPointAtDistance(*edge, minCollapseRadius(), disk) > 0.0 ) {
					if( debug ) {
						cerr << "Found cut point at " << disk.center_ << "(r=" << disk.radius_ << ')' << endl;
					}
					graph.insert(nullptr, edge, disk);
					to_add_in_component.push_back(edge->to());
					if( sharpCut ) continue;
				} else {
					if( debug ) {
						cerr << " No cut at minCollapseRadius=" << minCollapseRadius() << endl;
					}
				}
				markCollapsedAndPropagate(edge);
				continue;
			} else {
				if( debug ) {
					cerr << "big enough\n";
				}
			}
			if( sharpCut ) continue;
			// Here:
			// 		edge->from()->radius() > minCollapseRadius() - 1e-4
			// and	edge->to()->radius() >= minCollapseRadius() + 1e-4
			double goodGradientPos = -1; Disk goodGradientDisk;
			double maxRadiusPos; Disk maxRadiusDisk;
			double gradient = edge->getRadiusGradientAtStart();
			if( debug ) {
				cerr << "Gradient at " << edge->from()->pos() << " is " << gradient << ", goodGrad is " << goodGradient << endl;
			}
			if( gradient >= goodGradient ) {
				double rad = minToolRadius_*(2.0+1.0/gradient);
				if( EdgeType::EdgeEdge == edge->type ) rad *= 1.05;
				if( debug ) {
					cerr << "We have a big-enough gradient. Radius at " << edge->from()->pos() << " is " << edge->from()->radius() << " and goodRadius is " << rad << endl;
				}
				if( edge->from()->radius() >= rad ) {
					if( debug ) {
						cerr << "\"From\" is large enough, reverted to Normal label\n";
					}
					// edge->from() is all good. Revert it back to normal
					edge->from()->labelAsNormal();
					continue;
				} else {
					if( debug ) {
						cerr << "Finding good gradient...\n";
					}
					goodGradientPos = findPointAtDistance(*edge, rad, goodGradientDisk);
				}
			}
			checkTargets(edge->to());
			fromIsVerySmall = r0 <= maxCollapseRadius() - 1e-4;
			toIsSmall = r1 < maxCollapseRadius() + 1e-4;
			if( fromIsVerySmall || toIsSmall ) {
				if( debug ) { cerr << "~A~"; }
				maxRadiusPos = findPointAtDistance(*edge, maxCollapseRadius(), maxRadiusDisk);
			} else {
				if( debug ) { cerr << "~B~"; }
				maxRadiusPos = 0;
			}
			if( (goodGradientPos >= 0) && ( (maxRadiusPos < 0) || (goodGradientPos < maxRadiusPos) ) ) {
				if( debug ) {
					cerr << "Good gradient further from " << edge->from()->pos() << endl;
				}
				graph.insert(nullptr, edge, goodGradientDisk);
			} else if( maxRadiusPos >= 0.0 ) {
				if( debug ) { cerr << "~C~"; }
				if( maxRadiusPos > 0.0 ) {
					if( debug ) { cerr << "~MAX RAD REACHED~"; }
					graph.insert(nullptr, edge, maxRadiusDisk);
				} else {
					// from() has radius very close to maxCollapseRadius() and
					// to() has radius >= maxCollapseRadius() + 1e-4,
					// So we revert from() to Normal label and stop here.
					edge->from()->labelAsNormal();
					continue;
				}
			} else {
				if( debug ) { cerr << "markCollapsedAndProp.....\n"; }
				markCollapsedAndPropagate(edge);
				continue;
			}
			assert(edge->to()->radius() <= maxCollapseRadius());
			to_add_in_component.push_back(edge->to());

#if 0 // Add virtual targets_ (cheek precursors) for making as-large-as-possible cheeks.
	  // Commented-out as currently, it does not improve thing in any of the tests.
	  // On the contrary, it make some voids bigger (especially in test-6.txt)
	  // without removing other voids.
			if( edge->to()->radius() <= maxCollapseRadius() - 1e-4 ) continue;
			// Compute and insert cheeks as targets for subsequent collapse-trimming
			Vec2d p = edge->to()->circumcircle.center_;
			double targetRadius = edge->to()->circumcircle.radius_ / 2.0;
			Vec2d a;
			if( edge->site().is_point() ) a = edge->site().point();
			else a = Utils::project(p, edge->site().source(), edge->site().destination());
			Vec2d dir = (a - p).normalized();
			cheekPrecursors.emplace_back(p + targetRadius * dir, targetRadius);
			edge->to()->targets_.emplace_back(cheekPrecursors.back(), targetRadius);
			Vec2d b;
			MATedge * twin = edge->twin();
			if( twin->site().is_point() ) b = twin->site().point();
			else b = Utils::project(p, twin->site().source(), twin->site().destination());
			dir = (b - p).normalized();
			cheekPrecursors.emplace_back(p + targetRadius * dir, targetRadius);
			edge->to()->targets_.emplace_back(cheekPrecursors.back(), targetRadius);
#endif
		}
	}
	component_verts.insert(to_add_in_component.begin(), to_add_in_component.end());
	return hadSomeCollapse;
}

void
VariableWidthContouring::
clipAndShaveCollapsedParts(Components & components, bool noShaving) {
	initializeTargets(/*walk_on_fire*/true);
	list<EdgeIterator> shavingPropagations;

	vector<MATvert *> to_add_in_component;
	for ( Component & component : components ) {
		setCurrentComponent(component);
		for ( MATvert * vert : component.vertexSet ) {
			if ( ! vert->is_collapsed() ) continue;
			// Check for Normal neighbors for the clipping process.
			for( EdgeIterator edge = vert->edges_.begin(); edge != vert->edges_.end(); ++edge ) {
				MATvert * neighbor = edge->to();
				if( ! neighbor->is_normal() ) continue;
				// Compute |nextEdge| which support the ToTheRight boundary circle
				// out of which the collapsed part stems.
				EdgeIterator nextEdge = edge->nextNormalOrTwin(); // could be any edge whose from() is |neighbor|
				assert(nextEdge->from()->is_normal());
				// we don't assert(nextEdge->to()->is_normal()) because
				// it is possible that edge->to() == nextEdge->from() is a lonely
				// degree 0 Normal vertex. This happens for example in test-6.txt
				// --scaling 2.001
				MATvert * normalVert = edge->to();
				double targetRadius = normalVert->circumcircle.radius_; // no min() !
				Disk disk(normalVert->circumcircle.center_, targetRadius);
				Target clipTarget(disk, targetRadius);
				clipTarget.edgeSource = nextEdge;
				clipTarget.isClipping = true;
				clipTarget.clipDirection = edge->from();
				clip(edge->twin(), clipTarget, to_add_in_component);
				if( debug ) cerr << "Added one clipping target at " << disk.center_ << endl;
			}
		}
		component.vertexSet.insert(to_add_in_component.begin(), to_add_in_component.end());
		to_add_in_component.clear();
		if( noShaving ) continue;
		for ( MATvert * vert : component.vertexSet ) {
			//debug = (vert->pos().x() > -27)&&(vert->pos().x() < -26) &&
			//	(vert->pos().y() > 21)&&(vert->pos().y() < 22);
			if ( ! vert->is_collapsed() ) continue;
			if( debug ) cerr << "Looking at " << vert->circumcircle.center_ << " for targets... ";
			if( ! vert->targets_.empty() ) {
				if( vert->targets_[0].isClipping ) continue;
				if( debug ) cerr << "have some." << endl;
				EdgeIterator freeway;
				const int deg = degree(vert, freeway, /*walk_on_fire*/true);
				if(debug && 1 != deg) {
					cerr << "!!! " << deg << " at " << vert->pos() << endl;
				}
				assert(1 == deg);
				shavingPropagations.push_front(freeway);
			} else {
				if( debug ) cerr << "have None." << endl;
			}
		}
		shaveLeafs(shavingPropagations);
	}
	if( debug ) cerr << endl;
	// clean remaining non-clipping targets on degree-2 collapsed vertices
	for( MATvert & vert : graph.verts) {
		Targets & t = vert.targets_;
		if( t.empty() ) continue;
		if( t[0].isClipping ) {
			/*assert(1 == t.size());*/
			if( 1 != t.size() ) {
				cerr << "WEIRD THERE! vertex at " <<
				vert.pos() << " has " << t.size() << " targets with one clipping, has normal-degree "
				<< degree(&vert)
				<< endl;
			}
			continue;
		}
		t.clear();
	}
}

void
VariableWidthContouring::
clip(EdgeIterator edge, const Target & clippingTarget, vector<MATvert *> & to_add_in_component) {
	//bool debug = true;
	struct EIComp {
		bool operator()(const EdgeIterator & left, const EdgeIterator & right) const {
			return &(*left) < &(*right);
		};
	};
	set<EdgeIterator,EIComp> props;
	vector<Target> vt = {clippingTarget};
	const Disk & clipDisk = clippingTarget.medialAxisDisk;
	while( true ) {
		MATvert * from = edge->from();
		assert(from->is_normal() || from->is_shaved());
		MATvert * to = edge->to();
		assert(to->is_collapsed());
		bool clipped;
		Disk newDisk;
		trimArcToTargetForCollapse(*edge, vt, newDisk, clipped);
		if( debug ) {
			cerr << "CLIPPING FROM(" << from->pos() << ",R=" << from->radius() << ')';
		}
		if( clipped ) {
			if( debug ) {
				cerr << "CLIPPED At(" << newDisk.center_ << ",R=" << newDisk.radius_ << ')';
			}
			MATvert * inserted_vert;
			if( (clipDisk.center_ - edge->to()->pos()).squaredLength() < 1e-8 ) {
				inserted_vert = edge->to();
				if( debug ) cerr << "Forgetting about " << inserted_vert->pos() << endl;
			} else {
				inserted_vert = & graph.insert(nullptr, edge, newDisk);
				to_add_in_component.push_back(inserted_vert);
			}
			if( 10 == degree(inserted_vert, true) )
				inserted_vert->shave();
			else
				inserted_vert->collapse();
			inserted_vert->targets_.clear(); // in case it was an already existing vertex
			inserted_vert->targets_.push_back(clippingTarget);
		} else { // propagate
			if( debug ) {
				cerr << "CLIP-PROPAG At(" << to->pos() << ",R=" << to->radius() << ')';
			}
			to->shave();
			for ( EdgeIterator e = to->edges_.begin(); e != to->edges_.end(); ++e ) {
				if( ! e->to()->is_collapsed() ) continue; // already collapsed
				if( e == edge->twin() ) continue; // don't go back!
				props.insert(e);
			}
		}
		if( debug ) cerr << endl;
		if( props.empty() ) break;
		edge = *props.begin(); props.erase(props.begin());
	}
}

void
VariableWidthContouring::
shaveLeafs(list<EdgeIterator> & propagations) {
	if( debug ) {
		cerr << "We have " << propagations.size() << " propagations.\n";
		for( auto & l : propagations ) cerr << l->from()->circumcircle.center_ << " -- ";
		cerr << endl;
	}
	while( ! propagations.empty() ) {
		EdgeIterator edge = propagations.back(); propagations.pop_back();
		MATvert * vert = edge->from();
		MATvert * neighbor = edge->to();
		//bool debug = (vert->pos().x() > 21)&&(vert->pos().x() < 22) &&
		//		(vert->pos().y() > 1000)&&(vert->pos().y() < 1001);
		if( debug ) {
			cerr << "\n Shaving " << vert->pos() << "[R=" << vert->radius() << "]-->" << neighbor->pos()
				<< " (" << vert->status << ", " << vert->targets_.size() << " targets): ";
		}
		if ( vert->targets_.empty() ) { // Sanity check
			cerr << "NO TARGET!" << endl;
			continue;
		}
		if( ! vert->is_collapsed() ) {
			cerr << "COLLAPSED OOPS" << endl;
			exit(-1);
		}
		EdgeIterator freeway;
		if( 0 == statusDegree(vert, MatStatus::Collapsed, freeway) ) continue;

		bool shaved(false);
		Disk shaveDisk;
		const Target * firstClearedTarget;
		firstClearedTarget = trimArcToTargetForCollapse(*edge, vert->targets_, shaveDisk, shaved);
		if( shaved ) { // we succesfully found a position along the branch |edge| where we can shave the branch
			if( debug ) {
				cerr << "ShavedAt(" << shaveDisk.center_ << ",R=" << shaveDisk.radius_ << ')';
			}
			MATvert * inserted_vert;
			if( (shaveDisk.center_ - edge->to()->pos()).squaredLength() < 1e-8 ) {
				inserted_vert = edge->to();
			} else {
				inserted_vert = & graph.insert(currentComponent_, edge, shaveDisk);
			}
			inserted_vert->collapse();
			MATvert * v= firstClearedTarget->vertexSource;
			if( firstClearedTarget->medialAxisDisk.radius_ == firstClearedTarget->radius ) {
				v->fullyCleared_ = true;
				if( debug ) cerr << "[ShaveFullyCleared at "<<firstClearedTarget->medialAxisDisk.center_<<"]";
			}
			vert->shave();
		} else {
			if( debug ) cerr << "NOTShaved";
			if( degree(vert, /*walk_on_fire*/true) > 0 ) { vert->shave(); }
			if( ! neighbor->is_collapsed() ) goto clean_targets;
			if( debug ) cerr << "--nei is collapsed, ";
			EdgeIterator freeway;
			int neighbor_degree = statusDegree(neighbor, MatStatus::Collapsed, freeway);
			if( neighbor_degree >= 1 ) {
				Targets & nt = neighbor->targets_;
				nt.insert(nt.end(), vert->targets_.begin(), vert->targets_.end());
			}
			if( neighbor_degree == 1 ) {
				if( debug ) cerr << "--propagate";
				propagations.push_back(freeway);
			} else if( neighbor_degree == 0 ) {
				if( debug ) cerr << "--shave nei0?";
				for( const auto & t : neighbor->targets_ ) {
					if( ! t.isClipping ) continue;
					if( debug ) cerr << "yes";
					neighbor->shave();
					break;
				}
				if( debug && neighbor->is_collapsed() ) cerr << "no";
			}
		}
clean_targets:
		vert->targets_.clear();
		if( debug ) cerr << "-del_all_tgt" << endl;
	}
}

// Used for print-path computation with CGAL
void
VariableWidthContouring::
removeAllBranches() {
	list<MATvert *> leafs;
	for( MATvert& vert_ : graph.verts) {
		if( degree(& vert_) == 1 )
			leafs.push_back(& vert_ );
	}
	while( ! leafs.empty() ) {
		MATvert * vert = leafs.back();
		leafs.pop_back();
		EdgeIterator freeway;
		int deg = degree(vert, freeway); // degree == 1 in general so freeway != nullptr
		if( 1 < deg ) {
			cerr << "ERROR: having a leaf with degree " << deg << " > 1. BUG." << endl;
			exit(-1);
		}
		if( 0 == deg ) { // It means the neighbor leaf was shrunk and died.
			vert->trim();
			continue;
		}
		MATvert * neighbor = freeway->to();
		vert->trim();
		if( ! neighbor->is_destroyed() && degree(neighbor) == 1 ) {
			leafs.push_back(neighbor);
		}
	}
	graph.removeDestroyed();
}
