#ifndef VARYFILL_FOR_ICESL_
#define VARYFILL_FOR_ICESL_

#include "VariableWidthContouring.h"

#include <polyclipping/clipper.hpp>

using CLPoint = ClipperLib::IntPoint;
using CLPath = ClipperLib::Path;
using CLPaths = ClipperLib::Paths;

class VariableWidthContouringForIceSL {
protected:
    // various parameters
    bool debug_;
    bool verbose_;

    // work parameters
    double minBeadWidth_, maxBeadWidth_, simplificationThreshold_;

    // work data
    MATGraph mat_;
    BoundaryCircles maoi_[2];

public:
    VariableWidthContouringForIceSL();
    void go(const CLPaths & in, CLPaths & out, int numberOfBeads);
    void setDebug(bool d) { debug_ = d; }
    void setVerbose(bool d) { verbose_ = d; }
};

#endif // VARYFILL_FOR_ICESL_
