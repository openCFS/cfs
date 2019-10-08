/**
 * This is a C++ re-implementation of bspline.py from sfepy.
 */

#ifndef BSPLINE_HH_
#define BSPLINE_HH_

#include "Utils/Point.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {

//! \todo
//! Add support for cyclic splines
//! Add support for different knot_types (clamped|cyclic|userdef)
class BSpline
{
public:

  BSpline(const unsigned int degree, unsigned int numKnots, StdVector<double> knot_range = StdVector<double>());

  StdVector<double> GetKnots() { return this->knots; };

  /** Calculate Greville abscissae */
  StdVector<double> GrevilleAbscissae();

  Matrix<double> EvalBasis(StdVector<double>* t);

private:

  unsigned int degree;

  unsigned int numKnots;

  StdVector<double> knots;

  /** Generates a knot vector of equally spaced knots in the interval (0,1),
   * where the first and last knots are repeated degree+1 times (generates
   * clamped BSpline curve). */
  void MakeKnotVector(StdVector<double> knot_range);

  /** Generates all BSpline basis functions of degree zero evaluated at t */
  Matrix<double> BasisFuncDg0(StdVector<double>* t, StdVector<double> knots, unsigned int numIntervals);

  /** Generates all BSpline basis functions of a given degree evaluated at t
   * with Cox-de Boor recursion formula*/
  Matrix<double> BasisFuncDg(unsigned int degree, StdVector<double>* t, StdVector<double> knots, unsigned int numIntervals);

};


class BSplineCurve // : public BSpline
{
public:

  BSplineCurve(unsigned int degree, unsigned int numControlPoints);

  void SetControlPoint(int idx, Point coords);

  /** Set control points of the BSpline curve */
  void SetControlPoints(StdVector<Point> coords);

  /** Returns the control points of the BSpline curve */
  StdVector<Point> GetControlPoints() { return this->control_points; }

  Matrix<double> Eval(StdVector<double>* t);

private:

  BSpline* bspline;

  /** number of control points */
  unsigned int numCP;

  /** vector of control points */
  StdVector<Point> control_points;

};

} //end of namespace

#endif //BSPLINE_HH_
