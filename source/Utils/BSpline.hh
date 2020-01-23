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

  BSpline(const unsigned int degree, unsigned int numKnots, double min_knot = 0.0, double max_knot = 1.0);

  StdVector<double> GetKnots() { return this->knots_; };

  /** Calculate Greville abscissae */
  StdVector<double> GrevilleAbscissae();

  /** Values of all basis functions at points t
   * @return Matrix<double> of size t.GetSize() x number of basis functions
   */
  Matrix<double> Eval(StdVector<double>* t);

  /** For convenience. Same as above */
  Matrix<double> Eval(double t);

private:

  /** Generates a knot vector of equally spaced knots in the interval (0,1),
   * where the first and last knots are repeated degree+1 times (generates
   * clamped BSpline curve). */
  void MakeKnotVector(double min_knot = 0.0, double max_knot = 1.0);

  /** Generates all BSpline basis functions of degree zero evaluated at t
   * @return Matrix<double> of size t.GetSize() x numIntervals
   */
  Matrix<double> BasisFuncDg0(StdVector<double>* t, StdVector<double> knots, unsigned int numIntervals);

  /** Generates all BSpline basis functions of a given degree evaluated at t
   * with Cox-de Boor recursion formula*/
  Matrix<double> BasisFuncDg(unsigned int degree, StdVector<double>* t, StdVector<double> knots, unsigned int numIntervals);

  unsigned int degree_;

  unsigned int numKnots_;

  StdVector<double> knots_;

};


class BSplineCurve // : public BSpline
{
public:

  BSplineCurve(unsigned int degree, unsigned int numControlPoints);

  void SetControlPoint(int idx, Point coords);

  /** Set control points of the BSpline curve */
  void SetControlPoints(StdVector<Point> coords);

  /** Returns the control points of the BSpline curve */
  StdVector<Point> GetControlPoints() { return this->control_points_; }

  Matrix<double> Eval(StdVector<double>* t);

private:

  BSpline* bspline_;

  /** number of control points */
  unsigned int numCP_;

  /** vector of control points */
  StdVector<Point> control_points_;

};

} //end of namespace

#endif //BSPLINE_HH_
