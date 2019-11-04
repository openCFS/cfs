#include <boost/test/unit_test.hpp>
#include "Optimization/Design/SplineBoxDesign.hh"
#include "Utils/BSpline.hh"
#include "MatVec/Matrix.hh"

using namespace CoupledField;

BOOST_AUTO_TEST_CASE(bspline_test)
{
  StdVector<double> ref;
  Matrix<double> Ref;

  StdVector<double> knot_range;
  knot_range.Resize(2);
  knot_range[0] = 0.0;
  knot_range[1] = 10.0;
  BSpline bspline(3, 8, knot_range);

  StdVector<double> ga = bspline.GrevilleAbscissae();
  ref.Resize(4);
  ref[0] = 0; ref[1] = 10/3.; ref[2] = 20/3.; ref[3] = 10;
  // BOOST_TEST(ga == ref); // only works in boost version > 1.66
  for(unsigned int i = 0; i < 4; ++i) {
    BOOST_TEST((ga[i] == ref[i]));
  }

  StdVector<double> b(3);
  b[0] = 0; b[1] = 5; b[2] = 10;
  Matrix<double> basis = bspline.Eval(&b);
  ref.Resize(3);
  ref[0] = 1; ref[1] = 1/8.; ref[2] = 0;
  for(unsigned int i = 0; i < 3; ++i) {
    BOOST_TEST((std::abs(basis[i][0] - ref[i]) < 1e-14));
  }
  ref[0] = 0; ref[1] = 3/8.; ref[2] = 0;
  for(unsigned int i = 0; i < 3; ++i) {
    BOOST_TEST((std::abs(basis[i][1] - ref[i]) < 1e-14));
  }
  for(unsigned int i = 0; i < 3; ++i) {
    BOOST_TEST((std::abs(basis[i][2] - ref[i]) < 1e-14));
  }
  ref[0] = 0; ref[1] = 1/8.; ref[2] = 1;
  for(unsigned int i = 0; i < 3; ++i) {
    BOOST_TEST((std::abs(basis[i][3] - ref[i]) < 1e-14));
  }

  BSplineCurve curve(3,4);
  Point p;
  p = Point(1,0,1/3.);
  curve.SetControlPoint(1, p);
  p = Point(1,1,2/3.);
  curve.SetControlPoint(2, p);
  p = Point(0,1,1);
  curve.SetControlPoint(3, p);

  StdVector<double> t(3);
  t[0] = 0; t[1] = .5; t[2] = 1;
  Matrix<double> mtx = curve.Eval(&t);
  Ref.Resize(3,3);
  Ref[0][0] = 0; Ref[0][1] = 0; Ref[0][2] = 0;
  Ref[1][0] = 3/4.; Ref[1][1] = .5; Ref[1][2] = .5;
  Ref[2][0] = 0; Ref[2][1] = 1; Ref[2][2] = 1;
  for(unsigned int i = 0; i < 3; ++i) {
    for(unsigned int d = 0; d < 3; ++d) {
      BOOST_TEST((std::abs(mtx[i][d] - Ref[i][d]) < 1e-14));
    }
  }
}
