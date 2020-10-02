#include <boost/test/unit_test.hpp>
#include <numeric>

#include "Utils/BiCubicInterpolate.hh"
#include "MatVec/Matrix.hh"
#include "Utils/StdVector.hh"

using namespace CoupledField;

BOOST_AUTO_TEST_CASE(bspline_test)
{
  unsigned int m = 11;
  unsigned int n = 11;

  // generate equally spaced vectors (0,1,2,...,m-1)
  StdVector<double> x(m);
  std::iota(x.Begin(), x.End(), 0);
  StdVector<double> y(n);
  std::iota(y.Begin(), y.End(), 0);

  // generate linear data
  StdVector<double> data(m*n);
  for(unsigned int i = 0; i < m; ++i) {
    for(unsigned int j = 0; j < n; ++j) {
      data[j*m + i] = i+j;
    }
  }

  // construct interpolator
  BiCubicInterpolate* bci = new BiCubicInterpolate(data, x, y, false);
  // interpolate
  bci->CalcApproximation();

  // evaluate and check some example data
  double val;
  val = bci->EvaluateFunc(2, 5);
  BOOST_TEST(val == 2 + 5);
  val = bci->EvaluateFunc(7.1, 1.4);
  BOOST_TEST(val == 7.1 + 1.4);
  val = bci->EvaluateFunc(std::sqrt(2), M_PI);
  BOOST_TEST(val == std::sqrt(2) + M_PI);

  // generate nonlinear data of periodic function
  for(unsigned int i = 0; i < m; ++i) {
    for(unsigned int j = 0; j < n; ++j) {
      data[j*m + i] = std::sin(i / (double)m * 2 * M_PI);
    }
  }

  // construct interpolator for periodic function
  BiCubicInterpolate* bci2 = new BiCubicInterpolate(data, x, y, true);
  // interpolate
  bci2->CalcApproximation();

  // evaluate and check some example data
  val = bci2->EvaluateFunc(2, 5);
  BOOST_TEST( std::abs(val - std::sin(2 / (double)m * 2 * M_PI)) < 1e-2);
  val = bci2->EvaluateFunc(7.1, 1.4);
  BOOST_TEST( std::abs(val == std::sin(7.1 / (double)m * 2 * M_PI)) < 1e-2);
  val = bci2->EvaluateFunc(std::sqrt(2), M_PI);
  BOOST_TEST( std::abs(val == std::sin(std::sqrt(2) / (double)m * 2 * M_PI)) < 1e-2);
}
