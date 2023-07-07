#include <boost/test/unit_test.hpp>
#include <numeric>
#include <random>

#include "Utils/CubicInterpolate.hh"
#include "Utils/BiCubicInterpolate.hh"
#include "Utils/TriCubicInterpolate.hh"
#include "MatVec/Matrix.hh"
#include "Utils/StdVector.hh"

using namespace CoupledField;
namespace tt = boost::test_tools;

BOOST_AUTO_TEST_CASE(cubic_interpolation_test)
{
  StdVector<double> x;
  StdVector<double> y;
  StdVector<double> z;

  // vectors generated with random data
  x = {0.937834, 1.11794, 1.13992, 1.23913, 1.78653, 2.0868, 2.25648, 2.39954, 2.77077, 3.56925, 3.73325, 3.97086, 4.4773, 4.5256, 4.79547, 4.82073, 4.96648, 5.30127, 5.40927, 5.41837, 5.51375, 5.53045, 5.5902, 5.61671, 5.63746, 6.08254, 6.27103, 6.29292, 6.31776, 6.32885, 6.95834, 7.00418, 7.24709, 7.36147, 7.62384, 7.63857, 7.67445, 7.72959, 7.938, 7.95097, 8.32739, 8.52808, 8.5491, 8.60572, 8.75395, 9.35709, 9.47514, 9.57827, 9.7194};
  y = {2.04151, 2.19219, 2.40335, 2.45827, 2.4829, 2.51386, 2.55797, 2.69608, 2.90534, 3.00021, 3.10684, 3.21089, 3.42914, 3.46223, 3.52802, 3.57174, 3.5836, 3.64165, 3.9216, 3.9313, 3.97504, 4.1332, 4.22177, 4.24356, 4.37898, 4.38541, 4.4113, 4.67886, 4.73698, 4.76531, 4.86669, 4.94559, 5.16865, 5.28212, 5.44558, 5.53435, 5.56211, 5.5969, 5.72705, 5.98373, 6.00778, 6.10382, 6.14923, 6.28648, 6.50166, 6.65614, 6.72795};
  z = {2.79877, 3.05751, 3.68023, 3.91935, 3.93835, 4.14635, 4.3491, 4.37901, 4.64531, 4.68799, 5.06096, 5.4016, 5.75909, 6.03477, 6.07384, 6.37583, 6.45982, 6.66914, 6.86211, 8.77086, 9.19177, 9.72735, 9.78728, 9.87244, 10.2095, 10.2986, 10.7067, 10.7292, 10.809, 10.9975, 11.5525, 11.718, 12.3462, 13.1937, 13.4069, 13.7091, 13.9965, 14.1317, 14.1505, 14.8051, 15.111, 15.5572, 15.9827, 16.0957, 16.7998};

  unsigned int m = x.GetSize(); // = 49
  unsigned int n = y.GetSize(); // = 47
  unsigned int o = z.GetSize(); // = 45

  //////////////////////////////////// 1D ////////////////////////////////////
  // generate linear data
  StdVector<double> datac(m);
  for (unsigned int i = 0; i < m; ++i)
    datac[i] = 2.*x[i];

  // construct interpolator
  CubicInterpolate* ci = new CubicInterpolate(datac, x, false);
  // interpolate
  ci->CalcApproximation();

  // evaluate and check some example data
  // interpolation of linear data has to be exact
  double val;
  val = ci->EvaluateFunc(2);
  BOOST_TEST(val == 2.*2., tt::tolerance(1e-12));
  val = ci->EvaluateFunc(7.1);
  BOOST_TEST(val == 2.*7.1, tt::tolerance(1e-12));
  val = ci->EvaluateFunc(std::sqrt(2));
  BOOST_TEST(val == 2.*std::sqrt(2), tt::tolerance(1e-12));

  double dvalc;
  dvalc = ci->EvaluateDeriv(2);
  BOOST_TEST(dvalc == 2., tt::tolerance(1e-12));
  dvalc = ci->EvaluateDeriv(7.1);
  BOOST_TEST(dvalc == 2., tt::tolerance(1e-12));
  dvalc = ci->EvaluateDeriv(std::sqrt(2));
  BOOST_TEST(dvalc == 2., tt::tolerance(1e-12));

  double ddvalc;
  ddvalc = ci->EvaluateSecondDeriv(2);
  BOOST_TEST(ddvalc == 0., tt::tolerance(1e-12));
  ddvalc = ci->EvaluateSecondDeriv(7.1);
  BOOST_TEST(ddvalc == 0., tt::tolerance(1e-12));
  ddvalc = ci->EvaluateSecondDeriv(std::sqrt(2));
  BOOST_TEST(ddvalc == 0., tt::tolerance(1e-12));

  // generate nonlinear data of periodic function
  for (unsigned int i = 0; i < m; ++i)
    datac[i] = std::sin(x[i] / (x[m-1]-x[0]) * 2 * M_PI);

  // construct interpolator for periodic function
  CubicInterpolate* ci2 = new CubicInterpolate(datac, x, true);
  // interpolate
  ci2->CalcApproximation();

  // evaluate and check some example data
  val = ci2->EvaluateFunc(0.937834); // interpolation at data point has to be exact
  BOOST_TEST(val == std::sin(0.937834 / (x[m-1]-x[0]) * 2 * M_PI), tt::tolerance(1e-12));
  val = ci2->EvaluateFunc(2);
  BOOST_TEST(val == std::sin(2 / (x[m-1]-x[0]) * 2 * M_PI), tt::tolerance(0.0005));
  val = ci2->EvaluateFunc(-7.1);
  BOOST_TEST(val == std::sin(-7.1 / (x[m-1]-x[0]) * 2 * M_PI), tt::tolerance(0.006));
  val = ci2->EvaluateFunc(std::sqrt(2));
  BOOST_TEST(val == std::sin(std::sqrt(2) / (x[m-1]-x[0]) * 2 * M_PI), tt::tolerance(0.015));

  dvalc = ci2->EvaluateDeriv(2);
  BOOST_TEST(dvalc == std::cos(2 / (x[m-1]-x[0]) * 2 * M_PI) / (x[m-1]-x[0]) * 2 * M_PI, tt::tolerance(0.2));
  dvalc = ci2->EvaluateDeriv(-7.1);
  BOOST_TEST(dvalc == std::cos(-7.1 / (x[m-1]-x[0]) * 2 * M_PI) / (x[m-1]-x[0]) * 2 * M_PI, tt::tolerance(0.2));
  dvalc = ci2->EvaluateDeriv(std::sqrt(2));
  BOOST_TEST(dvalc == std::cos(std::sqrt(2) / (x[m-1]-x[0]) * 2 * M_PI) / (x[m-1]-x[0]) * 2 * M_PI, tt::tolerance(0.2));


  //////////////////////////////////// 2D ////////////////////////////////////

  // generate linear data
  StdVector<double> databc(m*n);
  for (unsigned int i = 0; i < m; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      databc[j*m + i] = 2.*x[i] + y[j];
    }
  }

  // construct interpolator
  BiCubicInterpolate* bci = new BiCubicInterpolate(databc, x, y, false);
  // interpolate
  bci->CalcApproximation();

  // evaluate and check some example data
  // linear interpolation has to be exact
  val = bci->EvaluateFunc(2, 5);
  BOOST_TEST(val == 2.*2. + 5., tt::tolerance(1e-12));
  val = bci->EvaluateFunc(7.1, 4.0);
  BOOST_TEST(val == 2.*7.1 + 4.0, tt::tolerance(1e-12));
  val = bci->EvaluateFunc(std::sqrt(20), M_PI);
  BOOST_TEST(val == 2.*std::sqrt(20) + M_PI, tt::tolerance(1e-12));

  Vector<double> dvalb(2, 0);
  Vector<double> refb(2);
  refb[0] = 2.;
  refb[1] = 1.;
  dvalb = bci->EvaluatePrime(2, 5);
  // tt::per_element() does not work, because there is no elem_op for Vector
  for (unsigned int i = 0; i < refb.GetSize(); ++i)
    BOOST_TEST(dvalb[i] == refb[i], tt::tolerance(1e-10));
  dvalb = bci->EvaluatePrime(7.1, 4.0);
  for (unsigned int i = 0; i < refb.GetSize(); ++i)
    BOOST_TEST(dvalb[i] == refb[i], tt::tolerance(1e-10));
  dvalb = bci->EvaluatePrime(std::sqrt(20), M_PI);
  for (unsigned int i = 0; i < refb.GetSize(); ++i)
    BOOST_TEST(dvalb[i] == refb[i], tt::tolerance(1e-10));

  // generate nonlinear data of periodic function
  for (unsigned int i = 0; i < m; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      databc[j*m + i] = std::sin(x[i] / (x[m-1]-x[0]) * 2 * M_PI) + std::cos(y[j] / (y[n-1]-y[0]) * 2 * M_PI);
    }
  }

  // construct interpolator for periodic function
  BiCubicInterpolate* bci2 = new BiCubicInterpolate(databc, x, y, true);
  // interpolate
  bci2->CalcApproximation();

  // evaluate and check some example data
  val = bci2->EvaluateFunc(0.937834, 2.19219); // interpolation at data point has to be exact
  BOOST_TEST(val == std::sin(0.937834 / (x[m-1]-x[0]) * 2 * M_PI) + std::cos(2.19219 / (y[n-1]-y[0]) * 2 * M_PI), tt::tolerance(1e-12));
  val = bci2->EvaluateFunc(2, 5);
  BOOST_TEST(val == std::sin(2 / (x[m-1]-x[0]) * 2 * M_PI) + std::cos(5 / (y[n-1]-y[0]) * 2 * M_PI), tt::tolerance(0.003));
  val = bci2->EvaluateFunc(-7.1, 0.0);
  BOOST_TEST(val == std::sin(-7.1 / (x[m-1]-x[0]) * 2 * M_PI) + std::cos(0.0 / (y[n-1]-y[0]) * 2 * M_PI), tt::tolerance(0.003));
  val = bci2->EvaluateFunc(std::sqrt(20), M_PI);
  BOOST_TEST(val == std::sin(std::sqrt(20) / (x[m-1]-x[0]) * 2 * M_PI) + std::cos(M_PI / (y[n-1]-y[0]) * 2 * M_PI), tt::tolerance(0.003));

  dvalb = bci2->EvaluatePrime(2, 4);
  BOOST_TEST(dvalb[0] ==  std::cos(2 / (x[m-1]-x[0]) * 2 * M_PI) / (x[m-1]-x[0]) * 2 * M_PI, tt::tolerance(0.2));
  BOOST_TEST(dvalb[1] == -std::sin(4 / (y[n-1]-y[0]) * 2 * M_PI) / (y[n-1]-y[0]) * 2 * M_PI, tt::tolerance(0.04));
  dvalb = bci2->EvaluatePrime(-7.1, 0.0);
  BOOST_TEST(dvalb[0] ==  std::cos(-7.1 / (x[m-1]-x[0]) * 2 * M_PI) / (x[m-1]-x[0]) * 2 * M_PI, tt::tolerance(0.2));
  BOOST_TEST(dvalb[1] == -std::sin( 0.0 / (y[n-1]-y[0]) * 2 * M_PI) / (y[n-1]-y[0]) * 2 * M_PI, tt::tolerance(0.1));
  dvalb = bci2->EvaluatePrime(std::sqrt(20), M_PI);
  BOOST_TEST(dvalb[0] ==  std::cos(std::sqrt(20) / (x[m-1]-x[0]) * 2 * M_PI) / (x[m-1]-x[0]) * 2 * M_PI, tt::tolerance(0.01));
  BOOST_TEST(dvalb[1] == -std::sin(M_PI          / (y[n-1]-y[0]) * 2 * M_PI) / (y[n-1]-y[0]) * 2 * M_PI, tt::tolerance(0.01));

  //////////////////////////////////// 3D ////////////////////////////////////

  // generate linear data
  StdVector<double> datatc(m*n*o);
  for (unsigned int i = 0; i < m; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      for (unsigned int k = 0; k < o; ++k) {
        datatc[k*n*m + j*m + i] = 2.*x[i] + y[j] + 3.*z[k];
      }
    }
  }

  // construct interpolator
  TriCubicInterpolate* tci = new TriCubicInterpolate(datatc, x, y, z, false);
  // interpolate
  tci->CalcApproximation();

  // evaluate and check some example data
  // linear interpolation has to be exact
  val = tci->EvaluateFunc(2, 5, 3);
  BOOST_TEST(val == 2.*2. + 5. + 3.*3., tt::tolerance(1e-12));
  val = tci->EvaluateFunc(7.1, 4.4, 10.0);
  BOOST_TEST(val == 2.*7.1 + 4.4 + 3.*10.0, tt::tolerance(1e-12));
  val = tci->EvaluateFunc(std::sqrt(20), M_PI, 4.5);
  BOOST_TEST(val == 2.*std::sqrt(20) + M_PI + 3.*4.5, tt::tolerance(1e-12));

  Vector<double> dvalt(3,0);
  Vector<double> reft(3,0);
  reft[0] = 2.;
  reft[1] = 1.;
  reft[2] = 3.;
  dvalb = tci->EvaluatePrime(2, 5, 3);
  // tt::per_element() does not work, because there is no elem_op for Vector
  for (unsigned int i = 0; i < reft.GetSize(); ++i)
    BOOST_TEST(dvalb[i] == reft[i], tt::tolerance(1e-10));
  dvalb = tci->EvaluatePrime(7.1, 4.4, 10.0);
  for (unsigned int i = 0; i < reft.GetSize(); ++i)
    BOOST_TEST(dvalb[i] == reft[i], tt::tolerance(1e-10));
  dvalb = tci->EvaluatePrime(std::sqrt(20), M_PI, 4.5);
  for (unsigned int i = 0; i < reft.GetSize(); ++i)
    BOOST_TEST(dvalb[i] == reft[i], tt::tolerance(1e-10));

  // generate nonlinear data of periodic function
  for (unsigned int i = 0; i < m; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      for (unsigned int k = 0; k < o; ++k) {
        datatc[k*n*m + j*m + i] = std::sin(x[i] / (x[m-1]-x[0]) * 2 * M_PI) + std::cos(y[j] / (y[n-1]-y[0]) * 2 * M_PI) * std::sin(z[k] / (z[o-1]-z[0]) * 2 * M_PI);
      }
    }
  }

  // construct interpolator for periodic function
  TriCubicInterpolate* tci2 = new TriCubicInterpolate(datatc, x, y, z, true);
  // interpolate
  tci2->CalcApproximation();

  // evaluate and check some example data
  val = tci2->EvaluateFunc(0.937834, 2.19219, 3.68023); // interpolation at data point has to be exact
  BOOST_TEST(val == std::sin(0.937834 / (x[m-1]-x[0]) * 2 * M_PI)
      + std::cos(2.19219 / (y[n-1]-y[0]) * 2 * M_PI) * std::sin(3.68023 / (z[o-1]-z[0]) * 2 * M_PI), tt::tolerance(1e-12));
  val = tci2->EvaluateFunc(2, 5, 3);
  BOOST_TEST(val == std::sin(2 / (x[m-1]-x[0]) * 2 * M_PI)
      + std::cos(5 / (y[n-1]-y[0]) * 2 * M_PI) * std::sin(3 / (z[o-1]-z[0]) * 2 * M_PI), tt::tolerance(0.002));
  val = tci2->EvaluateFunc(7., 1., 0.);
  BOOST_TEST(val == std::sin(7. / (x[m-1]-x[0]) * 2 * M_PI)
      + std::cos(1. / (y[n-1]-y[0]) * 2 * M_PI) * std::sin(0. / (z[o-1]-z[0]) * 2 * M_PI), tt::tolerance(0.0003));
  val = tci2->EvaluateFunc(std::sqrt(20), M_PI, 4.5);
  BOOST_TEST(val == std::sin(std::sqrt(20) / (x[m-1]-x[0]) * 2 * M_PI)
      + std::cos(M_PI / (y[n-1]-y[0]) * 2 * M_PI) * std::sin(4.5 / (z[o-1]-z[0]) * 2 * M_PI), tt::tolerance(0.0008));

  dvalt = tci2->EvaluatePrime(2, 4, 3);
  BOOST_TEST(dvalt[0] ==  std::cos(2 / (x[m-1]-x[0]) * 2 * M_PI) / (x[m-1]-x[0]) * 2 * M_PI, tt::tolerance(0.3));
  BOOST_TEST(dvalt[1] == -std::sin(4 / (y[n-1]-y[0]) * 2 * M_PI) / (y[n-1]-y[0]) * 2 * M_PI * std::sin(3 / (z[o-1]-z[0]) * 2 * M_PI), tt::tolerance(0.04));
  BOOST_TEST(dvalt[2] ==  std::cos(4 / (y[n-1]-y[0]) * 2 * M_PI) * std::cos(3 / (z[o-1]-z[0]) * 2 * M_PI) / (z[o-1]-z[0]) * 2 * M_PI, tt::tolerance(0.3));
  dvalt = tci2->EvaluatePrime(7., 1., 0.);
  BOOST_TEST(dvalt[0] ==  std::cos(7. / (x[m-1]-x[0]) * 2 * M_PI) / (x[m-1]-x[0]) * 2 * M_PI, tt::tolerance(0.3));
  BOOST_TEST(dvalt[1] == -std::sin(1. / (y[n-1]-y[0]) * 2 * M_PI) / (y[n-1]-y[0]) * 2 * M_PI * std::sin(0. / (z[o-1]-z[0]) * 2 * M_PI), tt::tolerance(1e-5));
  BOOST_TEST(dvalt[2] ==  std::cos(1. / (y[n-1]-y[0]) * 2 * M_PI) * std::cos(0. / (z[o-1]-z[0]) * 2 * M_PI) / (z[o-1]-z[0]) * 2 * M_PI, tt::tolerance(0.003));
  dvalt = tci2->EvaluatePrime(std::sqrt(20), M_PI, 4.5);
  BOOST_TEST(dvalt[0] ==  std::cos(std::sqrt(20) / (x[m-1]-x[0]) * 2 * M_PI) / (x[m-1]-x[0]) * 2 * M_PI, tt::tolerance(0.01));
  BOOST_TEST(dvalt[1] == -std::sin(M_PI / (y[n-1]-y[0]) * 2 * M_PI) / (y[n-1]-y[0]) * 2 * M_PI * std::sin(4.5 / (z[o-1]-z[0]) * 2 * M_PI), tt::tolerance(0.01));
  BOOST_TEST(dvalt[2] ==  std::cos(M_PI / (y[n-1]-y[0]) * 2 * M_PI) * std::cos(4.5 / (z[o-1]-z[0]) * 2 * M_PI) / (z[o-1]-z[0]) * 2 * M_PI, tt::tolerance(0.01));
}
