#define BOOST_TEST_MODULE cfs unittest
#include <boost/test/included/unit_test.hpp>
#include "boost/date_time/posix_time/posix_time.hpp"
#include "MatVec/Vector.hh"

#include "def_use_openmp.hh"
#ifdef USE_OPENMP
#include <omp.h>
#endif

// https://kuniganotas.wordpress.com/2011/01/14/measuring-time-with-boost-library/
typedef boost::posix_time::ptime Time;
typedef boost::posix_time::time_duration TimeDuration;
//using boost::posix_time::microsec_clock::local_time;

using boost::posix_time::microsec_clock;
using boost::posix_time::time_duration;


using namespace CoupledField;

// this has been initialized with plain old boost 1.58
// http://www.boost.org/doc/libs/1_58_0/libs/test/doc/html/index.html
// current boost has different macros.

BOOST_AUTO_TEST_CASE(first_test)
{
  int i = 1;
  BOOST_CHECK(i);
//  BOOST_CHECK(i == 2);
}

BOOST_AUTO_TEST_CASE(vector_compare)
{
  Vector<double> a(2);
  a[0] = 1;
  a[1] = 2;

  Vector<double> b(2);
  b[0] = 1;
  b[1] = 3;

  BOOST_CHECK(a == a);
  BOOST_CHECK(!(a == b));
  BOOST_CHECK(a != b);

  Vector<Complex> c(2);
  c[0] = Complex(1,2);
  c[1] = Complex(2,1);

  Vector<Complex> d = c;
  BOOST_CHECK(c == d);
  BOOST_CHECK(!(c != d));
}

BOOST_AUTO_TEST_CASE(push_back_serial)
{
  StdVector<double> a,b;
  size = 1e9;
  Time start(boost::posix_time::microsec_clock::local_time());

  for (int i = 0; i < size; i ++)
    a.Push_back(i);

  Time end(boost::posix_time::microsec_clock::local_time());
  TimeDuration dt = end - start;

  std::cout << "serial time:" << dt.total_seconds() << std::endl;
  #pragma omp parallel for critical
  for (int i = 0; i < size; i ++)
    b.Push_back(i);
}
