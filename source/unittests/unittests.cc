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
  int datSize = 1e9;
  int nRuns = 10;
  StdVector<double> a,b;
  TimeDuration dt1, dt2, dt3, dt4;

  for (int runs = 0; runs < nRuns; runs++) {
    a.Clear();
    b.Clear();
    Time start1(boost::posix_time::microsec_clock::local_time());
    for (int i = 0; i < datSize; i ++)
      a.Push_back(i);

    Time end1(boost::posix_time::microsec_clock::local_time());
    dt1 += end1 - start1;

    Time start2(boost::posix_time::microsec_clock::local_time());

    #pragma omp parallel num_threads(1)
    {
//      std::cout << "running with " << omp_get_num_threads() << " thread" << std::endl;
      #pragma omp for
      for (int i = 0; i < datSize; i ++) {
        #pragma omp critical
        b.Push_back(i);
      }
    }
    Time end2(boost::posix_time::microsec_clock::local_time());
    dt2 += end2 - start2;

    BOOST_CHECK(a == b);

    a.Clear();
    b.Clear();

    Time start3(boost::posix_time::microsec_clock::local_time());
    #pragma omp parallel
    {
  //    std::cout << "running with " << omp_get_num_threads() << std::endl;
      #pragma omp for
      for (int i = 0; i < datSize; i ++)
        #pragma omp critical
        a.Push_back(i);
    }
    Time end3(boost::posix_time::microsec_clock::local_time());
    dt3 += end3 - start3;

    a.Clear();

    Time start4(boost::posix_time::microsec_clock::local_time());
    #pragma omp parallel
    {
  //    std::cout << "running with " << omp_get_num_threads() << std::endl;
      #pragma omp for
      for (int i = 0; i < datSize; i ++)
        a.Push_back(i);
    }
    Time end4(boost::posix_time::microsec_clock::local_time());
    dt4 += end4 - start4;
  }

  std::cout << "Average run times for 10 loop iterations" << std::endl;
  std::cout << "serial time: " << dt1.total_seconds()/nRuns << " seconds" << std::endl;
  std::cout << "1 thread with critical :" << dt2.total_seconds()/nRuns << " seconds" << std::endl;
  std::cout << "multiple threads with critical :" << dt3.total_seconds()/nRuns << " seconds" << std::endl;
  std::cout << "multiple threads without critical :" << dt4.total_seconds()/nRuns << " seconds" << std::endl;
}
