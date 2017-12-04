#define BOOST_TEST_MODULE cfs unittest
#include <boost/test/included/unit_test.hpp>
#include "boost/date_time/posix_time/posix_time.hpp"
#include "MatVec/Vector.hh"

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
}



BOOST_AUTO_TEST_CASE(vector_compare_performance)
{
  const int ops = 1e8;

  for(unsigned size = 3; size < 1e7; size *= 100)
  {
    Vector<double> a(size);
    for(unsigned v = 0; v < a.GetSize(); v++)
      a[v] = v;
    Vector<double> b = a;

    Time t1(boost::posix_time::microsec_clock::local_time());
    for(unsigned int i = 0, n = ops/size; i < n; i++)
    {
      a[0] = i; // make sure the compiler does not cheat
      b[0] = i;
      if(!(a == b))
        std::cout << "shall not happen\n";
    }
    Time t2(boost::posix_time::microsec_clock::local_time());

    Time t3(boost::posix_time::microsec_clock::local_time());
    for(unsigned int i = 0, n = ops/size; i < n; i++)
    {
      a[0] = i;
      b[0] = i;
      if(a != b)
        std::cout << "shall not happen\n";
    }
    Time t4(boost::posix_time::microsec_clock::local_time());

    TimeDuration dt1 = t2 - t1;
    TimeDuration dt2 = t4 - t3;
    std::cout << "vector compare: size=" << size << " n=" << (ops/size) << " opt='==' dt=" << dt1 << std::endl;
    std::cout << "vector compare: size=" << size << " n=" << (ops/size) << " opt='!=' dt=" << dt2 << std::endl;
  }

}
