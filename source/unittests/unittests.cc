// define name "cfs unittest" of test module
#define BOOST_TEST_MODULE cfs unittest
//defines test framework including function main,
// wich will call the subsequently defined test cases
#include <boost/test/included/unit_test.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/mpl/list.hpp>

#include "MatVec/Vector.hh"

/* Example unit tests for CFS.
 *
 * List all available tests with
 * ./cfstestbin --list_content
 *
 * Run tests and see errors and std output:
 * ./cfstestbin --color_output --log_level=message
 *
 * More output with runtimes of individual test cases:
 * ./cfstestbin --color_output --log_level=unit_scope
 *
 * Even more output with passed tests:
 * ./cfstestbin --color_output --log_level=success
 */

// https://kuniganotas.wordpress.com/2011/01/14/measuring-time-with-boost-library/
typedef boost::posix_time::ptime Time;
typedef boost::posix_time::time_duration TimeDuration;
//using boost::posix_time::microsec_clock::local_time;

using boost::posix_time::microsec_clock;
using boost::posix_time::time_duration;


using namespace CoupledField;


// BOOST_AUTO_TEST_CASE declares a test case named first_test, which in turn
// will run the content of first_test inside the controlled testing environment
BOOST_AUTO_TEST_CASE(first_test)
{
  int i = 1;
  BOOST_TEST(i); // passes
  //BOOST_TEST(i == 2); // fails
}



// there is an (fixed) issue that clang on mac cannot catch exceptions
BOOST_AUTO_TEST_CASE(exception)
{
  try
  {
    throw Exception("test");
  }
  catch(...)
  {
    BOOST_TEST(true);
    return;
  }
  BOOST_TEST(false);
}


BOOST_AUTO_TEST_CASE(vector_compare)
{
  Vector<double> a(2);
  a[0] = 1;
  a[1] = 2;

  Vector<double> b(2);
  b[0] = 1;
  b[1] = 3;

  BOOST_TEST(a == a);
  BOOST_TEST(!(a == b));
  BOOST_TEST(a != b);

  Vector<Complex> c(2);
  c[0] = Complex(1,2);
  c[1] = Complex(2,1);

  Vector<Complex> d = c;
  BOOST_TEST(c == d);
  BOOST_TEST(!(c != d));
}

// example for performance test
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

// example for template test
typedef boost::mpl::list<int,long,unsigned char> test_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( my_test, T, test_types )
{
  BOOST_TEST( sizeof(T) == (unsigned)4 );
}
