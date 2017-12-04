#define BOOST_TEST_MODULE cfs unittest
#include <boost/test/included/unit_test.hpp>
#include "MatVec/Vector.hh"

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
