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


struct Foo
{
  int noRef() { return v; }
  int& withRef() { return v; }
  int v = 0;
};

BOOST_AUTO_TEST_CASE(reference_test)
{
  Foo foo;

  // working with a returned reference works only we have references on both sides.

  int a = foo.withRef();
  BOOST_TEST(a == 0);
  foo.v = 1;
  BOOST_TEST(a == 0);

  int& b = foo.withRef();
  b = 2;
  BOOST_TEST(foo.v == 2);

  // cannot compile: error: non-const lvalue reference to type 'int' cannot bind to a temporary of type 'int'
  // int& c = foo.noRef();

  // trivial
  int d = foo.noRef();
  if(d == 45) // we need to use the variable
    std::cout << "Hans";
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
    StdVector<TimeDuration> use;
    use.Push_back(dt1);
    use.Push_back(dt2);
    // std::cout << "vector compare: size=" << size << " n=" << (ops/size) << " opt='==' dt=" << dt1 << std::endl;
    // std::cout << "vector compare: size=" << size << " n=" << (ops/size) << " opt='!=' dt=" << dt2 << std::endl;
  }
}


/** copy function for copy_performance_test */
void loop_copy(double* source, double* target, int size)
{
  for(int i = 0; i < size; i++)
    target[i] = source[i];
}

void par_loop_copy(double* source, double* target, int size)
{
  #pragma omp parallel for num_threads(4)
  for(int i = 0; i < size; i++)
    target[i] = source[i];
}


/** measure performance of looped copy vs. std::copy */
BOOST_AUTO_TEST_CASE(loop_vs_copy)
{
   int max_p = 6;
   StdVector<Vector<double> > source(max_p + 1);
   StdVector<Vector<double> > target(max_p + 1);

   // our exponent from 8 * 10^1 to 8 * 10^p bytes
   for(int p = 1; p <= max_p; p++) {
     source[p].Resize((int) pow(10,p));
     target[p].Resize((int) pow(10,p));
   }

   long values = (long) 1e7;

   std::cout << "#set logscale y; set xlabel \"datasize in doubles 10^p\"; set ylabel \"copy " << values << " doubles in s\"\n";
   std::cout << "# plot \"loop.dat\" u 1:4 t \"loop\" w lp, \"loop.dat\" u 1:5 t \"par (4) loop\" w lp, \"loop.dat\" u 1:6 t \"std::copy\" w lp, \"loop.dat\" u 1:7 t \"std::memcopy\" w lp\n";
   std::cout << "#(1) p \t(2) n \t(3) loops \t(4) loop \t(5) par_4_loop \t(6) std::copy \t(7) std::memcpy\n";
   for(int p = 1; p <= max_p; p++) {
     Vector<double>& s = source[p];
     Vector<double>& t = target[p];
     unsigned int n = s.GetSize();

     int loops = values / n;

     Time t1(boost::posix_time::microsec_clock::local_time());
     for(int l = 0; l < loops; l++)
       loop_copy(s.GetPointer(), t.GetPointer(), n);
     Time t2(boost::posix_time::microsec_clock::local_time());

     Time t3(boost::posix_time::microsec_clock::local_time());
     for(int l = 0; l < loops; l++)
       par_loop_copy(s.GetPointer(), t.GetPointer(), n);
     Time t4(boost::posix_time::microsec_clock::local_time());


     Time t5(boost::posix_time::microsec_clock::local_time());
     for(int l = 0; l < loops; l++)
      std::copy_n(s.GetPointer(), n, t.GetPointer());
     Time t6(boost::posix_time::microsec_clock::local_time());


     Time t7(boost::posix_time::microsec_clock::local_time());
     for(int l = 0; l < loops; l++)
       std::memcpy(t.GetPointer(), s.GetPointer(), sizeof(double) * n);
     Time t8(boost::posix_time::microsec_clock::local_time());

     TimeDuration loop = t2 - t1;
     TimeDuration parl = t4 - t3;
     TimeDuration copy = t6 - t5;
     TimeDuration memcpy = t8 - t7;

     std::cout << p << " \t" << n << " \t"  << loops << " \t" <<  std::setprecision(8) <<  loop.total_microseconds() / 1.0e6
               << " \t" << parl.total_microseconds() / 1.0e6 << " \t" << copy.total_microseconds() / 1.0e6
               << " \t" << memcpy.total_microseconds() / 1.0e6 << std::endl;
   }
}



// example for template test
typedef boost::mpl::list<int, unsigned int, float> test_types; // fails for long and char
BOOST_AUTO_TEST_CASE_TEMPLATE( my_test, T, test_types )
{
  BOOST_TEST( sizeof(T) == (unsigned) 4 );
}

BOOST_AUTO_TEST_CASE(desctructor_call)
{
  int cnt = 0;
  struct Hans
  {
    Hans() {/*std::cout << "Hans()\n";*/};
    Hans(int* cnt) {/*std::cout << "Hans(int*)\n";*/Init(cnt);}
    ~Hans() {/*std::cout << "~Hans()\n";*/ if(cnt) (*cnt)--; }
    void Init(int* p){this->cnt = p; (*cnt)++;}
    int* cnt = NULL;
  };

  BOOST_TEST(cnt == 0);
  //std::cout << "test StdVector<Hans*>\n";
  // class 1 constrcutor but 0 destructors
  {
    StdVector<Hans*> test(1);
    test[0] = new Hans(&cnt);
  }
  BOOST_TEST(cnt == 1); // destructor is not called

  //std::cout << "test StdVector<Hans>\n";
  // calls 2 constrcutors(!) and 2 destrcutors. No idea where the two constructors come from?!
  {
    StdVector<Hans> test(1);
    test[0].Init(&cnt);

    BOOST_TEST(cnt == 2);
  }
  BOOST_TEST(cnt == 1); // destructor is called
}

