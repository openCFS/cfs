// define name "cfs unittest" of test module
#define BOOST_TEST_MODULE cfs unittest

// defines test framework including function main,
// which will call the subsequently defined test cases

#pragma clang diagnostic push
#pragma GCC diagnostic push
// boost 1.78: boost/test/impl/debug.ipp:674:25: warning: variable 'junk' set but not used
#pragma clang diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"

#include <boost/test/included/unit_test.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/mpl/list.hpp>

#pragma clang diagnostic pop
#pragma GCC diagnostic pop

#undef max // prevent Windows issue
#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"

/* Example unit tests for CFS.
 *
 * List all available tests with
 * ./cfstest --list_content
 *
 * Run tests and see errors and std output:
 * ./cfstest --color_output --log_level=message
 *
 * More output with runtimes of individual test cases:
 * ./cfstest --color_output --log_level=ut_scope
 *
 * Even more output with passed tests:
 * ./cfstest --color_output --log_level=success
 */

// https://kuniganotas.wordpress.com/2011/01/14/measuring-time-with-boost-library/
typedef boost::posix_time::ptime Time;
typedef boost::posix_time::time_duration TimeDuration;
//using boost::posix_time::microsec_clock::local_time;

using boost::posix_time::microsec_clock;
using boost::posix_time::time_duration;


using namespace CoupledField;

//#define TOL boost::unit_test::tolerance( boost::test_tools::fpc::percent_tolerance(0.01) )

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

   std::cout << "loop_vs_copy" << std::endl;
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
     std::cout << std::endl;
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

BOOST_AUTO_TEST_CASE(small_direct_solver)
{
  // original matrix
  Matrix<double> O(3,3);
  O[0][0] = -484.829;
  O[0][1] = -484.336;
  O[0][2] = 0.13088;
  O[1][0] = -484.336;
  O[1][1] = -484.829;
  O[1][2] = 0.13088;
  O[2][0] = 0.13088;
  O[2][1] = 0.13088;
  O[2][2] = -0.0399633;

  Vector<double> b(3);
  b[0] = 0.00182783;
  b[1] = 0.00182783;
  b[2] = 0.0324977;

  Matrix<double> A(O);
  Vector<double> x; // automatic resize
  A.DirectSolve(x,b);
  std::cout << "small_direct_solver" << std::endl;
  std::cout << "A=" << A.ToString() << std::endl;
  std::cout << "b=" << b.ToString() << std::endl;
  std::cout << "x=" << x.ToString() << std::endl;
  std::cout << std::endl;
  BOOST_TEST(x.GetSize() == 3);
}

BOOST_AUTO_TEST_CASE(choleksy_lapack_solver)
{
  Matrix<double> A(2,2);
  A[0][0] = 1;
  A[0][1] = 2;
  A[1][0] = 2;
  A[1][1] = 100;

  Vector<double> b(2, 1.0); // all ones

  Matrix<double> chol;
  Vector<double> x;

  A.CholeskySolveLapack(chol,x,b,true);
  BOOST_TEST(x.GetSize() == 2);

  A[0][0] = -1;
  bool ok = A.CholeskySolveLapack(chol,x,b,false);
  BOOST_TEST(!ok);

  Matrix<double> H(3,3);
  H[0][0] = 3.2709520e+01;
  H[0][1] = 2.9368747e+01;
  H[0][2] = -8.4678084e-01;
  H[1][0] = 2.9368747e+01;
  H[1][1] = 3.2709520e+01;
  H[1][2] = -8.4678084e-01;
  H[2][0] = -8.4678084e-01;
  H[2][1] = -8.4678084e-01;
  H[2][2] = 2.3101089e-02;

  Vector<double> grad(3);
  grad[0]=1.2132635e+00;
  grad[1]=1.2132635e+00;
  grad[2]=-2.3976500e-01;
  ok = H.CholeskySolveLapack(chol,x,grad,false);
  std::cout << "choleksy_lapack_solver" << std::endl;
  std::cout << "H=" << H.ToString(TS_PYTHON) << std::endl;
  std::cout << "grad=" << grad.ToString(TS_PYTHON) << std::endl;
  std::cout << "x=" << x.ToString(TS_PYTHON);
  if (ok)
    std::cout << "-> ok" << std::endl;
  else
    std::cout << "-> not ok" << std::endl;
  std::cout << std::endl;

  BOOST_TEST(!ok);
}

BOOST_AUTO_TEST_CASE(Vector_ToString)
{
  Vector<double> vd(2, 3.0);
  BOOST_TEST(vd.ToString() == "[3, 3]");
  BOOST_TEST(vd.ToString(TS_PLAIN) == "3 3");
  BOOST_TEST(vd.ToString(TS_MATLAB) == "[3, 3]");
  BOOST_TEST(vd.ToString(TS_NONZEROS) == "0:3, 1:3");
  BOOST_TEST(vd.ToString(TS_PYTHON,", ",10) == "[3.0000000000e+00, 3.0000000000e+00]");

  Vector<int> vi(2, -3);
  BOOST_TEST(vd.ToString(TS_MATLAB, ", ",2) == "[3.00e+00, 3.00e+00]"); // digits are ignored in the int case

  Vector<std::complex<double> > vc(2, -3);
  BOOST_TEST(vc.ToString() == "[-3+0j, -3+0j]"); // Python style == e-tech style :)
  BOOST_TEST(vc.ToString(TS_MATLAB) == "[-3+0i, -3+0i]");

}


BOOST_AUTO_TEST_CASE(Matrix_ToString)
{
  Matrix<double> mr(2,2);
  mr.InitValue(2.0);
  mr[0][0] = 0.0;

  BOOST_TEST(mr.ToString(TS_PLAIN) == "0 2\n2 2");
  BOOST_TEST(mr.ToString(TS_MATLAB) == "[0 2\n2 2]");
  BOOST_TEST(mr.ToString(TS_MATLAB, "; ") == "[0 2; 2 2]");
  BOOST_TEST(mr.ToString(TS_PYTHON) == "[[0, 2],[2, 2]]");
  BOOST_TEST(mr.ToString(TS_INFO) == "rows=2 cols=2 nnz=3 min=0 max=2");
  BOOST_TEST(mr.ToString(TS_NONZEROS, ", ") == "row=0 1:2, row=1 0:2 1:2");
}

BOOST_AUTO_TEST_CASE(ConvertToVec_UpperTriangular)
{
  Matrix<int> m(6,6);
  for(int i = 0; i < 6; i++)
    for(int j= 0; j < 6; j++)
      m[i][j] = (i+1) * 10 + (j+1);

  //std::cout << "m(6,6) = " << m.ToString() << std::endl;
  Vector<int> v;
  m.ConvertToVec_UpperTriangular(v);
  std::cout << "v(m) = " << v.ToString() << std::endl;
}
