#include <boost/test/unit_test.hpp>
#include <boost/tokenizer.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <boost/container/small_vector.hpp>

#include "General/EnvironmentTypes.hh"
#include "MatVec/Vector.hh"
#include "Utils/Timer.hh"

using namespace CoupledField;

/** the purpose of this file is to provide an easy test bed for algorithms and other code */

BOOST_AUTO_TEST_CASE(StringParse)
{
  std::stringstream ss( "1,.1 .2  .3\t.4 .5   ");
  // somehow the tokenizer does not work when using ss.str() directly
  std::string str = ss.str();

  boost::char_separator<char> sep(", \t;");
  boost::tokenizer<boost::char_separator<char> > tokens(str, sep);

  // for(const auto& t: tokens)
  //   std::cout  << "'"<< t << "'" << std::endl;

  // tokenizer has no size() :(

  auto iter = tokens.begin();
  BOOST_TEST((iter != tokens.end()));
  BOOST_TEST((*iter == "1"));
  BOOST_TEST(std::stoi(*iter) == 1);
  iter++;
  BOOST_TEST((iter != tokens.end()));
  BOOST_TEST((*iter == ".1"));
  BOOST_TEST(std::stod(*iter) == .1);
  iter++;
  BOOST_TEST((iter != tokens.end()));
  BOOST_TEST((*iter == ".2"));
  iter++;
  BOOST_TEST((iter != tokens.end()));
  BOOST_TEST((*iter == ".3"));
  iter++;
  BOOST_TEST((iter != tokens.end()));
  BOOST_TEST((*iter == ".4"));
  iter++;
  BOOST_TEST((iter != tokens.end()));
  BOOST_TEST((*iter == ".5"));
  iter++;
  BOOST_TEST((iter == tokens.end()));
}


// https://stackoverflow.com/questions/7728478/c-template-class-function-with-arbitrary-container-type-how-to-define-it
template <typename Iter>
void copy(Iter it, Iter end, int size = -1) {
    for(; it!=end; ++it)
      std::cout << *it << std::endl;
}


BOOST_AUTO_TEST_CASE(container_copy)
{
  std::list<int> l;
  l.push_back(1);
  l.push_back(2);
  l.push_back(3);

  copy(l.begin(), l.end());
}

double func_array() {  
  std::array<double, 20> a; 
  a.fill(0.0);
  return a[0];
}
double func_vector() {  
  Vector<double> v(20); // inits to 0.0
  return v[0];
}
double func_array_vector() { 
  std::array<double, 20> a;  
  Vector<double> v(a);
  // alternative Vector<double> v(a.size(), a.data(), false); // false means we're a wrapper
  v.Init();
  return v[0];
}
double func_small_vector() { 
  boost::container::small_vector<double, 20> s;  
  s.resize(20);
  std::fill(s.begin(), s.end(), 0.0); // much faster than assign()
  return s[0];
}
double func_small_vector_assign() { 
  boost::container::small_vector<double, 20> s;  
  s.assign(20, 0.0); // assign is very expensive?!
  return s[0];
}
double func_small_resize() { 
  boost::container::small_vector<double, 10> s;  
  s.resize(20);
  std::fill(s.begin(), s.end(), 0.0);
  return s[0];
}



BOOST_AUTO_TEST_CASE(heap_vs_stack)
{
  int N = 1e6;
  double t = 0.0;

  Timer array("", false, true);
  for(int i = 0; i < N; i++) t += func_array();
  array.Stop();

  Timer vector("", false, true);
  for(int i = 0; i < N; i++) t += func_vector();
  vector.Stop();

  Timer array_vector("", false, true);
  for(int i = 0; i < N; i++) t += func_array_vector();
  array_vector.Stop();

  // assign is very expensive
  Timer small_vector_assign("", false, true);
  for(int i = 0; i < N; i++) t += func_small_vector_assign();
  small_vector_assign.Stop();

  Timer small_vector("", false, true);
  for(int i = 0; i < N; i++) t += func_small_vector();
  small_vector.Stop();

  Timer small_resize("", false, true);
  for(int i = 0; i < N; i++) t += func_small_resize();
  small_resize.Stop();

  std::cout << "heap_vs_stack N=" << N << "s array=" << array.GetWallTime() << "s vector=" << vector.GetWallTime() << "s array_vector=" 
            << array_vector.GetWallTime() << "s small_vector_assign=" << small_vector_assign.GetWallTime() << "s small_vector=" << small_vector.GetWallTime() << "s small_resize=" << small_resize.GetWallTime() << "s" << " t=" << t << std::endl;
  // MacBook M1 Pro
  // heap_vs_stack N=1000000 no_fill=5.84e-07s array=4.58e-07s vector=0.0268434s array_vector=5.41e-07s small_vector_no_fill=3.75e-07s small_vector_assign=0.00480396s small_vector=3.75e-07s small_resize=0.0275428s t=nan

  // i9-14900KF
  // heap_vs_stack N=1000000 no_fill=0.00909007s array=0.000882648s vector=0.0106425s array_vector=0.000884599s small_vector_no_fill=0.000892204s small_vector_assign=0.000900859s small_vector=0.0008941s small_resize=0.0140576s t=6.95299e-304
         
}
