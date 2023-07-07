#include <boost/test/unit_test.hpp>
#include <boost/tokenizer.hpp>
#include <string>
#include <iostream>
#include <algorithm>

#include "General/EnvironmentTypes.hh"

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

