#ifndef CFS_BOOST_SERIALIZATION
#define CFS_BOOST_SERIALIZATION

// Define various i/o-archives
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

// Include support for standards STL-containers, which cna be readily serialized
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/list.hpp>

#include "General/environment.hh"

// Define Conversion of Complex data
namespace boost {
  namespace serialization {
    template<class Archive>
    void serialize(Archive & ar, CoupledField::Complex & c, const unsigned int version) {
      ar & c.real();
      ar & c.imag();
    }
  }
}

#endif
