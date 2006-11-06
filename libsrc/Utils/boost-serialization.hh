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

// Additional includes
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/tracking.hpp>

#include "General/environment.hh"

// Define Conversion of Complex data
namespace boost {
  namespace serialization {
    template<class Archive>
    void load(Archive & ar, CoupledField::Complex & c, const unsigned int version) {
      Double real, imag;
      ar >> real;
      ar >>  imag;
      c = Complex( real, imag );
    }
    template<class Archive>
    void save(Archive & ar, const CoupledField::Complex & c, const unsigned int version) {
      const Double real = c.real();
      const Double imag = c.imag();
      ar << real;
      ar << imag;
    }
    
  }
}
BOOST_CLASS_TRACKING( CoupledField::Complex, boost::serialization::track_never );
BOOST_SERIALIZATION_SPLIT_FREE( CoupledField::Complex );


#endif
