// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef CFS_BOOST_SERIALIZATION
#define CFS_BOOST_SERIALIZATION

// Define various i/o-archives
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"

// Include support for standards STL-containers, which can be readily serialized
// if this the following 4 lines are commented, much more can be compiled with SUSE 11.0
#include "boost/serialization/vector.hpp"
#include "boost/serialization/map.hpp"
#include "boost/serialization/set.hpp"
#include "boost/serialization/list.hpp"

// Additional includes
#include "boost/serialization/split_free.hpp"
#include "boost/serialization/tracking.hpp"
#include "boost/serialization/nvp.hpp"

#include "General/defs.hh"

// Define Conversion of Complex data
namespace boost {
  namespace serialization {
    template<class Archive>
    void load(Archive & ar, CoupledField::Complex & c, const unsigned int version) {
      CoupledField::Double real, imag;
      ar >> real;
      ar >>  imag;
      c = CoupledField::Complex( real, imag );
    }
    template<class Archive>
    void save(Archive & ar, const CoupledField::Complex & c, const unsigned int version) {
      const CoupledField::Double real = c.real();
      const CoupledField::Double imag = c.imag();
      ar << real;
      ar << imag;
    }
    
  }
}


BOOST_CLASS_TRACKING( CoupledField::Complex, boost::serialization::track_never )
BOOST_SERIALIZATION_SPLIT_FREE( CoupledField::Complex )

#endif
