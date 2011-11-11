#ifndef FILE_VECTOR_SERIALIZATION_HH
#define FILE_VECTOR_SERIALIZATION_HH

#include "SingleVector.hh"
#include "Vector.hh"
#include "Utils/boost-serialization.hh"

#include <boost/version.hpp>
#include <boost/serialization/split_member.hpp>
#if BOOST_VERSION < 103600
# include <boost/serialization/is_abstract.hpp>
#else
# include <boost/serialization/assume_abstract.hpp>
#endif


namespace boost {
  namespace serialization {
    

  // ************************************
  // Implement stuff for SingleVector
  // ************************************

    // Save method for vector class
    template<class ARCHIVE> 
    void serialize( ARCHIVE & ar, CoupledField::SingleVector &v, 
                    const unsigned int version ) {
    }

  }
}

//! Indicate that this is a pure virtual class
#if BOOST_VERSION < 103600
  BOOST_IS_ABSTRACT( CoupledField::SingleVector )
#else
  BOOST_SERIALIZATION_ASSUME_ABSTRACT( CoupledField::SingleVector )
#endif

namespace boost {
  namespace serialization {
    // ************************************
    // Implement stuff for Vector<T>
    // ************************************

    // Save method for vector class
    template<class T, class ARCHIVE> 
    void save( ARCHIVE & ar, const CoupledField::Vector<T> &v, 
               const unsigned int version ) {
      ar & boost::serialization::base_object<SingleVector>(v);
      UInt size = v.GetSize();
      ar << size;
      for( UInt i = 0; i < size; i++ ) 
        ar << v[i];
    }

    // Load method for vector class
    template<class T, class ARCHIVE>
    void load( ARCHIVE & ar, CoupledField::Vector<T>& v, 
               const unsigned int version ) {
      ar & boost::serialization::base_object<SingleVector>(v);   

      UInt size;
      ar >> size;

      v.Resize( size );
      for( UInt i = 0; i < size; i++ ) {
        ar >> v[i];
      }
    
    }
    // split non-intrusive serialization function member into separate
    // non intrusive save/load member functions
    template<class T, class ARCHIVE>
    inline void serialize(
        ARCHIVE & ar,
        CoupledField::Vector<T> &t,
        const unsigned int file_version
    ){
        boost::serialization::split_free(ar, t, file_version);
    }
    

  } //end of namespace
} //end of namespace

#endif // header guard
