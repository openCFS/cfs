// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_STDVECTOR_SERIALIZE_HH
#define FILE_STDVECTOR_SERIALIZE_HH

#include "Utils/boost-serialization.hh"
#include "StdVector.hh"

// Define Serialization for StdVector<T> class
namespace boost {
  namespace serialization {
    // ************************************
    // Implement stuff for StdVector<T>
    // ************************************

    // Save method for vector class
    template<class T, class ARCHIVE> 
    void save( ARCHIVE & ar, const CoupledField::StdVector<T> &v, 
               const unsigned int version ) {
      UInt size = v.GetSize();
      ar << size;
      for( UInt i = 0; i < size; i++ ) 
        ar << v[i];
    }

    // Load method for vector class
    template<class T, class ARCHIVE>
    void load( ARCHIVE & ar, CoupledField::StdVector<T>& v, 
               const unsigned int version ) {

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
         CoupledField::StdVector<T> &t,
         const unsigned int file_version
     ){
         boost::serialization::split_free(ar, t, file_version);
     }

  } //end of namespace
 } //end of namespace
 
#endif // header guard
