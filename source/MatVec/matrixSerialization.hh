// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MATRIX_SERIALIZATION_HH
#define FILE_MATRIX_SERIALIZATION_HH

#include "matrix.hh"
#include "Utils/boost-serialization.hh"

// Define Serialization for Matrix<T> class
namespace boost {
  namespace serialization {
    
    // Save method for dense matrix class
    template<class T, class ARCHIVE> 
    void save( ARCHIVE & ar, const CoupledField::Matrix<T> &m, const unsigned int version ) {

      // save own members
      UInt sizeRow = m.GetNumRows();
      UInt sizeCol = m.GetNumCols();
      ar << sizeRow;
      ar << sizeCol;

      for( UInt i = 0; i < sizeRow; i++ ) {
        for( UInt j = 0; j < sizeCol; j++ ) {
          ar << m[i][j];
        }
      }
    }

    template<class T, class ARCHIVE>
    void load( ARCHIVE & ar, CoupledField::Matrix<T>& m, const unsigned int version ) {
      UInt sizeCol, sizeRow;
      ar >> sizeRow;
      ar >> sizeCol;

      // create storage for data to read in
      m.Resize( sizeRow, sizeCol );
      for( UInt i = 0; i < sizeRow; i++ ) {
        for( UInt j = 0; j < sizeCol; j++ ) { 
          ar >> m[i][j];
        }
      }

    }
    
    // split non-intrusive serialization function member into separate
    // non intrusive save/load member functions
    template<class T, class ARCHIVE>
    inline void serialize(
        ARCHIVE & ar,
        CoupledField::Matrix<T> &t,
        const unsigned int file_version
    ){
        boost::serialization::split_free(ar, t, file_version);
    }
  } //end of namespace
} //end of namespace

#endif  // header guard
