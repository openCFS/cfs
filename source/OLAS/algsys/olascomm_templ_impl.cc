// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/basematrix.hh"
#include "General/environment.hh"

#include "olascomm.hh"

namespace CoupledField {


  // ************************************
  //   Set Enumeration value for a key
  // ************************************
  template <class enumT>
  void OLAS_BaseComm::SetValue( const std::string key, enumT value ) {
    enumPool_[key] = (enumT) value;
  }


  // ***************************************
  //   Query value of a key from enum pool
  // ***************************************
  template <class enumT>
  void OLAS_BaseComm::GetEnumValue( const std::string key, enumT &value ) {
  
  
    // look for key in pool
    enumMap::iterator position = enumPool_.find(key);
  
    // If key is in pool, return value
    if( position != enumPool_.end() ) {
      value = (enumT) enumPool_[key];
    }
  
    // Otherwise, report an error
    else {
      std::cerr << " Current pool contents:" << std::endl;
      ShowPool( INT_POOL, std::cerr );
      ShowPool( DOUBLE_POOL, std::cerr ); 
      ShowPool( BOOLEAN_POOL, std::cerr );
      ShowPool( STRING_POOL, std::cerr );
      ShowPool( ENUM_POOL, std::cerr ); 
      EXCEPTION( "Cannot find key \"" << key << "\" in enumPool_");
    }
  }

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  
#define INSTANT_SETVAL(TYPE)                           \
  template void OLAS_BaseComm::                       \
  SetValue<TYPE>( const std::string key, TYPE value ); \
  template void OLAS_BaseComm::                       \
  GetEnumValue<TYPE>( const std::string key, TYPE &value )

  
  INSTANT_SETVAL( ReorderingType );
  INSTANT_SETVAL( FEMatrixType );
  INSTANT_SETVAL( PrecondType );
  INSTANT_SETVAL( SolverType );
  INSTANT_SETVAL( EigenSolverType );
  INSTANT_SETVAL( CycleType );
  INSTANT_SETVAL( BaseMatrix::StructureType );
  INSTANT_SETVAL( BaseMatrix::EntryType );
  INSTANT_SETVAL( BaseMatrix::StorageType );
  INSTANT_SETVAL( StopCritType );
  INSTANT_SETVAL( AMGInterpolationType );
  INSTANT_SETVAL( AMGSmootherType );
#undef INSTANT_SETVAL

#endif
  
}
