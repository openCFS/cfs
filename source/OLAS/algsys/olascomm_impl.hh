#include "algsys/olascomm.hh"
#include "utils/utils.hh"

// this file contains the implementation of 
// the template functions to allow including 
// it whereever it might be needed

namespace OLAS {


  // ************************************
  //   Set Enumeration value for a key
  // ************************************
  template <class enumT>
  void OLAS_BaseComm::SetValue( const std::string key, enumT value ) {
    ENTER_FCN( "OLAS_BaseComm::SetValue" );
    enumPool_[key] = (enumT) value;
  }


  // ***************************************
  //   Query value of a key from enum pool
  // ***************************************
  template <class enumT>
  void OLAS_BaseComm::GetEnumValue( const std::string key, enumT &value ) {
  
    ENTER_FCN( "OLAS_BaseComm::GetEnumValue" );
  
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
      std::string errmsg = "Cannot find key \"";
      errmsg += key;
      errmsg += "\" in enumPool_";
      Error( errmsg.c_str(), __FILE__, __LINE__ );
    }
  }

}
