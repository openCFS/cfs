// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/basematrix.hh"

#include "olascomm.hh"

namespace CoupledField {


  // *******************
  //   Deep Destructor
  // *******************
  OLAS_BaseComm::~OLAS_BaseComm() {
    stringPool_.clear();
    intPool_.clear();
    doublePool_.clear();

  }

  // ******************************
  //   Set string value for a key
  // ******************************
  void OLAS_BaseComm::SetValue( const std::string key, std::string value ) {
    stringPool_[key] = value;
  }


  // *******************************
  //   Set Integer value for a key
  // *******************************
  void OLAS_BaseComm::SetValue( const std::string key, Integer value ) {
    intPool_[key] = value;
  }


  // ******************************
  //   Set Double value for a key
  // ******************************
  void OLAS_BaseComm::SetValue( const std::string key, Double value ) {
    doublePool_[key] = value;
  }


  // ******************************
  //   Set Bool value for a key
  // ******************************
  void OLAS_BaseComm::SetValue( const std::string key, bool value ) {
    booleanPool_[key] = value;
  }  


  // *****************************************
  //   Query value of a key from string pool
  // *****************************************
  std::string OLAS_BaseComm::GetStringValue( std::string key ) {


    // Look for key in pool
    stringMap::iterator position = stringPool_.find(key);
  
    // If key is not in pool, report an error
    if( position == stringPool_.end() ) {
      std::cerr << " Current pool contents:" << std::endl;
      ShowAll( std::cerr );
      std::string errmsg = "Cannot find key ";
      errmsg += key;
      errmsg += " in stringPool_";
      EXCEPTION( errmsg.c_str() );
    }

    // Else return key value
    return stringPool_[key];
  }

  bool OLAS_BaseComm::HasContent(const std::string& value)
  {
    
     for(stringMap::iterator iter = stringPool_.begin(); iter != stringPool_.end(); iter++)
     {
        if(value == iter->second) return true;
     }

     return false;
  }


  bool OLAS_BaseComm::HasKey(const std::string key)
  {
     
     // check all pools!
     if(intPool_.find(key)    != intPool_.end()) return true;
     if(stringPool_.find(key) != stringPool_.end()) return true;
     if(doublePool_.find(key) != doublePool_.end()) return true;
     if(booleanPool_.find(key)!= booleanPool_.end()) return true;
     if(enumPool_.find(key)   != enumPool_.end()) return true;
     
     return false; // in no pools
  }

  // ******************************************
  //   Query value of a key from Integer pool
  // ******************************************
  Integer OLAS_BaseComm::GetIntValue( const std::string key ) {


    // look for key in pool
    intMap::iterator position = intPool_.find(key);
  
    // If key is not in pool, report an error
    if( position == intPool_.end() ) {
      std::cerr << " Current pool contents:" << std::endl;
      ShowAll( std::cerr );
      EXCEPTION( "Cannot find key " << key << " in intPool_");
    }

    // Else return key value
    return intPool_[key];
  }


  // *****************************************
  //   Query value of a key from double pool
  // *****************************************
  Double OLAS_BaseComm::GetDoubleValue( const std::string key ) {


    // look for key in pool
    doubleMap::iterator position = doublePool_.find(key);
  
    // If key is not in pool, report an error
    if( position == doublePool_.end() ) {
      std::cerr << " Current pool contents:" << std::endl;
      ShowAll( std::cerr );
      EXCEPTION( "Cannot find key " << key << " in doublePool_");
    }

    // Else return key value
    return doublePool_[key];
  }


  // *****************************************
  //   Query value of a key from boolean pool
  // *****************************************
  bool OLAS_BaseComm::GetBoolValue( const std::string key ) {


    // look for key in pool
    booleanMap::iterator position = booleanPool_.find(key);
  
    // If key is not in pool, report an error
    if( position == booleanPool_.end() ) {
      std::cerr << " Current pool contents:" << std::endl;
      ShowAll( std::cerr );
      EXCEPTION( "Cannot find key " << key << " in booleanPool_");
    }

    // Else return key value
    return booleanPool_[key];
  }


  // *****************************
  //   Return the size of a pool
  // *****************************
  Integer OLAS_BaseComm::GetSize( const poolType pool ) const {
    switch( pool ) {
    case STRING_POOL:
      return (Integer) stringPool_.size();
      break;
    case INT_POOL:
      return (Integer) intPool_.size();
      break;
    case DOUBLE_POOL:
      return (Integer) doublePool_.size();
      break; 
    case BOOLEAN_POOL:
      return (Integer) booleanPool_.size();
      break;
    case ENUM_POOL:
      return (Integer) enumPool_.size();
      break;
    default:
      EXCEPTION( "Wrong pool type" );
    }
    return 0;
  }


  // ****************
  //   Flush a pool
  // ****************
  void OLAS_BaseComm::FlushPool( const poolType pool ) {
    switch( pool ) {
    case STRING_POOL:
      stringPool_.clear();
      break;
    case INT_POOL:
      intPool_.clear();
      break;
    case DOUBLE_POOL:
      doublePool_.clear();
      break;
    case BOOLEAN_POOL:
      booleanPool_.clear();
      break;
    case ENUM_POOL:
      enumPool_.clear();
      break;
    default:
      EXCEPTION( "Wrong pool type" );
    }
  }


  // ***************************
  //   Show contents of a pool
  // ***************************
  void OLAS_BaseComm::ShowPool( const poolType pool,
				std::ostream &liststream ) {

    switch( pool ) {

    case STRING_POOL:
      for ( stringMap::iterator it = stringPool_.begin();
	    it != stringPool_.end(); it++ ) {
	liststream << " stringPool_:\t(key) " << it->first << "\t(value) "
		   << it->second << std::endl;
      }
      break;

    case INT_POOL:
      for ( intMap::iterator it = intPool_.begin();
	    it != intPool_.end(); it++ ) {
	liststream << " intPool_:\t(key) " << it->first << "\t(value) "
		   << it->second << std::endl;
      }
      break;

    case DOUBLE_POOL:
      for ( doubleMap::iterator it = doublePool_.begin();
	    it != doublePool_.end(); it++ ) {
	liststream << " doublePool_:\t(key) " << it->first << "\t(value) "
		   << it->second << std::endl;
      }
      break;

    case BOOLEAN_POOL:
      for ( booleanMap::iterator it = booleanPool_.begin();
	    it != booleanPool_.end(); it++ ) {
	liststream << " booleanPool_:\t(key) " << it->first << "\t(value) "
		   << it->second << std::endl;
      }
      break;
      
    case ENUM_POOL:
      for ( enumMap::iterator it = enumPool_.begin();
	    it != enumPool_.end(); it++ ) {
	liststream << " enumPool_:\t(key) " << it->first << "\t(value) "
		   << it->second << std::endl;
      }
      break;

    default:
      EXCEPTION( "Wrong pool type" );
    }
  }


  // ******************************
  //   Show contents of all pools
  // ******************************
  void OLAS_BaseComm::ShowAll( std::ostream &liststream ) {
    ShowPool( INT_POOL    , liststream );
    ShowPool( DOUBLE_POOL , liststream );
    ShowPool( BOOLEAN_POOL, liststream );
    ShowPool( STRING_POOL , liststream );
    ShowPool( ENUM_POOL   , liststream );
  }

}
