#include "algsys/olascomm.hh"
#include "utils/utils.hh"

//this file contains only the implementation of non-template
//functions. In the case of g++ 3.3 upwards we include the 
//template implementation as well
#if __GNUC_PREREQ(3,3) || defined (__INTEL_COMPILER)
#include "algsys/olascomm_impl.hh"
#endif

namespace OLAS {


  // *******************
  //   Deep Destructor
  // *******************
  OLAS_BaseComm::~OLAS_BaseComm() {
    ENTER_FCN( "OLAS_BaseComm::~OLAS_BaseComm" );
    stringPool_.clear();
    intPool_.clear();
    doublePool_.clear();

  }


  // ******************************
  //   Set string value for a key
  // ******************************
  void OLAS_BaseComm::SetValue( const std::string key, std::string value ) {
    ENTER_FCN( "OLAS_BaseComm::SetValue" );
    stringPool_[key] = value;
  }


  // *******************************
  //   Set Integer value for a key
  // *******************************
  void OLAS_BaseComm::SetValue( const std::string key, Integer value ) {
    ENTER_FCN( "OLAS_BaseComm::SetValue" );
    intPool_[key] = value;
  }


  // ******************************
  //   Set Double value for a key
  // ******************************
  void OLAS_BaseComm::SetValue( const std::string key, Double value ) {
    ENTER_FCN( "OLAS_BaseComm::SetValue" );
    doublePool_[key] = value;
  }


  // ******************************
  //   Set Bool value for a key
  // ******************************
  void OLAS_BaseComm::SetValue( const std::string key, bool value ) {
    ENTER_FCN( "OLAS_BaseComm::SetValue" );
    booleanPool_[key] = value;
  }  


  // *****************************************
  //   Query value of a key from string pool
  // *****************************************
  std::string OLAS_BaseComm::GetStringValue( std::string key ) {

    ENTER_FCN( "OLAS_BaseComm::GetStringValue" );

    // Look for key in pool
    stringMap::iterator position = stringPool_.find(key);
  
    // If key is not in pool, report an error
    if( position == stringPool_.end() ) {
      std::cerr << " Current pool contents:" << std::endl;
      ShowAll( std::cerr );
      std::string errmsg = "Cannot find key ";
      errmsg += key;
      errmsg += " in stringPool_";
      Error( errmsg.c_str(), __FILE__, __LINE__ );
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
     ENTER_FCN( "OLAS_BaseComm::HasValue" );
     
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

    ENTER_FCN( "OLAS_BaseComm::GetIntValue" );

    // look for key in pool
    intMap::iterator position = intPool_.find(key);
  
    // If key is not in pool, report an error
    if( position == intPool_.end() ) {
      std::cerr << " Current pool contents:" << std::endl;
      ShowAll( std::cerr );
      std::string errmsg = "Cannot find key ";
      errmsg += key;
      errmsg += " in intPool_";
      Error( errmsg.c_str(), __FILE__, __LINE__ );
    }

    // Else return key value
    return intPool_[key];
  }


  // *****************************************
  //   Query value of a key from double pool
  // *****************************************
  Double OLAS_BaseComm::GetDoubleValue( const std::string key ) {

    ENTER_FCN( "OLAS_BaseComm::GetDoubleValue" );

    // look for key in pool
    doubleMap::iterator position = doublePool_.find(key);
  
    // If key is not in pool, report an error
    if( position == doublePool_.end() ) {
      std::cerr << " Current pool contents:" << std::endl;
      ShowAll( std::cerr );
      std::string errmsg = "Cannot find key ";
      errmsg += key;
      errmsg += " in doublePool_";
      Error( errmsg.c_str(), __FILE__, __LINE__ );
    }

    // Else return key value
    return doublePool_[key];
  }


  // *****************************************
  //   Query value of a key from boolean pool
  // *****************************************
  bool OLAS_BaseComm::GetBoolValue( const std::string key ) {

    ENTER_FCN( "OLAS_BaseComm::GetBoolValue" );

    // look for key in pool
    booleanMap::iterator position = booleanPool_.find(key);
  
    // If key is not in pool, report an error
    if( position == booleanPool_.end() ) {
      std::cerr << " Current pool contents:" << std::endl;
      ShowAll( std::cerr );
      std::string errmsg = "Cannot find key ";
      errmsg += key;
      errmsg += " in booleanPool_";
      Error( errmsg.c_str(), __FILE__, __LINE__ );
    }

    // Else return key value
    return booleanPool_[key];
  }


  // *****************************
  //   Return the size of a pool
  // *****************************
  Integer OLAS_BaseComm::GetSize( const poolType pool ) const {
    ENTER_FCN( "OLAS_BaseComm::GetSize" );
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
      Error( "Wrong pool type", __FILE__, __LINE__ );
    }
    return 0;
  }


  // ****************
  //   Flush a pool
  // ****************
  void OLAS_BaseComm::FlushPool( const poolType pool ) {
    ENTER_FCN( "OLAS_BaseComm::FlushPool" );
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
      Error( "Wrong pool type", __FILE__, __LINE__ );
    }
  }


  // ***************************
  //   Show contents of a pool
  // ***************************
  void OLAS_BaseComm::ShowPool( const poolType pool,
				std::ostream &liststream ) {
    ENTER_FCN( "OLAS_BaseComm::ShowPool" );

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
      Error( "Wrong pool type", __FILE__, __LINE__ );
    }
  }


  // ******************************
  //   Show contents of all pools
  // ******************************
  void OLAS_BaseComm::ShowAll( std::ostream &liststream ) {
    ENTER_FCN( "OLAS_BaseComm::ShowAll" );
    ShowPool( INT_POOL    , liststream );
    ShowPool( DOUBLE_POOL , liststream );
    ShowPool( BOOLEAN_POOL, liststream );
    ShowPool( STRING_POOL , liststream );
    ShowPool( ENUM_POOL   , liststream );
  }


  // ****************************************
  //   Instantiate set/getEnumValue methods
  // ****************************************
  void OLAS_Params::EnumInterfaces() {

    // ReorderingType
    {
      SetValue( "Instantiation", NOREORDERING );
      ReorderingType val;
      GetEnumValue( "Instantiation", val );
    }
    // FEMatrixType
    {
      SetValue( "Instantiation", NOTYPE );
      FEMatrixType val;
      GetEnumValue( "Instantiation", val );
    }
    // PrecondType
    {
      SetValue( "Instantiation", NOPRECOND );
      PrecondType val;
      GetEnumValue( "Instantiation", val );
    }
    // SolverType
    {
      SetValue( "Instantiation", NOSOLVER );
      SolverType val;
      GetEnumValue( "Instantiation", val );
    }
    // EigenSolveType
    {
      SetValue( "Instantiation", NOEIGENSOLVER );
      EigenSolverType val;
      GetEnumValue( "Instantiation", val );
    }
    // CycleType
    {
      SetValue( "Instantiation", NOCYCLE );
      CycleType val;
      GetEnumValue( "Instantiation", val );
    }
    // MatrixStructureType
    {
      SetValue( "Instantiation", NOSTRUCTURETYPE );
      MatrixStructureType val;
      GetEnumValue( "Instantiation", val );
    }
    // MatrixEntryType
    {
      SetValue( "Instantiation", NOENTRYTYPE );
      MatrixEntryType val;
      GetEnumValue( "Instantiation", val );
    }
    // MatrixStorageType
    {
      SetValue( "Instantiation", NOSTORAGETYPE );
      MatrixStorageType val;
      GetEnumValue( "Instantiation", val );
    }
    // StopCrit
    {
      SetValue( "Instantiation", NOSTOPCRITTYPE );
      StopCritType val;
      GetEnumValue( "Instantiation", val );
    }
    // AMG interpolation
    {
      SetValue( "Instantiation", AMG_INTERPOLATION_CONSTANT );
      AMGInterpolationType val;
      GetEnumValue( "Instantiation", val );
    }
    // AMG smoother
    {
      SetValue( "Instantiation", AMG_SMOOTHER_GAUSSSEIDEL );
      AMGSmootherType val;
      GetEnumValue( "Instantiation", val );
    }
  }

  // Explicit template instantiation in the case of gnu g++,
  // as the compiler seems to be not satisfied with the old way the
  // templates are instantiated. This seems to happen only in g++
  // version < 4.0
#if __GNUC_PREREQ(4,0) || defined( __INTEL_COMPILER)
#define INSTANT_SETVAL( TYPE)                        \
  template void OLAS_BaseComm::                      \
  SetValue<TYPE>( const std::string key, TYPE value ) 
  
  INSTANT_SETVAL( ReorderingType );
  INSTANT_SETVAL( FEMatrixType );
  INSTANT_SETVAL( PrecondType );
  INSTANT_SETVAL( SolverType );
  INSTANT_SETVAL( EigenSolverType );
  INSTANT_SETVAL( CycleType );
  INSTANT_SETVAL( MatrixStructureType );
  INSTANT_SETVAL( MatrixEntryType );
  INSTANT_SETVAL( MatrixStorageType );
  INSTANT_SETVAL( StopCritType );
  INSTANT_SETVAL( AMGInterpolationType );
  INSTANT_SETVAL( AMGSmootherType );
#undef INSTANT_SETVAL
#endif

}
