// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_OLASCOMM_HH
#define OLAS_OLASCOMM_HH

#include <map>
#include <string>
#include <iostream>
#include "utils/defs.hh"

namespace OLAS {

  //! Base class for the objects for passing information to and from OLAS
  class OLAS_BaseComm {

  public:

    //! The different pool types.

    //! This enumeration data type is used to distinguish the different
    //! pools. It offers the following values:
    //! - STRING_POOL
    //! - INT_POOL
    //! - DOUBLE_POOL
    //! - BOOLEAN_POOL
    //! - ENUM_POOL
    typedef enum { STRING_POOL, INT_POOL, DOUBLE_POOL, BOOLEAN_POOL,
		   ENUM_POOL } poolType;

    //! Default Destructor

    //! The default destructor needs to be deep.
    virtual ~OLAS_BaseComm();

    //! Set the value of a key in the string pool

    //! This method is used to change the value the key is connected to in the
    //! string pool. An existing value for this key ist over-written. If the
    //! key does not exist, it is added to the pool.
    void SetValue( const std::string key, std::string value );

    //! Set the value of a key in the Integer pool

    //! This method is used to change the value the key is connected to in the
    //! Integer pool. An existing value for this key ist over-written. If the
    //! key does not exist, it is added to the pool.
    void SetValue( const std::string key, Integer value );

    //! Set the value of a key in the Double pool

    //! This method is used to change the value the key is connected to in the
    //! Double pool. An existing value for this key ist over-written. If the
    //! key does not exist, it is added to the pool.
    void SetValue( const std::string key, Double value );

    //! Set the value of a key in the Boolean pool

    //! This method is used to change the value the key is connected to in the
    //! Boolean pool. An existing value for this key ist over-written. If the
    //! key does not exist, it is added to the pool.
    void SetValue( const std::string key, bool value );

    //! Set the value of a key in the Enum pool

    //! This method is used to change the value the key is connected to in the
    //! Enum pool. An existing value for this key ist over-written. If the
    //! key does not exist, it is added to the pool.
    template <class enumT>
    void SetValue( const std::string key, enumT value );

    //! Query the value of a key

    //! This method is used to query the value of a key connected to a string.
    //! If the key is not found in the pool the method will issue an error.
    std::string GetStringValue( std::string key );

    //! Query the value of a key

    //! This method is used to query the value of a key connected to an
    //! Integer.
    //! If the key is not found in the pool the method will issue an error.
    Integer GetIntValue( const std::string key );

    //! Query the value of a key

    //! This method is used to query the value of a key connected to a Double.
    //! If the key is not found in the pool the method will issue an error.
    Double GetDoubleValue( const std::string key );

    //! Query the value of a key

    //! This method is used to query the value of a key connected to a Bool.
    //! If the key is not found in the pool the method will issue an error.
    bool GetBoolValue( const std::string key );

    //! Query the value of a key

    //! This method is used to query the value of a key connected to an
    //! enumeration data type. If the key is not found in the pool the method
    //! will issue an error.
    template <class enumT>
    void GetEnumValue( const std::string key, enumT &value );

    //! Query the size of a pool

    //! The method returns the size of the pool specified by the pool input
    //! parameter.
    Integer GetSize( const poolType pool ) const;

    //! List all key/value pairs in a pool

    //! The method will print a list of all keys and their current values in
    //! the pool designated by the pool input parameter to the ouptut stream
    //! specified by liststream.
    void ShowPool( const poolType pool, std::ostream& liststream );

    //! List all key/value pairs in all pools

    //! The method will print a list of all keys and their current values in
    //! all pools to the ouptut stream specified by liststream.
    void ShowAll( std::ostream& liststream );

    //! Remove all key/value pairs from a pool
    void FlushPool( const poolType pool );

    /** This is a checker methof for the getter */
    bool HasKey(const std::string key);

    /** To check if we have a certain value on the map on the rhs */
    bool HasContent(const std::string& value); 

  protected:

    //! Shortcut for a type consisting of an STL map for strings

    //! A type that consists of an STL map which connects a string (key) to
    //! a string (value).
    typedef std::map<std::string, std::string> stringMap;

    //! Shortcut for a type consisting of an STL map for Integers

    //! A type that consists of an STL map which connects a string (key) to
    //! an Integer (value).
    typedef std::map<std::string, Integer> intMap;

    //! Shortcut for a type consisting of an STL map for Doubles

    //! A type that consists of an STL map which connects a string (key) to
    //! a Double (value).
    typedef std::map<std::string, Double> doubleMap;

    //! Shortcut for a type consisting of an STL map for Boolean

    //! A type that consists of an STL map which connects a string (key) to
    //! a Boolean (value).
    typedef std::map<std::string, bool> booleanMap;

    //! Shortcut for a type consisting of an STL map for Enums

    //! A type that consists of an STL map which connects a string (key) to
    //! an enum (value).
    typedef std::map<std::string, Integer> enumMap;
    
    //! Pool of string/string pairs
    stringMap stringPool_;

    //! Pool of string/Integer pairs
    intMap intPool_;

    //! Pool of string/Double pairs
    doubleMap doublePool_;

    //! Pool of string/Bool pairs
    booleanMap booleanPool_;

    //! Pool of string/Enum pairs
    enumMap enumPool_;

  };

  //! Class for objects that return solution information from OLAS

  //! This class implements a communication object that can be used to obtain
  //! information about the OLAS solution process.
  class OLAS_Report : public OLAS_BaseComm {
  };
  
}

//! if we are using icc or g++ 3.2 and older we include the
//! implementation of the template functions in the header
#ifndef DOXYGEN
#if(!__GNUC_PREREQ(3,3))
//#include "algsys/olascomm_impl.hh"
#endif
#endif
#endif
