// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "resultInfo.hh"


namespace CoupledField {

  ResultInfo::ResultInfo() {

    ENTER_FCN( "ResultInfo::ResultInfo");
    resultType = NO_SOLUTION_TYPE;
    dofNames = "";
    unit = "";
    complexFormat = AMPLITUDE_PHASE;
    entryType = UNKNOWN;
    definedOn = FREE;
    
    

  }


  UInt ResultInfo::GetDofIndex( const std::string & dof ) const {
    ENTER_FCN( "ResultInfo::GetDofIndex");

    Integer pos = dofNames.Find( dof );
    if( pos < 0  ) {
      (*error) << "Dof with name '" << dof << "' not found!";
      Error(  __FILE__, __LINE__ );
    }
    return (UInt) pos+1;
    
  }


  std::string ResultInfo::GetDofName( const UInt dof ) const {
    ENTER_FCN( "ResultInfo::GetDofName");
    if( dof <= 0 || dof > dofNames.GetSize()+1 ) {
      *error << "'dof' must be in the range of [1.." 
             << dofNames.GetSize()+1 << "]!";
      Error( __FILE__, __LINE__ );
    }
    
    return dofNames[dof-1];
  }


  ResultInfo& ResultInfo::operator=( const ResultInfo& data ) {
    ENTER_FCN( "ResultInfo::operator=");

    resultType = data.resultType;
    dofNames = data.dofNames;
    unit = data.unit;
    complexFormat = data.complexFormat;
    entryType = data.entryType;
    definedOn = data.definedOn;
    fctType = data.fctType;
    
    *warning << "In operator ResultInfo::operator=\n";
    Warning( __FILE__, __LINE__ );

  }


  bool operator==( const ResultInfo& a, const ResultInfo& b ) {
    
    if (a.dofNames != b.dofNames) {
      std::cerr << "-> dofNames are different\n";
    }

    if ( a.resultType != b.resultType ||
         a.dofNames != b.dofNames ||
         a.definedOn != b.definedOn ||
         a.entryType != b.entryType ) {
 
      return false; 
    } else {
      return true;
    }
         
  }

  bool operator<( const ResultInfo& a, const ResultInfo& b ) {
    //return !(a==b);
    return ( a.resultType  < b.resultType );
  }
}
