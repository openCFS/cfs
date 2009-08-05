// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "resultInfo.hh"
#include <iostream>

#include "Domain/entityList.hh"
#include "General/environment.hh"

namespace CoupledField {

  ResultInfo::ResultInfo() {

    resultType = NO_SOLUTION_TYPE;
    dofNames.Clear();
    unit = "";
    complexFormat = AMPLITUDE_PHASE;
    entryType = UNKNOWN;
    definedOn = FREE;
  }


  UInt ResultInfo::GetDofIndex( const std::string & dof ) const {

    Integer pos = dofNames.Find( dof );
    if( pos < 0  ) {
      EXCEPTION( "Dof with name '" << dof << "' not found!" );
    }
    return (UInt) (pos-1);
    
  }


  std::string ResultInfo::GetDofName( const UInt dof ) const {
    if( dof <= 0 || dof > dofNames.GetSize()+1 ) {
      EXCEPTION( "'dof' must be in the range of [1.." 
                 << dofNames.GetSize()+1 << "]!" );
    }
    
    return dofNames[dof-1];
  }

  std::string ResultInfo::ToString() const 
  {
    std::ostringstream os;
    os << "name = " << resultName << " dofs " << dofNames.GetSize();
    return os.str();
  }



  ResultInfo& ResultInfo::operator=( const ResultInfo& data ) {

    resultType = data.resultType;
    resultName = data.resultName;
    dofNames = data.dofNames;
    unit = data.unit;
    complexFormat = data.complexFormat;
    entryType = data.entryType;
    definedOn = data.definedOn;
    fctType = data.fctType;
    
    EXCEPTION( "In operator ResultInfo::operator=\n" );
  }


  void ResultInfo::Enum2String(EntityUnknownType in, 
                               std::string& out ) {
    
    switch( in ) {
      
    case NODE: 
      out = "node";
      break;
    case EDGE:
      out = "edge";
      break;
    case FACE:  
      out = "face";
      break;
    case ELEMENT:  
      out = "element";
      break;
    case SURF_ELEM:
      out = "surfElement";
      break;
    case PFEM:
      out = "pfem";
      break;
    case REGION:
      out = "region";
      break;
    case SURF_REGION:  
      out = "surfRegion";
      break;
    case NODELIST:
      out = "nodeList";
      break;
    case COIL:
      out = "coil";
      break;
    case FREE:
      out = "free";
      break;
    default:
      EXCEPTION( "Conversion of " << in
                 << " to EntityUnkwownType not possible" );
    }
    
  }
  
  void  ResultInfo::String2Enum( const std::string& in,
                                  EntityUnknownType& out ) {

    if( in == "node")
      out =  NODE;
      else if( in == "edge")
        out = EDGE;
      else if( in == "face")
        out = FACE;
      else if( in == "element")
        out = ELEMENT; 
      else if( in == "surfElement")
        out = SURF_ELEM;
      else if( in == "pfem")
        out = PFEM;
      else if( in == "region")
        out = REGION;
      else if( in == "surfRegion")
        out = SURF_REGION;
      else if( in == "nodeList")
        out = NODELIST;
      else if( in == "coil")
        out = COIL;
      else if( in == "free")
        out = FREE;
      else {
        EXCEPTION( "Can not convert '" << in << "' to EntityUnknownType");
      }
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
