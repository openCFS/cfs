#include "resultDof.hh"


namespace CoupledField {

  ResultDof::ResultDof() {

    ENTER_FCN( "ResultDof::ResultDof");
  }


  UInt ResultDof::GetDofIndex( const std::string & dof ) const {
    ENTER_FCN( "ResultDof::GetDofIndex");

    Integer pos = dofNames.Find( dof );
    if( pos < 0  ) {
      (*error) << "Dof with name '" << dof << "' not found!";
      Error(  __FILE__, __LINE__ );
    }
    return (UInt) pos+1;
    
  }


  std::string ResultDof::GetDofName( const UInt dof ) const {
    ENTER_FCN( "ResultDof::GetDofName");
    if( dof <= 0 || dof > dofNames.GetSize()+1 ) {
      *error << "'dof' must be in the range of [1.." 
             << dofNames.GetSize()+1 << "]!";
      Error( __FILE__, __LINE__ );
    }
    
    return dofNames[dof-1];
  }


  ResultDof& ResultDof::operator=( const ResultDof& data ) {
    ENTER_FCN( "ResultDof::operator=");

    *warning << "In operator ResultDof::operator=\n";
    Warning( __FILE__, __LINE__ );

  }


  bool operator==( const ResultDof& a, const ResultDof& b ) {
    
    if (a.dofNames != b.dofNames) {
      std::cerr << "-> dofNames are different\n";
    }

    if ( a.resultType != b.resultType ||
         a.dofNames != b.dofNames ||
         a.definedOn != b.definedOn) {
 
      return false; 
    } else {
      return true;
    }
         
  }

  bool operator<( const ResultDof& a, const ResultDof& b ) {
    //return !(a==b);
    return ( a.resultType  < b.resultType );
  }
}
