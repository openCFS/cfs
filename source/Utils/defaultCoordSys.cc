// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "defaultCoordSys.hh"
#include <cmath>

namespace CoupledField{

  DefaultCoordSystem::DefaultCoordSystem(Grid * ptGrid ) 
    : CoordSystem(std::string("default") , ptGrid, PtrParamNode() ) {
    
   // initialize rotation matrix
    rotationMat_.Resize(3, 3);
    rotationMat_[0][0] = 1.0;
    rotationMat_[1][1] = 1.0;
    rotationMat_[2][2] = 1.0;
    invRotationMat_ = rotationMat_;
  }
  
  DefaultCoordSystem::~DefaultCoordSystem(){
  }

  void DefaultCoordSystem::Local2GlobalCoord( Vector<Double> & glob, 
                                          const Vector<Double> & loc ) const {
    glob = loc;
}
  
  void DefaultCoordSystem::Global2LocalCoord( Vector<Double> & loc, 
                                          const Vector<Double> & glob ) const {
    loc = glob;
  }
  
  
  void DefaultCoordSystem::
  Local2GlobalVector( Vector<Double> & globVec, 
                      const Vector<Double> & locVec, 
                      const Vector<Double> & globModelPoint ) const { 

    globVec = locVec;
  }

  void DefaultCoordSystem::
  Local2GlobalVector( Vector<Complex> & globVec, 
                      const Vector<Complex> & locVec, 
                      const Vector<Double> & globModelPoint ) const { 

    globVec = locVec;
  }
  
  void  DefaultCoordSystem::
  GetGlobRotationMatrix( Matrix<Double> & mat,
                         const Vector<Double>& point ) const {
   
    // no need to calculate anything
    mat = invRotationMat_;
  }


  UInt DefaultCoordSystem::GetVecComponent( const std::string & dof )  const {

    
    UInt component = 0;
  
    if ( dof == "x" )
      component = 1;
    if ( dof == "y" )
      component = 2;
    if ( dof == "z" )
      component = 3;

    if ( component == 0 ) {
      EXCEPTION( "DefaultCoordSystem:GetVecComponent:\n"
                 << "The component with name '" << dof 
                 << "' is not known in the global cartesian coordinate system" );
    }
    
    return component;
  }

  
  const std::string DefaultCoordSystem::GetDofName( const UInt dof ) const {
    
    std::string ret = "";
    
    switch (dof) {
    case 1:
      ret = "x";
      break;
    case 2:
      ret = "y";
      break;
    case 3:
      ret = "z";
      break;
    default:
      EXCEPTION( "DefaultCoordSystem::GetDofName:\n"
                 << "The component number " << dof << " does not exist in a "
                 << "global cartesian coordinate system!" );
    }

    return ret;
  }


} // end of namespace
