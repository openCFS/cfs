#include "DefaultCoordSystem.hh"
#include <cmath>

namespace CoupledField{

  DefaultCoordSystem::DefaultCoordSystem(Grid * ptGrid ) 
    : CoordSystem(std::string("default") , ptGrid, PtrParamNode() ) {
    
   // initialize rotation matrix
    rotationMat_.Resize(dim_, dim_);
    rotationMat_.Init();
    rotationMat_[0][0] = 1.0;
    rotationMat_[1][1] = 1.0;
    if( dim_ == 3 ) {
      rotationMat_[2][2] = 1.0;
    }
    invRotationMat_ = rotationMat_;
    
    // "calculate" full inverse rotation matrix
    invRotationMatFull_.Resize(3,3);
    invRotationMatFull_.Init();
    invRotationMatFull_.SetSubMatrix( invRotationMat_, 0, 0);
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

  void  DefaultCoordSystem::
  GetFullGlobRotationMatrix( Matrix<Double> & mat,
                             const Vector<Double>& point ) const {

    // no need to calculate anything
    mat = invRotationMatFull_;
  }

  UInt DefaultCoordSystem::GetVecComponent( const std::string& dof, bool silent )  const 
  {
    assert(ptGrid_ != nullptr);

    if(dof == "x" || (ptGrid_->IsAxi() && dof == "r"))
      return 1;
    if(dof == "y")
      return 2;
    if(dof == "z" && dim_ == 3 )
      return 3;
    if(dof == "z" && ptGrid_->IsAxi())
      return 2; 

    if(!silent)
      EXCEPTION( "DefaultCoordSystem:GetVecComponent: invalid dof '" << dof  << "' for dimension " << dim_);
    
    return INVALID_DOF;
  }

  
  const std::string DefaultCoordSystem::GetDofName( const UInt dof ) const {
    
    if( dof == 1 )
      return "x";

    if( dof == 2 )
      return "y";

    if( dof == 3 && dim_ == 3 )
      return "z";

    EXCEPTION( "DefaultCoordSystem::GetDofName:\n"
        << "The component number " << dof << " does not exist in a "
        << "global Cartesian coordinate system of dimension "
        << dim_ << "!" );
    
    return "";
  }


} // end of namespace
