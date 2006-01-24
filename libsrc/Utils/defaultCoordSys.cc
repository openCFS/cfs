#include "defaultCoordSys.hh"
#include <cmath>

namespace CoupledField{

  DefaultCoordSystem::DefaultCoordSystem(Grid * ptGrid) 
    : CoordSystem(std::string("default") , ptGrid) {
    
    ENTER_FCN("DefaultCoordSystem::DefaultCoordSystem");
   
  }
  
  DefaultCoordSystem::~DefaultCoordSystem(){
    ENTER_FCN("DefaultCoordSystem::~DefaultCoordSystem");
  }

  void DefaultCoordSystem::Local2GlobalCoord( Vector<Double> & glob, 
                                          const Vector<Double> & loc ) const {
    ENTER_FCN("DefaultCoordSystem::Local2GlobalCoord");
    glob = loc;
}
  
  void DefaultCoordSystem::Global2LocalCoord( Vector<Double> & loc, 
                                          const Vector<Double> & glob ) const {
    ENTER_FCN("DefaultCoordSystem:: Global2LocalCoord");
    loc = glob;
  }
  
  void DefaultCoordSystem::
  Local2GlobalVector( Vector<Double> & globVec, 
                      const Vector<Double> & locVec, 
                      const Vector<Double> & globModelPoint ) const { 
    ENTER_FCN("DefaultCoordSystem::Local2GlobalVector");

    globVec = locVec;
  }

  UInt DefaultCoordSystem::GetVecComponent( const std::string & dof )  const {
    ENTER_FCN( "DefaultCoordSystem::GetVecComponent" );

    
    UInt component = 0;
    
    //for mechanic PDE
    if ( dof == "ux" )
      component = 1;
    if ( dof == "uy" )
      component = 2;
    if ( dof == "uz" )
      component = 3;

    //for splitted acoustic PDE
    if ( dof == "px" )
      component = 1;
    if ( dof == "py" )
      component = 2;
    if ( dof == "pz" )
      component = 3;

    //for stokesFluid PDE
    if ( dof == "vx" )
      component = 1;
    if ( dof == "vy" )
      component = 2;
    if ( dof == "vz" )
      component = 3;
    if ( dof == "phi" )
      component = 4;
    if ( dof == "omegax" )
      component = 5;
    if ( dof == "omegay" )
      component = 6;
    if ( dof == "omegaz" )
      component = 7;
    if ( dof == "p" )
      component = 8;

    // --- DOF for piezo - HARD CODED -
    if ( dof == "ep" )
      component = ptGrid_->GetDim() + 1;
    
    if ( component == 0 ) {
      (*error) << "DefaultCoordSystem:GetVecComponent:\n"
               << "The component with name '" << dof 
               << "' is not known in the global cartesian coordinate system";
      Error( __FILE__, __LINE__ );
    }

    return component;
  }

  
  const std::string DefaultCoordSystem::GetDofName( const UInt dof ) const {
    ENTER_FCN( "DefaultCoordSystem::GetDofName" );
    
    std::string ret = "";
    
    switch (dof) {
    case 1:
      ret = "ux";
      break;
    case 2:
      ret = "uy";
      break;
    case 3:
      ret = "uz";
      break;
    default:
      (*error) << "DefaultCoordSystem::GetDofName:\n"
               << "The component number " << dof << " does not exist in a "
               << "global cartesian coordinate system!";
      Error( __FILE__, __LINE__ );
    }

    return ret;
  }

} // end of namespace
