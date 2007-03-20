#include "defaultCoordSys.hh"
#include <cmath>

namespace CoupledField{

  DefaultCoordSystem::DefaultCoordSystem(Grid * ptGrid) 
    : CoordSystem(std::string("default") , ptGrid, NULL ) {
    
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

  void DefaultCoordSystem::
  Local2GlobalVector( Vector<Complex> & globVec, 
                      const Vector<Complex> & locVec, 
                      const Vector<Double> & globModelPoint ) const { 
    ENTER_FCN("DefaultCoordSystem::Local2GlobalVector");

    globVec = locVec;
  }

  UInt DefaultCoordSystem::GetVecComponent( const std::string & dof )  const {
    ENTER_FCN( "DefaultCoordSystem::GetVecComponent" );

    
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
    ENTER_FCN( "DefaultCoordSystem::GetDofName" );
    
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
