// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "polCoordSys.hh"
#include <cmath>
#include <sstream>

#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField{

  PolarCoordSystem::PolarCoordSystem( const std::string & name,
                                      Grid * ptGrid, 
                                      ParamNode * myParamNode ) 
    : CoordSystem( name, ptGrid, myParamNode ) {
    
    ENTER_FCN("PolarCoordSystem::PolarCoordSystem");
    
    // check, if grid has three dimensions
    if ( ptGrid_->GetDim() != 2 ) {
      EXCEPTION( "A polar coordinate system can only be used in "
                 << "a two-dimensional grid!" );
    }

    // get origin point of coordinate system
    GetPoint( origin_, myParam_->Get( "origin" ) );

    // get second point for defining the r-axis
    GetPoint( rAxis_, myParam_->Get( "rAxis" ) );

    // calculate the rotation matrix and the invers
    CalcRotationMatrix();

    // print information to .info-file
    PrintInfo();
    
  }
  
  PolarCoordSystem::~PolarCoordSystem(){
    ENTER_FCN("PolarCoordSystem::~PolarCoordSystem");
  }

  void PolarCoordSystem::Local2GlobalCoord( Vector<Double> & glob, 
                                            const Vector<Double> & loc ) const {
    ENTER_FCN("PolarCoordSystem::Local2GlobalCoord");
    
    Vector<Double> temp(2);

    if ( loc.GetSize() != 2 ) {
      EXCEPTION( "PolarCoordSystem::Local2GlobalCoord: Coordinate system only "
                 << "defined for 2-dimensional coordinates!" );
    }
    
    // transform cylindrical coords into local cartesian ones
    temp[0] = loc[0] * std::cos(loc[1]/180*PI);
    temp[1] = loc[0] * std::sin(loc[1]/180*PI);
    
     // rotate local cartesian coordinate system to global one
    glob.Resize(2);
    invRotationMat_.Mult(temp,glob);
    
    // add global coordinate midpoint
    glob += origin_;
    
  }
  
  void PolarCoordSystem::Global2LocalCoord( Vector<Double> & loc, 
                                          const Vector<Double> & glob ) const {
    ENTER_FCN("PolarCoordSystem:: Global2LocalCoord");

    Vector<Double> temp(2);
    if ( glob.GetSize() != 2 ) {
      EXCEPTION( "PolarCoordSystem::Local2GlobalCoord: Coordinate system only "
                 << "defined for 2-dimensional coordinates!" );
    }
    
    // calculate differential vector   
    Vector<Double> d(2);
    d = glob - origin_; 

    // rotate global cartesian coordinate system to local one
    rotationMat_.Mult(d,temp);

    // transform local cartesian nodes to polar ones
    loc.Resize(2);
    loc[0] = std::sqrt(temp[0] * temp[0] + temp[1] * temp[1]);
    loc[1] = std::atan2(temp[1],temp[0])/PI*180;
   
  }

  void  PolarCoordSystem::
  GetGlobRotationAngles( Vector<Double> & angles,
                         const Vector<Double>& point ) const {
     ENTER_FCN( "PolarCoordSystem::GetGlobRotationAngles" );

     Vector<Double> loc;
     Global2LocalCoord( loc, point );

     // Note: As this coordinate system can only be defined
     //       in a 2-dimensional plane, the rotation angle
     //       is always defined as rotation around the global
     //       z-axis, which is why we set it as the 3rd component.
     angles.Resize(3);
     angles.Init();
     angles[2] = loc[1] + (rotationAng_[0] * 180 / PI);
      
  }
  
  void PolarCoordSystem::
  Local2GlobalVector( Vector<Double> & globVec, 
                      const Vector<Double> & locVec, 
                      const Vector<Double> & globModelPoint ) const { 
    Local2GlobalVectorInt<Double>( globVec, locVec, globModelPoint );
  }

  void PolarCoordSystem::
  Local2GlobalVector( Vector<Complex> & globVec, 
                      const Vector<Complex> & locVec, 
                      const Vector<Double> & globModelPoint ) const { 
    Local2GlobalVectorInt<Complex>( globVec, locVec, globModelPoint );
  }


  template <class TYPE>
  void PolarCoordSystem::
  Local2GlobalVectorInt( Vector<TYPE> & globVec, 
                         const Vector<TYPE> & locVec, 
                         const Vector<Double> & globModelPoint ) const { 
    ENTER_FCN("PolarCoordSystem::Local2GlobalVectorInt");

    Double phi, r;
    Vector<Double> localPoint(2), d(2);
    Vector<TYPE> temp(2);

    // Calculate the distance and angle to the point
    d =  globModelPoint - origin_;

    // Transform global cartesian model point into local
    // cartesian one
    rotationMat_.Mult(d, localPoint);

    r = std::sqrt(localPoint[0] * localPoint[0] 
                  + localPoint[1] * localPoint[1]);
    phi = std::atan2(localPoint[1], localPoint[0]);

    // calculate global coordinate of locVec, by applying the
    // standard conversion routines
    globVec.Resize(2);
    temp.Resize(2);

    // Note: the second term in temp[0] and temp[1] accounts for the
    // tangential component of the locVector
    temp[0] = locVec[0] * std::cos(phi) - locVec[1]*std::sin(phi);
    temp[1] = locVec[0] * std::sin(phi) + locVec[1]*std::cos(phi);
    
    // Now we have the vector in cartesian coordinates for the
    // LOCAL cartesian system. To get the cartesian representation for
    // the GLOBAL one, we have to apply the inverse rotation matrix.
    globVec = invRotationMat_ * temp;

  }
    
  void PolarCoordSystem::CalcRotationMatrix() {
    ENTER_FCN( "CalcRotationMatrix" );
        
    Vector<Double> x(2), y(2);


    // 1) In order to calculate the rotation matrix, we have to find out
    //    how the local x' and y' axes are defined, because in a 'global'
    //    polar coordinate system we assume, that the y-axi is aligned
    //    with the axial axis and that the x-axis is aligned with the
    //    radial axis (where phi=0°).
    //    Therefore we can calculate the 'local' cartesian axes as follows:
    //    x': difference of 'rAxis-origin' and normalized
    //    y': rotation of x-axis by 90°

    // Vector<Double> temp;
//     temp = rAxis_- origin_;
//     if( temp.NormL2() < EPS ) {
//       EXCEPTION( "The pointing vector for the r_axis and the origin coincide "
//                   << "in the polar coordinate system '" << name_ << "'!" );
//     }

    x = rAxis_ - origin_;
    x /= x.NormL2();
    
    y[0] = -x[1];
    y[1] =  x[0];


    // Check if there were coincident points for defining the different axes
    if (x.NormL2() < EPS ||
        y.NormL2() < EPS ) {
      EXCEPTION( "At least two of your points for origin and  r-Axis "
                 << "coincide in the coordinate system '" << name_
                 << "'.\nPlease correct your parameter file!" );
    }


    // 2) Calculation of the rotation matrix, which defines the mapping
    //    from the global to the local coordinate system
    //    (ref. Bronstein: Taschenbuch der Mathematik, p. 217f)
    // Note: in order to avoid dividing by zero, an additional check
    //       is performed, if the x/y-component is in the order of
    //       machine precision.

    rotationMat_.Resize(2,2);
    rotationMat_[0][0] = (std::abs(x[0]) < EPS ) ? 0.0 : (x[0]);
    rotationMat_[0][1] = (std::abs(x[1]) < EPS ) ? 0.0 : (x[1]);

    rotationMat_[1][0] = (std::abs(y[0]) < EPS ) ? 0.0 : (y[0]);
    rotationMat_[1][1] = (std::abs(y[1]) < EPS ) ? 0.0 : (y[1]);

    // 3) Calculate transposed inverse rotation matrix, which defines 
    //    mapping from  local to global cartesian coordinate system
    rotationMat_.Invert(invRotationMat_);


    // Now calculate the related kardan angles for forward and 
    // backward transformation
    rotationAng_.Resize(1);
    rotationAng_[0] = GetAngle( x[1], x[0] );

    invRotationAng_.Resize(1);
    invRotationAng_[0] = -rotationAng_[0];
    
  }


  UInt PolarCoordSystem::GetVecComponent( const std::string & dof ) const{
    ENTER_FCN( "PolarCoordSystem::GetVecComponent" );
    
    
    UInt component = 0;
    
    if ( dof == "r" )
      component = 1;
    if ( dof == "phi" )
      component = 2;
    
    if ( component == 0 ) {
      EXCEPTION( "PolarSystem:GetVecComponent:\n"
                 << "The component with name '" << dof 
                 << "' is not known in the global cylinder coordinate system '"
                 << name_ << "'!" );
    }

    return component;
  }  

  const std::string PolarCoordSystem::GetDofName( const UInt dof ) const {
    ENTER_FCN( "PolarCoordSystem::GetDofName" );
    
    std::string ret = "";
    
    switch (dof) {
    case 1:
      ret = "r";
      break;
    case 2:
      ret = "phi";
      break;
    default:
      EXCEPTION( "PolarCoordSystem::GetDofName:\n"
                 << "The component number " << dof << " does not exist in a "
                 << "global cartesian coordinate system!" );
    }

    return ret;
  }

  void PolarCoordSystem::PrintInfo() {
    ENTER_FCN( "PolarCoordSystem::PrintInfo") ;

    std::ostringstream out;
    out << "\n--- local coordinate system ---\n"
        << "    name:\t" << name_ << std::endl
        << "    type:\tpolar" << std::endl
        << "  origin:\t" << origin_[0] << "," << origin_[1] << std::endl
        << "  r-axis:\t" << rAxis_[0] << "," << rAxis_[1] << std::endl;
    out << "   angle:\t" << rotationAng_[0]/PI*180  << "\n\n";
    Info->PrintF(std::string(), out.str().c_str());

  }


} // end of namespace
