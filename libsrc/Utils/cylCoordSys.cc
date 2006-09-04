#include "cylCoordSys.hh"
#include <cmath>
#include <sstream>

namespace CoupledField{

  CylCoordSystem::CylCoordSystem(const std::string & name,
                                 Grid * ptGrid) 
    : CoordSystem(name, ptGrid) {
    
    ENTER_FCN("CylCoordSystem::CylCoordSystem");
    
    // check, if grid has three dimensions
    if ( ptGrid_->GetDim() != 3 ) {
      (*error) << "A cylindrical coordinate system can only be used in "
               << "a three dimensional grid!";
    }

    // get origin point of coordinate system
    StdVector<std::string>  keyVec, attrVec, valVec;
    keyVec = "domain", "coordSys", "cylindric", "origin";
    attrVec = "", "", "name";
    valVec = "", "", name;
    GetPoint(origin_, keyVec, attrVec, valVec);

    // get second point for defining the h-axis
    keyVec = "domain", "coordSys", "cylindric", "hAxis";
    attrVec = "", "", "name";
    valVec = "", "", name;
    GetPoint(hAxis_, keyVec, attrVec, valVec);

    // get second point for defining the r-axis
    keyVec = "domain", "coordSys", "cylindric", "rAxis";
    attrVec = "", "", "name";
    valVec = "", "", name;
    GetPoint(rAxis_, keyVec, attrVec, valVec);

    // calculate the rotation matrix and the invers
    CalcRotationMatrix();

    // print information to .info-file
    PrintInfo();
    
  }
  
  CylCoordSystem::~CylCoordSystem(){
    ENTER_FCN("CylCoordSystem::~CylCoordSystem");
  }

  void CylCoordSystem::Local2GlobalCoord( Vector<Double> & glob, 
                                          const Vector<Double> & loc ) const {
    ENTER_FCN("CylCoordSystem::Local2GlobalCoord");
    
    Vector<Double> temp(3);

    if ( loc.GetSize() != 3 ) {
      (*error) << "CylCoordSystem::Local2GlobalCoord: Coordinate system only "
               << "defined for 3-dimensional coordinates!";
      Error( __FILE__, __LINE__ );
    }
    
    // transform cylindrical coords into local cartesian ones
    temp[0] = loc[0] * std::cos(loc[1]);
    temp[1] = loc[0] * std::sin(loc[1]);
    temp[2] = loc[2];
    
    // rotate local cartesian coordinate system to global one
    glob.Resize(3);
    invRotationMat_.Mult(temp,glob);
    
  }
  
  void CylCoordSystem::Global2LocalCoord( Vector<Double> & loc, 
                                          const Vector<Double> & glob ) const {
    ENTER_FCN("CylCoordSystem:: Global2LocalCoord");

    Vector<Double> temp(3);
    if ( glob.GetSize() != 3 ) {
      (*error) << "CylCoordSystem::Local2GlobalCoord: Coordinate system only "
               << "defined for 3-dimensional coordinates!";
      Error( __FILE__, __LINE__ );
    }

    // rotate global cartesian coordinate system to local one
    rotationMat_.Mult(glob,temp);

    // transform local cartesian nodes to cylindrical ones
    loc.Resize(3);
    loc[0] = std::sqrt(temp[0] * temp[0] + temp[1] * temp[1]);
    loc[1] = std::atan2(temp[1],temp[0])/PI*180;
    loc[2] = temp[2];
                        
  }

  void  CylCoordSystem::
  GetGlobRotationAngles( Vector<Double> & angles,
                         const Vector<Double>& point ) const {
    ENTER_FCN( "CylCoordSystem::GetGlobRotationAngles" );
    Vector<Double> loc, anglesLoc, anglesGlob(3);
    Global2LocalCoord( loc, point );

    // Note: currently we simply return the 'phi' part, as this 
    // is the only changing one.
    anglesLoc.Resize(3);
    anglesLoc.Init();
    anglesLoc[2] = loc[1];

    // Now add to the point rotation angles the 
    // angles to rotate the angles back to the global system
    angles = invRotationAng_;
    angles *= 180 / PI;
    angles += anglesLoc;
      
  }

  
  void CylCoordSystem::
  Local2GlobalVector( Vector<Double> & globVec, 
                      const Vector<Double> & locVec, 
                      const Vector<Double> & globModelPoint ) const { 
    ENTER_FCN("CylCoordSystem::Local2GlobalVector");

    Double phi, r;
    Vector<Double> locModelPoint(3), temp(3), d(3);

    // Transform global cartesian model point into local
    // cartesian one
    rotationMat_.Mult(globModelPoint, locModelPoint);

    // Calculate the distance and angle to the point
    d =  locModelPoint;
    d -= origin_;

    r = std::sqrt(d[0] * d[0] + d[1] *d[1]);
    phi = std::atan2(d[1],d[0]);

    // calculate global coordinate of locVec, by applying the
    // standard conversion routines
    globVec.Resize(3);
    temp.Resize(3);

    // Note: the second term in temp[0] and temp[1] accounts for the
    // tangential component of the locVector
    temp[0] = locVec[0] * std::cos(phi) - locVec[1]*std::sin(phi);
    temp[1] = locVec[0] * std::sin(phi) + locVec[1]*std::cos(phi);
    temp[2] = locVec[2];
    
    // Now we have the vector in cartesian coordinates for the
    // LOCAL cartesian system. To get the cartesian representation for
    // the GLOBAL one, we have to apply the inverse rotation matrix.
    invRotationMat_.Mult(temp,globVec);

  }
    
  void CylCoordSystem::CalcRotationMatrix() {
    ENTER_FCN( "CalcRotationMatrix" );
        
    Vector<Double> x(3), y(3), z(3);
    Double r_dot;


    // 1) In order to calculate the rotation matrix, we have to find out,
    //    how the local x', y' and 'z' axes are defined, because in a 'global'
    //    cylindric coordinate system we assume, that the z-axis is aligned with
    //    the h-axis and that the axis for rho=0 is equal to the x-axis. Therefore
    //    we can compute the 'local' cartesian axes as follows
    //    z': defined by vector from origin to hAxis point
    //    x': is the normal part of the vector 'origin_-rAxis' w.r.t
    //        z'-Axis
    //    y': implicitly given by cross product of z' and  x' (right hand rule)
    z = hAxis_ - origin_;
    z /= z.NormL2();

    Vector<Double> temp;
    temp = rAxis_- origin_;
    r_dot = (temp * z) / hAxis_.NormL2();
    x = temp  - ((hAxis_ *r_dot) / hAxis_.NormL2());
    x /= x.NormL2();

    y[0] = z[1]*x[2] - z[2]*x[1];
    y[1] = z[2]*x[0] - z[0]*x[2];
    y[2] = z[0]*x[1] - z[1]*x[0];

    // Check if there were coincident points for defining the different axes
    if (x.NormL2() < EPS ||
        y.NormL2() < EPS ||
        z.NormL2() < EPS ) {
      (*error) << "At least two of your points for origin, h-Axis and "
               << "r-Axis coincide in the coordinate system '" << name_
               << "'.\nPlease correct your parameter file!";
      Error( __FILE__, __LINE__ );
    }


    // 2) Calculation of the rotation matrix, which defines the mapping
    //    from the global to the local coordinate system
    //    (ref. Bronstein: Taschenbuch der Mathematik, p. 217f)
    // Note: in order to avoid dividing by zero, an additional check
    //       is performed, if the x/y/z-component is in the order of
    //       machine precision.

    rotationMat_.Resize(3,3);
    rotationMat_[0][0] = (std::abs(x[0]) < EPS ) ? 0.0 : (x[0]);
    rotationMat_[0][1] = (std::abs(x[1]) < EPS ) ? 0.0 : (x[1]);
    rotationMat_[0][2] = (std::abs(x[2]) < EPS ) ? 0.0 : (x[2]);

    rotationMat_[1][0] = (std::abs(y[0]) < EPS ) ? 0.0 : (y[0]);
    rotationMat_[1][1] = (std::abs(y[1]) < EPS ) ? 0.0 : (y[1]);
    rotationMat_[1][2] = (std::abs(y[2]) < EPS ) ? 0.0 : (y[2]);

    rotationMat_[2][0] = (std::abs(z[0]) < EPS ) ? 0.0 : (z[0]);
    rotationMat_[2][1] = (std::abs(z[1]) < EPS ) ? 0.0 : (z[1]);
    rotationMat_[2][2] = (std::abs(z[2]) < EPS ) ? 0.0 : (z[2]);

    // 3) Calculate transposed inverse rotation matrix, which defines 
    //    mapping from  local to global cartesian coordinate system
    Matrix<Double> tempInvers;
    rotationMat_.Invert(invRotationMat_);

    // Now calculate the related kardan angles for forward and 
    // backward transformation
    CalcKardanAngles( rotationAng_, rotationMat_ );
    CalcKardanAngles( invRotationAng_, invRotationMat_ );
  }


  UInt CylCoordSystem::GetVecComponent( const std::string & dof ) const{
    ENTER_FCN( "CylCoordSystem::GetVecComponent" );
    
    
    UInt component = 0;
    
    if ( dof == "urad" )
      component = 1;
    if ( dof == "uphi" )
      component = 2;
    if ( dof == "uax" )
      component = 3;
    
    if ( component == 0 ) {
      (*error) << "CylSystem:GetVecComponent:\n"
               << "The component with name '" << dof 
               << "' is not known in the global cylinder coordinate system '"
               << name_ << "'!";
      Error( __FILE__, __LINE__ );
    }

    return component;
  }  

  const std::string CylCoordSystem::GetDofName( const UInt dof ) const {
    ENTER_FCN( "CylCoordSystem::GetDofName" );
    
    std::string ret = "";
    
    switch (dof) {
    case 1:
      ret = "urad";
      break;
    case 2:
      ret = "uphi";
      break;
    case 3:
      ret = "uax";
      break;
    default:
      (*error) << "CylCoordSystem::GetDofName:\n"
               << "The component number " << dof << " does not exist in a "
               << "global cartesian coordinate system!";
      Error( __FILE__, __LINE__ );
    }

    return ret;
  }

  void CylCoordSystem::PrintInfo() {
    ENTER_FCN( "CylCoordSystem::PrintInfo") ;

    std::ostringstream out;
    out << "\n--- local coordinate system ---\n"
        << "    name:\t" << name_ << std::endl
        << "    type:\tcylindrical" << std::endl
        << "  origin:\t" << origin_[0] << "," << origin_[1] << "," 
        << origin_[2] << std::endl
        << "  h-axis:\t" << hAxis_[0] << "," << hAxis_[1] << ","
        << hAxis_[2] << std::endl
        << "  r-axis:\t" << rAxis_[0] << "," << rAxis_[1] << ","
        << rAxis_[2] << "\n";
    out << "  angles:\t" << rotationAng_[0]/PI*180 << ","
        << rotationAng_[1]/PI*180 << "," << rotationAng_[2]/PI*180 << "\n\n";
    Info->PrintF(std::string(), out.str().c_str());

  }

} // end of namespace
