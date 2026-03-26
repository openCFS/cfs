#include "CylCoordSystem.hh"
#include <cmath>
#include <sstream>

#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField{

  CylCoordSystem::CylCoordSystem(const std::string & name,
                                 Grid * ptGrid,
                                 PtrParamNode myParamNode ) 
    : CoordSystem( name, ptGrid, myParamNode ) {
    
    
    // check, if grid has three dimensions
    if ( ptGrid_->GetDim() != 3 ) {
      EXCEPTION( "A cylindrical coordinate system can only be used in "
                 << "a three dimensional grid!" );
    }

    // get origin point of coordinate system
    GetPoint( origin_, myParam_->Get( "origin" ) );

    // get second point for defining the z-axis
    GetPoint( zAxis_, myParam_->Get( "zAxis" ) );

    // get second point for defining the r-axis
    GetPoint( rAxis_, myParam_->Get( "rAxis" ) );

    // calculate the rotation matrix and the invers
    CalcRotationMatrix();

  }
  
  CylCoordSystem::~CylCoordSystem(){
  }

  void CylCoordSystem::Local2GlobalCoord( Vector<Double> & glob, 
                                          const Vector<Double> & loc ) const {
    
    Vector<Double> temp(3);

    if ( loc.GetSize() != 3 ) {
      EXCEPTION( "CylCoordSystem::Local2GlobalCoord: Coordinate system only "
                 << "defined for 3-dimensional coordinates!" );
    }
    
    // transform cylindrical coords into local cartesian ones
    temp[0] = loc[0] * std::cos(loc[1]/180*M_PI);
    temp[1] = loc[0] * std::sin(loc[1]/180*M_PI);
    temp[2] = loc[2];

    
    // rotate local cartesian coordinate system to global one
    glob.Resize(3);
    invRotationMat_.Mult(temp,glob);


    // add global coordinate midpoint
    glob += origin_;

  }
  
  void CylCoordSystem::Global2LocalCoord( Vector<Double> & loc, 
                                          const Vector<Double> & glob ) const {

    Vector<Double> temp(3);
    if ( glob.GetSize() != 3 ) {
      EXCEPTION( "CylCoordSystem::Local2GlobalCoord: Coordinate system only "
                 << "defined for 3-dimensional coordinates!" );
    }

   // calculate differential vector   
   Vector<Double> d(3);
   d = glob - origin_; 

   // rotate global cartesian coordinate system to local one
   rotationMat_.Mult(d,temp);

   // transform local cartesian nodes to cylindrical ones
   loc.Resize(3);
   loc[0] = std::sqrt(temp[0] * temp[0] + temp[1] * temp[1]);
   loc[1] = std::atan2(temp[1],temp[0])/M_PI*180;
   loc[2] = temp[2];
   
  }

  void  CylCoordSystem::
  GetGlobRotationMatrix( Matrix<Double> & mat,
                         const Vector<Double>& point ) const {


    Vector<Double> loc;
    Global2LocalCoord( loc, point );

    // calculate directional sine / cosine
    Double s = sin( loc[1] / 180 * M_PI );
    Double c = cos( loc[1] / 180 * M_PI );

    mat.Resize(3,3);
    const Matrix<Double> & a = invRotationMat_;

    // perform explicit multiplication of a * r
    // where r is the rotation matrix around the z-axis
    //     ( c  -s  0 )
    // r = ( s   c  0 )
    //     ( 0   0  1 )

    mat[0][0] =  c * a[0][0] + s * a[0][1];
    mat[0][1] = -s * a[0][0] + c * a[0][1];
    mat[0][2] = a[0][2];

    mat[1][0] =  c * a[1][0] + s * a[1][1];
    mat[1][1] = -s * a[1][0] + c * a[1][1];
    mat[1][2] = a[1][2];

    mat[2][0] =  c * a[2][0] + s * a[2][1];
    mat[2][1] = -s * a[2][0] + c * a[2][1];
    mat[2][2] = a[2][2];
  }
  
  void  CylCoordSystem::
  GetFullGlobRotationMatrix( Matrix<Double> & mat,
                             const Vector<Double>& point ) const {
    // Here we can simply call the normal GetGlobRoationMatrix() method,
    // as a cylindric coordinate system can be only defined in 3D
    GetGlobRotationMatrix(mat, point);
  }
  
  void CylCoordSystem::
  Local2GlobalVector( Vector<Double> & globVec, 
                      const Vector<Double> & locVec, 
                      const Vector<Double> & globModelPoint ) const { 
    Local2GlobalVectorInt<Double>( globVec, locVec, globModelPoint );
  }

  void CylCoordSystem::
  Local2GlobalVector( Vector<Complex> & globVec, 
                      const Vector<Complex> & locVec, 
                      const Vector<Double> & globModelPoint ) const { 
    Local2GlobalVectorInt<Complex>( globVec, locVec, globModelPoint );
  }


  template <class TYPE>
  void CylCoordSystem::
  Local2GlobalVectorInt( Vector<TYPE> & globVec, 
                         const Vector<TYPE> & locVec, 
                         const Vector<Double> & globModelPoint ) const { 

    Double phi;
    // Double r;
    Vector<Double> localPoint(3), d(3);
    Vector<TYPE> temp(3);

    // Calculate the distance and angle to the point
    d =  globModelPoint - origin_;
    
    // Transform global cartesian model point into local
    // cartesian one
    rotationMat_.Mult(d, localPoint);


    // r = std::sqrt(localPoint[0] * localPoint[0] 
    //               + localPoint[1] * localPoint[1]);
    phi = std::atan2(localPoint[1], localPoint[0]);
    
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
    globVec = invRotationMat_ * temp;

  }
    
  void CylCoordSystem::CalcRotationMatrix() {
        
    Vector<Double> x(3), y(3), z(3);

    // 1) In order to calculate the rotation matrix, we have to find out,
    //    how the local x', y' and 'z' axes are defined, because in a 'global'
    //    cylindric coordinate system we assume, that the z-axis is aligned with
    //    the h-axis and that the axis for rho=0 is equal to the x-axis. Therefore
    //    we can compute the 'local' cartesian axes as follows
    //    z': defined by vector from origin to hAxis point
    //    x': is the normal part of the vector 'origin_-rAxis' w.r.t
    //        z'-Axis
    //    y': implicitly given by cross product of z' and  x' (right hand rule)
    //z = hAxis_ - origin_;
    z = zAxis_;
    z /= z.NormL2();



    //Vector<Double> temp;
    //temp = rAxis_- origin_;
    // if( temp.NormL2() < EPS ) {
//       EXCEPTION( "The pointing vector for the r_axis and the origin coincide "
//                  << "in the cylindrical coordinate system '"
//                  << name_ << "'!" );
//     }
    
    x = rAxis_;
    x /= x.NormL2();


    y[0] = z[1]*x[2] - z[2]*x[1];
    y[1] = z[2]*x[0] - z[0]*x[2];
    y[2] = z[0]*x[1] - z[1]*x[0];


    // Check if there were coincident points for defining the different axes
    if (x.NormL2() < EPS ||
        y.NormL2() < EPS ||
        z.NormL2() < EPS ) {
      EXCEPTION( "At least two of the vectors for origin, z-Axis and "
                 << "r-Axis coincide in the coordinate system '" << name_
                 << "'.\nPlease correct your parameter file!" );
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

    // Perform on check on rotation matrix
    CheckRotationMat( rotationMat_ );

    
    // 3) Calculate transposed inverse rotation matrix, which defines 
    //    mapping from  local to global cartesian coordinate system
    Matrix<Double> tempInvers;
    rotationMat_.Invert(invRotationMat_);

    
//    CalcKardanAngles( rotationAng_, rotationMat_ );
//    CalcKardanAngles( invRotationAng_, invRotationMat_ );
  }


  UInt CylCoordSystem::GetVecComponent( const std::string & dof, bool silent ) const
  {
    if ( dof == "r" )
      return 1;
    if ( dof == "phi" )
      return 2;
    if ( dof == "z" )
      return 3;

    if(!silent)
      EXCEPTION( "CylSystem:GetVecComponent: invalid dof '" << dof  << "' in local cylindrical coordinate system " << name_);
    
    return INVALID_DOF;
  }    
    
  const std::string CylCoordSystem::GetDofName( const UInt dof ) const {
    
    std::string ret = "";
    
    switch (dof) {
    case 1:
      ret = "r";
      break;
    case 2:
      ret = "phi";
      break;
    case 3:
      ret = "z";
      break;
    default:
      EXCEPTION( "CylCoordSystem::GetDofName:\n"
                 << "The component number " << dof << " does not exist in a "
                 << "local cylindric coordinate system!" );
    }

    return ret;
  }

  void CylCoordSystem::ToInfo( PtrParamNode in ) {

    in = in->Get("cylindrical");
    
    in->Get("id")->SetValue(name_);
    PtrParamNode originNode = in->Get("origin");
    originNode->Get("x")->SetValue(origin_[0]);
    originNode->Get("y")->SetValue(origin_[1]);
    originNode->Get("z")->SetValue(origin_[2]);

    PtrParamNode rNode = in->Get("rAxis");
    rNode->Get("x")->SetValue(rAxis_[0]);
    rNode->Get("y")->SetValue(rAxis_[1]);
    rNode->Get("z")->SetValue(rAxis_[2]);
    
    PtrParamNode zNode = in->Get("zAxis");
    zNode->Get("x")->SetValue(zAxis_[0]);
    zNode->Get("y")->SetValue(zAxis_[1]);
    zNode->Get("z")->SetValue(zAxis_[2]);

    
//    PtrParamNode rotNode = in->Get("rotationAngles");
//    rotNode->Get("alpha")->SetValue(Rad2Grad(rotationAng_[0]));
//    rotNode->Get("beta")->SetValue(Rad2Grad(rotationAng_[1]));
//    rotNode->Get("gamma")->SetValue(Rad2Grad(rotationAng_[2]));

  }


} // end of namespace
