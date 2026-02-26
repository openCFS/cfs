#include <cmath>
#include <sstream>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "CartesianCoordSystem.hh"


namespace CoupledField{

  CartesianCoordSystem::CartesianCoordSystem(const std::string & name,
                                             Grid * ptGrid,
                                             PtrParamNode myParamNode ) 
  : CoordSystem( name, ptGrid, myParamNode )
  {
    Vector<Double> originTemp, xAxisTemp, yAxisTemp;

    // retrieve dimension
    
    
    // get origin point of coordinate system
    GetPoint( originTemp, myParam_->Get( "origin" ) );

    // get second point for defining the z-axis
    GetPoint( xAxisTemp, myParam_->Get( "xAxis" ) );

    
    // get second point for defining the r-axis
    if( dim_ == 3 ) {
      GetPoint( yAxisTemp, myParam_->Get( "yAxis" ) );
    }

    origin_.Resize(dim_);
    std::fill(&origin_[0], &origin_[0]+dim_, 0);
    xAxis_.Resize(dim_);
    std::fill(&xAxis_[0], &xAxis_[0]+dim_, 0);
    yAxis_.Resize(dim_);
    std::fill(&yAxis_[0], &yAxis_[0]+dim_, 0);
    
    for(UInt i=0; i<originTemp.GetSize(); i++)
    {
      origin_[i] = originTemp[i];
      xAxis_[i]  = xAxisTemp[i];
      if( dim_ == 3 ) {
        yAxis_[i]  = yAxisTemp[i];
      }
    }

    // calculate the rotation matrix and the invers
    CalcRotationMatrix();

  }
  
  CartesianCoordSystem::~CartesianCoordSystem(){
  }

  void CartesianCoordSystem::Local2GlobalCoord( Vector<Double> & glob, 
                                                const Vector<Double> & loc ) const {
    Vector<Double> temp(dim_);
    UInt n=loc.GetSize();
    if (dim_ < n) {
      n = dim_;
    }
    
    for(UInt i=0; i < n; i++)
      temp[i] = loc[i];
    
    // rotate local cartesian coordinate system to global one
    glob.Resize(dim_);
    invRotationMat_.Mult(temp,glob);

    // add global coordinate midpoint
    glob += origin_;
  }
  
  void CartesianCoordSystem::Global2LocalCoord( Vector<Double> & loc, 
                                                const Vector<Double> & glob ) const {
    Vector<Double> temp(dim_);
    UInt n=glob.GetSize();
    if (dim_ < n) {
      n = dim_;
    }
    loc.Resize(dim_);
    
    for(UInt i=0; i < n; i++)
      temp[i] = glob[i];
    
    // calculate differential vector   
    Vector<Double> d(dim_);
    d = temp - origin_; 
    
    // rotate global cartesian coordinate system to local one
    rotationMat_.Mult(d,loc);
  }
  
  void  CartesianCoordSystem::
  GetGlobRotationMatrix( Matrix<Double> & mat,
                         const Vector<Double>& point ) const {
    mat = invRotationMat_;
  }  
  
  void  CartesianCoordSystem::
  GetFullGlobRotationMatrix( Matrix<Double> & mat,
                             const Vector<Double>& point ) const {
    mat = invRotationMatFull_;
  }
  
  void CartesianCoordSystem::
  Local2GlobalVector( Vector<Double> & globVec, 
                      const Vector<Double> & locVec, 
                      const Vector<Double> & globModelPoint ) const { 
    Local2GlobalVectorInt<Double>( globVec, locVec, globModelPoint );
  }

  void CartesianCoordSystem::
  Local2GlobalVector( Vector<Complex> & globVec, 
                      const Vector<Complex> & locVec, 
                      const Vector<Double> & globModelPoint ) const { 
    Local2GlobalVectorInt<Complex>( globVec, locVec, globModelPoint );
  }


  template <class TYPE>
  void CartesianCoordSystem::
  Local2GlobalVectorInt( Vector<TYPE> & globVec, 
                         const Vector<TYPE> & locVec, 
                         const Vector<Double> & globModelPoint ) const { 
    globVec.Resize(3);
    
    // We have the vector in cartesian coordinates for the
    // LOCAL cartesian system. To get the cartesian representation for
    // the GLOBAL one, we have to apply the inverse rotation matrix.
    globVec = invRotationMat_ * locVec;

  }
    
  void CartesianCoordSystem::CalcRotationMatrix() {
    Vector<Double> x(dim_), y(dim_), z(dim_), ytemp(dim_);
    Double fac;
    

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
    x = xAxis_ - origin_;
    x /= x.NormL2();
    
    // === 2D Section ==
    if( dim_ == 2 ) {
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
    }
    
    // === 3D Section ==


    if( dim_ == 3 ) {
      ytemp = yAxis_ - origin_;

      fac = ytemp[0]*x[0] + ytemp[1]*x[1] + ytemp[2]*x[2];

      y[0] = ytemp[0] - (fac * x[0]);
      y[1] = ytemp[1] - (fac * x[1]);
      y[2] = ytemp[2] - (fac * x[2]);

      y /= y.NormL2();


      z[0] = x[1]*y[2] - x[2]*y[1];
      z[1] = x[2]*y[0] - x[0]*y[2];
      z[2] = x[0]*y[1] - x[1]*y[0];

      z /= z.NormL2();


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
    }
    
    // Perform on check on rotation matrix
    CheckRotationMat( rotationMat_ );
    
    // 3) Calculate transposed inverse rotation matrix, which defines 
    //    mapping from  local to global cartesian coordinate system
    Matrix<Double> tempInvers;
    rotationMat_.Invert(invRotationMat_);
    
    // "calculate" full inverse rotation matrix
    invRotationMatFull_.Resize(3,3);
    invRotationMatFull_.SetSubMatrix( invRotationMat_, 0, 0);
  }


  UInt CartesianCoordSystem::GetVecComponent( const std::string & dof, bool silent ) const
  {
    if ( dof == "x" )
      return 1;
    if ( dof == "y" )
      return 2;
    if ( dof == "z" && dim_ == 3 )
      return 3;

    if(!silent)
      EXCEPTION( "CartesianCoordSystem:GetVecComponent: invalid of '" << dof  << "' in a (rotated) Cartesian coordinate system for dimension " << dim_);
    
    return INVALID_DOF;
  }    

  const std::string CartesianCoordSystem::GetDofName( const UInt dof ) const {
    if( dof == 1 )
         return "x";

       if( dof == 2 )
         return "y";

       if( dof == 3 && dim_ == 3 )
         return "x";

       EXCEPTION( "CartesianCoordSystem::GetDofName:\n"
           << "The component number " << dof << " does not exist in a "
           << "(rotated) Cartesian coordinate system of dimension "
           << dim_ << "!" );
       
       return "";
  }

  void CartesianCoordSystem::ToInfo( PtrParamNode in ) {
    
    in = in->Get("cartesian");
    
    in->Get("id")->SetValue(name_);
    PtrParamNode originNode = in->Get("origin");
    originNode->Get("x")->SetValue(origin_[0]);
    originNode->Get("y")->SetValue(origin_[1]);
    if( dim_ == 3 )
      originNode->Get("z")->SetValue(origin_[2]);

    PtrParamNode xNode = in->Get("xAxis");
    xNode->Get("x")->SetValue(xAxis_[0]);
    xNode->Get("y")->SetValue(xAxis_[1]);
    if( dim_ == 3 )
      xNode->Get("z")->SetValue(xAxis_[2]);

    PtrParamNode yNode = in->Get("yAxis");
    yNode->Get("x")->SetValue(yAxis_[0]);
    yNode->Get("y")->SetValue(yAxis_[1]);
    if( dim_ == 3 )
      yNode->Get("z")->SetValue(yAxis_[2]);
    
//    PtrParamNode rotNode = in->Get("rotationAngles");
//    rotNode->Get("alpha")->SetValue(Rad2Grad(rotationAng_[0]));
//    rotNode->Get("beta")->SetValue(Rad2Grad(rotationAng_[1]));
//    rotNode->Get("gamma")->SetValue(Rad2Grad(rotationAng_[2]));
  }


} // end of namespace
