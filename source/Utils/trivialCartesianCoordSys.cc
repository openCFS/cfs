// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <cmath>
#include <sstream>

#include "DataInOut/ParamHandling/ParamNode.hh"

#include "trivialCartesianCoordSys.hh"


namespace CoupledField{

  TrivialCartesianCoordSystem::TrivialCartesianCoordSystem(const std::string & name,
                                 Grid * ptGrid,
                                 ParamNode * myParamNode ) 
    : CoordSystem( name, ptGrid, myParamNode ) {
    Vector<Double> originTemp, xAxisTemp, yAxisTemp;
    std::vector< std::string > axisMapTemp;
    axisMapTemp.push_back("x");
    axisMapTemp.push_back("y");
    axisMapTemp.push_back("z");

    // get origin point of coordinate system
    GetPoint( originTemp, myParam_->Get( "origin" ) );

    // get axis map
    ParamNode* node = myParam_->Get("axisMap", false);

    if(node)
    {
      node->Get("x", axisMapTemp[0], true);
      node->Get("y", axisMapTemp[1], true);
      node->Get("z", axisMapTemp[2], true);
    }

    origin_.Resize(3);
    axisFactors_.Resize(3);
    axisMap_.Resize(3);

    for(UInt i=0; i<originTemp.GetSize(); i++)
    {
      origin_[i] = originTemp[i];

      //      std::cout << "origin_[" << i << "]: " << origin_[i] << std::endl;
    }
    
    for(UInt i=0; i<3; i++)
    {
      const std::string& axis = axisMapTemp[i];
      
      if(axis == "x")
      {
        axisMap_[i] = 0;
        axisFactors_[i] = 1;
      } else if(axis == "-x")
      {
        axisMap_[i] = 0;
        axisFactors_[i] = -1;
      } else if(axis == "y")
      {
        axisMap_[i] = 1;
        axisFactors_[i] = 1;
      } else if(axis == "-y")
      {
        axisMap_[i] = 1;
        axisFactors_[i] = -1;
      } else if(axis == "z")
      {
        axisMap_[i] = 2;
        axisFactors_[i] = 1;
      } else if(axis == "-z")
      {
        axisMap_[i] = 2;
        axisFactors_[i] = -1;
      } else
      {
        EXCEPTION("Error while trying to build up trivial cartesian coord. sys.");
      }
    }

    // calculate the rotation matrix and the invers
    CalcRotationMatrix();

    // print information to .info-file
    PrintInfo();
    
  }
  
  TrivialCartesianCoordSystem::~TrivialCartesianCoordSystem(){
  }

  void TrivialCartesianCoordSystem::Local2GlobalCoord( Vector<Double> & glob, 
                                                       const Vector<Double> & loc ) const {
    Vector<Double> temp(3);
    UInt n=loc.GetSize();
    
    for(UInt i=0; i < n; i++)
      temp[i] = loc[i];
    
    // rotate local cartesian coordinate system to global one
    glob.Resize(3);
    invRotationMat_.Mult(temp,glob);

    // add global coordinate midpoint
    glob += origin_;
  }
  
  void TrivialCartesianCoordSystem::Global2LocalCoord( Vector<Double> & loc, 
                                                       const Vector<Double> & glob ) const {
    Vector<Double> temp(3);
    UInt n=glob.GetSize();
    loc.Resize(3);
    
    for(UInt i=0; i < n; i++)
      temp[i] = glob[i];
    
    // calculate differential vector   
    Vector<Double> d(3);
    d = temp - origin_; 
    
    // rotate global cartesian coordinate system to local one
    loc[0] = axisFactors_[0] * d[axisMap_[0]];
    loc[1] = axisFactors_[1] * d[axisMap_[1]];
    loc[2] = axisFactors_[2] * d[axisMap_[2]];
  }

  void  TrivialCartesianCoordSystem::
  GetGlobRotationAngles( Vector<Double> & angles,
                         const Vector<Double>& point ) const {
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
  
  void TrivialCartesianCoordSystem::
  Local2GlobalVector( Vector<Double> & globVec, 
                      const Vector<Double> & locVec, 
                      const Vector<Double> & globModelPoint ) const { 
    Local2GlobalVectorInt<Double>( globVec, locVec, globModelPoint );
  }

  void TrivialCartesianCoordSystem::
  Local2GlobalVector( Vector<Complex> & globVec, 
                      const Vector<Complex> & locVec, 
                      const Vector<Double> & globModelPoint ) const { 
    Local2GlobalVectorInt<Complex>( globVec, locVec, globModelPoint );
  }


  template <class TYPE>
  void TrivialCartesianCoordSystem::
  Local2GlobalVectorInt( Vector<TYPE> & globVec, 
                         const Vector<TYPE> & locVec, 
                         const Vector<Double> & globModelPoint ) const { 
    globVec.Resize(3);
    
    // We have the vector in cartesian coordinates for the
    // LOCAL cartesian system. To get the cartesian representation for
    // the GLOBAL one, we have to apply the inverse rotation matrix.
    globVec = invRotationMat_ * locVec;

  }
    
  void TrivialCartesianCoordSystem::CalcRotationMatrix() {
    Vector<Double> x(3), y(3), z(3), ytemp(3);
    Double fac;
    
    x[axisMap_[0]] = axisFactors_[0];

    //    std::cout << "x: " << x[0] << " " << x[1] << " " << x[2] << std::endl;

    y[axisMap_[1]] = axisFactors_[1];

    //    std::cout << "y: " << y[0] << " " << y[1] << " " << y[2] << std::endl;

    z[axisMap_[2]] = axisFactors_[2];

    //    std::cout << "z: " << z[0] << " " << z[1] << " " << z[2] << std::endl;

    // Check if there were coincident points for defining the different axes
    if (x.NormL2() < EPS ||
        y.NormL2() < EPS ||
        z.NormL2() < EPS ) {
      EXCEPTION( "At least two of the base vectors coincide"
                 << " in coordinate system '" << name_
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

    // 3) Calculate transposed inverse rotation matrix, which defines 
    //    mapping from  local to global cartesian coordinate system
    Matrix<Double> tempInvers;
    rotationMat_.Invert(invRotationMat_);

    //Now calculate the related kardan angles for forward and 
    //backward transformation
    CalcKardanAngles( rotationAng_, rotationMat_ );
    CalcKardanAngles( invRotationAng_, invRotationMat_ );

  }


  UInt TrivialCartesianCoordSystem::GetVecComponent( const std::string & dof ) const{
    UInt component = 0;
    
    if ( dof == "x" )
      component = 1;
    if ( dof == "y" )
      component = 2;
    if ( dof == "z" )
      component = 3;
    
    if ( component == 0 ) {
      EXCEPTION( "TrivialCartesianCoordSystem:GetVecComponent:\n"
                 << "The component with name '" << dof 
                 << "' is not known in the global cylinder coordinate system '"
                 << name_ << "'!" );
    }

    return component;
  }  

  const std::string TrivialCartesianCoordSystem::GetDofName( const UInt dof ) const {
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
      EXCEPTION( "TrivialCartesianCoordSystem::GetDofName:\n"
                 << "The component number " << dof << " does not exist in a "
                 << "global cartesian coordinate system!" );
    }

    return ret;
  }

  void TrivialCartesianCoordSystem::PrintInfo() {
    std::ostringstream out;
    out << "\n--- local coordinate system ---\n"
        << "    name:\t" << name_ << std::endl
        << "    type:\ttrivial cartesian" << std::endl
        << "  origin:\t" << origin_[0] << "," << origin_[1] << "," 
        << origin_[2] << std::endl
        << "  axisFactors:\t" << axisFactors_[0] << "," << axisFactors_[1] << ","
        << axisFactors_[2] << std::endl
        << "  axisMap:\t" << axisMap_[0] << "," << axisMap_[1] << ","
        << axisMap_[2] << "\n";
    out << "  angles:\t" << rotationAng_[0]/PI*180 << ","
        << rotationAng_[1]/PI*180 << "," << rotationAng_[2]/PI*180 << "\n\n";
    Info->PrintF(std::string(), out.str().c_str());

  }


} // end of namespace
