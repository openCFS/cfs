#include <cmath>
#include <sstream>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "TrivialCartesianCoordSystem.hh"


namespace CoupledField{

  TrivialCartesianCoordSystem::TrivialCartesianCoordSystem(const std::string & name,
                                                           Grid * ptGrid,
                                                           PtrParamNode myParamNode ) 
    : CoordSystem( name, ptGrid, myParamNode ) {
    Vector<Double> originTemp, xAxisTemp, yAxisTemp;
    std::vector< std::string > axisMapTemp;
    axisMapTemp.push_back("x");
    axisMapTemp.push_back("y");
    axisMapTemp.push_back("z");

    originTemp.Resize(3);
    originTemp.Init();
    
    if(myParam_) 
    {
      // get origin point of coordinate system
      GetPoint( originTemp, myParam_->Get( "origin" ) );
      
      // get axis map
      PtrParamNode node = myParam_->Get("axisMap", ParamNode::PASS );
      
      if(node)
      {
        node->GetValue("x", axisMapTemp[0], ParamNode::PASS);
        node->GetValue("y", axisMapTemp[1], ParamNode::PASS);
        node->GetValue("z", axisMapTemp[2], ParamNode::PASS);
      }
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

    // fast path: if origin is zero and axes are not remapped/flipped,
    // Local2GlobalCoord and Global2LocalCoord reduce to identity (no-op)
    isIdentity_ = (origin_[0] == 0.0 && origin_[1] == 0.0 && origin_[2] == 0.0 &&
                   axisMap_[0] == 0 && axisMap_[1] == 1 && axisMap_[2] == 2 &&
                   axisFactors_[0] == 1.0 && axisFactors_[1] == 1.0 && axisFactors_[2] == 1.0);
  }
  
  TrivialCartesianCoordSystem::~TrivialCartesianCoordSystem(){
  }

  void TrivialCartesianCoordSystem::Local2GlobalCoord( Vector<Double> & glob, 
                                                       const Vector<Double> & loc ) const {
    glob.Resize(3);
    if (isIdentity_)
    {
      for (auto i = 0; i < loc.GetSize(); i++) glob[i] = loc[i];
      for (auto i = loc.GetSize(); i < 3; i++) glob[i] = 0.0;
      return;
    }

    // rotate local cartesian coordinate system to global one
    if (loc.GetSize() == 3) 
    {
      invRotationMat_.Mult(loc, glob); // common 3D path: no temp allocation
    } 
    else 
    {
      std::array<double,3> tmp_array; // avoid dynamic allocation
      Vector<Double> temp(tmp_array.size(), tmp_array.data(), false); // don't copy data but wrap it
      for(unsigned int i = 0; i < loc.GetSize(); i++) 
        temp[i] = loc[i];
      invRotationMat_.Mult(temp, glob);
    }

    // add global coordinate midpoint
    glob += origin_;
  }
  
  void TrivialCartesianCoordSystem::Global2LocalCoord( Vector<Double> & loc, 
                                                       const Vector<Double> & glob ) const {
    loc.Resize(3);
    if (isIdentity_)
    {
      for (auto i = 0u; i < glob.GetSize(); i++) loc[i] = glob[i];
      for (auto i = glob.GetSize(); i < 3u; i++) loc[i] = 0.0;
      return;
    }
   
    // calculate differential vector (stack-allocated; zero-pads missing glob components)
    // we want to have d = glob - origin_, but this works only for 3D.
    // for 2D we want to have d[3] = -origin_[3], but glob[3] does not exist.
    assert(glob.GetSize() <= 3);
    std::array<double, 3> d{-origin_[0], -origin_[1], -origin_[2]};
    for(auto i = 0; i < glob.GetSize(); i++) 
      d[i] += glob[i];

    // rotate global cartesian coordinate system to local one
    loc[0] = axisFactors_[0] * d[axisMap_[0]];
    loc[1] = axisFactors_[1] * d[axisMap_[1]];
    loc[2] = axisFactors_[2] * d[axisMap_[2]];
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

  UInt TrivialCartesianCoordSystem::GetVecComponent( const std::string & dof ) const
  {
    if ( dof == "x" )
      return 1;
    if ( dof == "y" )
      return 2;
    if ( dof == "z" )
      return 3;
    
   EXCEPTION( "TrivialCartesianCoordSystem:GetVecComponent:\n"
                 << "The component with name '" << dof 
                 << "' is not known in the global cylinder coordinate system '"
                 << name_ << "'!" );
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

  void TrivialCartesianCoordSystem::ToInfo( PtrParamNode in ) {
    
    in = in->Get("trivialCartesianCoordinateSystem");
    
    in->Get("identity")->SetValue(isIdentity_);
    
    in->Get("id")->SetValue(name_);
    PtrParamNode originNode = in->Get("origin");
    originNode->Get("x")->SetValue(origin_[0]);
    originNode->Get("y")->SetValue(origin_[1]);
    originNode->Get("z")->SetValue(origin_[2]);
    
    PtrParamNode factorNode = in->Get("axisFactors");
    factorNode->Get("x")->SetValue(axisFactors_[0]);
    factorNode->Get("y")->SetValue(axisFactors_[1]);
    factorNode->Get("z")->SetValue(axisFactors_[2]);
    
    PtrParamNode mapNode = in->Get("axisMap");
    mapNode->Get("x")->SetValue(axisMap_[0]);
    mapNode->Get("y")->SetValue(axisMap_[1]);
    mapNode->Get("z")->SetValue(axisMap_[2]);
  }


} // end of namespace
