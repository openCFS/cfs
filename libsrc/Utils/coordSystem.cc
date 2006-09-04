#include "coordSystem.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField{


  CoordSystem::CoordSystem(const std::string & name,
                           Grid * ptGrid) {
    ENTER_FCN("CoordSystem::CoordSystem");

    name_ = name;
    ptGrid_ = ptGrid;

    dim_ = ptGrid_->GetDim();
    
  }

  CoordSystem::~CoordSystem(){
    ENTER_FCN("CoordSystem::~CoordSystem");
  }

  void CoordSystem::GetPoint(Vector<Double> & vec,
                             StdVector<std::string> & keyVec,
                             StdVector<std::string> & attrVec,
                             StdVector<std::string> & valVec 
                             ) {
    
    ENTER_FCN( "CoordSystem::GetVector");
    
   
    std::string nodeName;
    StdVector<std::string> mKeyVec, coordNames;
    StdVector<UInt> nodes;
    coordNames = "x", "y", "z";

    // check, if node is given by name and eventually get it
    // from the grid object
    attrVec.Push_back(std::string());
    valVec.Push_back(std::string());    
    mKeyVec = keyVec;
    mKeyVec.Push_back("node");

    params->Get(mKeyVec, attrVec, valVec, nodeName);

    vec.Resize(dim_);
    if ( nodeName != "none" ) {
      ptGrid_->GetNodesByName(nodes,nodeName);

      // check if more than one node is defined by this name
      if ( nodes.GetSize() != 1 ) {
        (*error) << "CoordinateSystem: There are more than 1 nodes defined "
                 << "defined by '" << nodeName << "'.\nTherefore it is "
                 << "impossible to determine ONE coordinate location. Please "
                 << "ensure, that only ONE node is defined by this name!";
        Error( __FILE__, __LINE__ );
      }

      ptGrid_->GetNodeCoordinate(vec, nodes[0]);

    } else {

      // if no node name was given, read in global x,y and z coordinate
      for (UInt i=0; i<vec.GetSize(); i++) {
        mKeyVec = keyVec;
        mKeyVec.Push_back(coordNames[i]);    
        params->Get(mKeyVec, attrVec, valVec, vec[i]);
      }
    }
      
  }

  void CoordSystem::CalcKardanAngles( Vector<Double>& angles,
                                      Matrix<Double>& rotMat ) {
    ENTER_FCN( "CoordSystem::CalcKardanAngles" );

    angles.Resize(3);

    // Safety check: T_33 must not be 1!
    if ( (std::abs(std::abs (rotMat[0][2]) - 1.0)) < EPS ) {
      *error << "Rotation angle beta for coordinate system '"
             << name_ << "' must not be 90°!";
      Error( __FILE__, __LINE__ );
    }

    // Calculate  beta
    Double cos_beta = std::sqrt( 1 - rotMat[0][2]*rotMat[0][2] );
    Double sin_beta = rotMat[0][2];
    Double beta = GetAngle( sin_beta, cos_beta );
    
    // Calculate alpha
    Double cos_alpha = rotMat[2][2] / cos_beta;
    Double sin_alpha = -rotMat[1][2] / cos_beta;
    Double alpha = GetAngle( sin_alpha, cos_alpha );

    // Calculate gamma
    Double cos_gamma = rotMat[0][0] / cos_beta;
    Double sin_gamma = -rotMat[0][1] / cos_beta;
    Double gamma = GetAngle( sin_gamma, cos_gamma );

    // Fill angles into vector
    angles[0] = alpha;
    angles[1] = beta;
    angles[2] = gamma;
  }

  Double CoordSystem::GetAngle( Double sinAlpha, Double cosAlpha ) {
    ENTER_FCN( "CoordSystem::GetAngle ");
    
    // Calculate absolute value of angle ( 0 < alpha < pi/2)
    Double angle = std::abs(std::acos( cosAlpha ) );
    
    // Determine correct sign for angle
    if ( sinAlpha < 0 ) {
      angle *= -1.0;
    }

    return angle;
         
  }


} // end of namespace
