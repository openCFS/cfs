// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "coordSystem.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"


namespace CoupledField{


  CoordSystem::CoordSystem(const std::string & name,
                           Grid * ptGrid,
                           ParamNode * myParamNode ) {

    name_ = name;
    ptGrid_ = ptGrid;
    myParam_ = myParamNode;

    if(ptGrid_)
      dim_ = ptGrid_->GetDim();
    else
      dim_ = 3;
  }

  CoordSystem::~CoordSystem(){
  }

  void CoordSystem::GetPoint(Vector<Double> & vec,
                             ParamNode * pointNode ) {

    
   
    std::string nodeName;
    StdVector<std::string> coordNames;
    StdVector<UInt> nodes;
    coordNames = "x", "y", "z";

    // check, if node is given by name and eventually get it
    // from the grid object
    pointNode->Get( "node", nodeName );

    vec.Resize(dim_);
    if ( nodeName != "none" ) {
      ptGrid_->GetNodesByName(nodes,nodeName);

      // check if more than one node is defined by this name
      if ( nodes.GetSize() != 1 ) {
        EXCEPTION( "CoordinateSystem: There are more than 1 nodes defined "
                   << "defined by '" << nodeName << "'.\nTherefore it is "
                   << "impossible to determine ONE coordinate location. Please "
                   << "ensure, that only ONE node is defined by this name!" );
      }

      ptGrid_->GetNodeCoordinate(vec, nodes[0]);

    } else {

      // if no node name was given, read in global x,y and z coordinate
      for (UInt i=0; i<vec.GetSize(); i++) {
        pointNode->Get( coordNames[i], vec[i] );
      }
    }
      
  }

  void CoordSystem::CalcKardanAngles( Vector<Double>& angles,
                                      Matrix<Double>& rotMat ) {

    angles.Resize(3);

    // Ref.: C. Woernle, "Skript: Dynamik von Mehrkörpersystemen,
    // Kapitel 2 "Grundlagen der Kinematik", S. 12, Univ. Rostock
    // http://iamserver.fms.uni-rostock.de/studium/mehrkoerpersysteme/unterlagen.htm
    

    // Safety check: T_33 must not be 1!
    if ( (std::abs(std::abs (rotMat[0][2]) - 1.0)) < EPS ) {
      EXCEPTION( "Rotation angle beta for coordinate system '"
                 << name_ << "' must not be 90°!" );
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
    
    // Calculate absolute value of angle ( 0 < alpha < pi/2)
    Double angle = std::abs(std::acos( cosAlpha ) );
    
    // Determine correct sign for angle
    if ( sinAlpha < 0 ) {
      angle *= -1.0;
    }

    return angle;
         
  }


} // end of namespace
