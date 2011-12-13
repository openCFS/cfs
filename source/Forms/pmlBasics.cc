// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <cmath>
#include <complex>

#include "Domain/domain.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/exprt/xpr1.hh"
#include "Utils/coordSystem.hh"
#include "Utils/tools.hh"
#include "pmlBasics.hh"

namespace CoupledField
{

  PMLBasics::PMLBasics( std::string dampingTypePML, Double damp, std::string type)
  {
    formsType_ = type;

    //check correct damping type
    if ( dampingTypePML == "constant" ) {
      dampingTypePML_ = dampingTypePML;
    }
    else if ( dampingTypePML == "quadDist" ) {
      dampingTypePML_ = dampingTypePML;
    }
    else if ( dampingTypePML == "inverseDist" ) {
      dampingTypePML_ = dampingTypePML;
    }
    else if ( dampingTypePML == "smoothDamp" ) {
      dampingTypePML_ = dampingTypePML;
    }    else {
      EXCEPTION("Damping type for PML not known");
    }

    dampingFactor_ = damp;
  }


  PMLBasics::~PMLBasics()
  {
  }


  void PMLBasics::ComputeFactorPML(Vector<Complex>& factorsPML, Complex& pmlDet, 
                                   Vector<Double> & globPos, Double omega)
  {

    UInt numVals = globPos.GetSize();
    Vector<Complex> factors(numVals);
    Double omegaInv = 1.0 / omega;
    Double imagVal;

    // new approach: perform addition of damping factors
    Vector<Double> tempFactors(numVals), realFactors(numVals);
    realFactors.Init();

    // get local coordinate representation
    Vector<Double> locPos (globPos.GetSize());
    coordSys_->Global2LocalCoord( locPos, globPos );
    
    // Calculate damping factors sigma_x, sigma_y, sigma_z,
    // which are always related to the local coordinate system
    for( UInt i = 0; i < numVals; ++i ) {
      if ( locPos[i] < minLoc_[i] || locPos[i] > maxLoc_[i] ) {
        ComputeDampingFactor(tempFactors, globPos, locPos, i);
        realFactors += tempFactors;
      }      
    }
    
    // set entries for nu_x, nu_y, nu_z and calculate
    // complex valued determinant
    pmlDet = Complex(1.0,0);
    for( UInt i = 0; i < numVals; ++i ) {
      imagVal = std::abs(realFactors[i]) * omegaInv;
      factors[i] = Complex(1.0, -imagVal);
      pmlDet *= factors[i];
    }

    // compute complex valued factors
    factorsPML.Resize(numVals);
    if ( formsType_ == "laplaceInt" || 
         formsType_ == "mixedPMLMassPPInt" || 
         formsType_ == "mixedPMLMassVVInt" ) {
      for( UInt i = 0; i < numVals; ++i ) {
        factorsPML[i] = Complex(1.0,0) / factors[i]; 
      } // for
    } // forms
  }

  //factors for almost PML
  void PMLBasics::ComputeFactorAPML(Vector<Complex>& factorsPML,  
                                    Vector<Double> & globPos, 
                                    Double omega)
  {

    UInt numVals = globPos.GetSize();

    // new approach: perform addition of damping factors
    Vector<Double> tempFactors(numVals), realFactors(numVals);
    realFactors.Init();

    // get local coordinate representation
    Vector<Double> locPos (globPos.GetSize());
    coordSys_->Global2LocalCoord( locPos, globPos );

    // Calculate damping factors sigma_x, sigma_y, sigma_z,
    // which are always related to the local coordinate system
    for( UInt i = 0; i < numVals; ++i ) {
      if ( locPos[i] < minLoc_[i] || locPos[i] > maxLoc_[i] ) {
        ComputeDampingFactor(tempFactors, globPos, locPos, i);
        realFactors += tempFactors;
      }      
    }

    factorsPML.Resize(numVals);
    factorsPML.Init();

    if ( formsType_ == "laplaceIntAPML") {
      //x-part
      factorsPML[0] = Complex(realFactors[1]+realFactors[2],omega) 
        / Complex( realFactors[0], omega);
//         + Complex(realFactors[1]*realFactors[2],0) 
//         / Complex(-omega*omega, omega* realFactors[0]); 

      //y-part
      factorsPML[1] = Complex(realFactors[0]+realFactors[2],omega) 
        / Complex( realFactors[1], omega);
//         +  Complex(realFactors[0]*realFactors[2],0) 
//         / Complex(-omega*omega, omega* realFactors[1]);
 
      //z-part
      factorsPML[2] = Complex(realFactors[0]+realFactors[1],omega) 
        / Complex( realFactors[2], omega);
//         + Complex(realFactors[0]*realFactors[1],0) 
//         / Complex(-omega*omega, omega* realFactors[2]); 
    }
    else {
      //PML factor for mass integrator
      factorsPML[0] = Complex(1, (-1.0/omega)*(realFactors[0]+realFactors[1]+realFactors[2]) )
        + Complex( -1/(omega*omega)* (  realFactors[0]*realFactors[1] 
                                        + realFactors[0]*realFactors[2]
                                        + realFactors[1]*realFactors[2] ), 0 );
      //      + Complex(0,1/(omega*omega*omega) * realFactors[0]*realFactors[1]*realFactors[2]);
    }
  }

  void PMLBasics::ComputeTimeFactorPML(Vector<Double>& factorsPML, 
                                       Vector<Double>& globPos )
  {

    UInt numVals = globPos.GetSize();
    
    Vector<Double> tempFactors(numVals);
    factorsPML.Resize(numVals);
    factorsPML.Init();

    // get local coordinate representation
    Vector<Double> locPos (globPos.GetSize());
    coordSys_->Global2LocalCoord( locPos, globPos );


    // Calculate damping factors sigma_x, sigma_y, sigma_z,
    // which are always related to the local coordinate system
    for( UInt i = 0; i < numVals; ++i ) {
      if ( locPos[i] < minLoc_[i] || locPos[i] > maxLoc_[i] ) {
        ComputeDampingFactor(tempFactors, globPos, locPos, i);
        factorsPML += tempFactors;
      }      
    }
  }


  void PMLBasics::ComputeTimeFactorPML(Vector<Double>& factorsPML, 
                                       Vector<Double>& pos,
                                       Vector<Double> center)
  {
    EXCEPTION("This method is not adjusted yet to arbitrary local "
              "coordinate systems. Please contact M. Kaltenbacher "
              "or A. Hauck ");
//
//    UInt numVals = pos.GetSize();
//
//    factorsPML.Resize(numVals);
//    factorsPML.Init();
//
//    if ( (pos[0] < minX_ || pos[0] > maxX_) || (center[0] < minX_ || center[0] > maxX_) ) {
//      factorsPML[0] = ComputeDampingFactor(pos, X);
//
//      if ( (pos[1] < minY_ || pos[1] > maxY_) || (center[1] < minY_ || center[1] > maxY_)) {
//	      factorsPML[1] = ComputeDampingFactor(pos, Y);
//
//	      //check for 3D
//	      if (numVals == 3) {
//	        if  ((pos[2] < minZ_ || pos[2] > maxZ_) || (center[2] < minZ_ || center[2] > maxZ_) ) {
//	          //compute z-value
//	          factorsPML[2] = ComputeDampingFactor(pos, Z);
//	        }
//	      }
//      }
//      
//      else {
//	      //check for 3D
//	      if (numVals == 3) {
//	        if  ((pos[2] < minZ_ || pos[2] > maxZ_) || (center[2] < minZ_ || center[2] > maxZ_) ) {
//	          //compute z-value
//	          factorsPML[2]  = ComputeDampingFactor(pos, Z);
//	        }
//	      }
//      }
//    } else {
//      if ( (pos[1] < minY_ || pos[1] > maxY_) || (center[1] < minY_ || center[1] > maxY_)) {
//	      //compute y-value
//	      factorsPML[1]  = ComputeDampingFactor(pos, Y);
//
//	      //check for 3D
//	      if (numVals == 3) {
//	        if  ((pos[2] < minZ_ || pos[2] > maxZ_) || (center[2] < minZ_ || center[2] > maxZ_) ) {
//	          //compute z-value
//	          factorsPML[2]  = ComputeDampingFactor(pos, Z);
//	        }
//	      }
//      } else {
//	      //check for 3D
//	      if (numVals == 3) {
//	        if  ((pos[2] < minZ_ || pos[2] > maxZ_) || (center[2] < minZ_ || center[2] > maxZ_) ) {
//	          factorsPML[2]  = ComputeDampingFactor(pos, Z);
//	        }
//	      }
//      }
//    }
    return;
  }

  void PMLBasics::ComputeDampingFactor( Vector<Double>& factors, 
                                        const Vector<Double>& globPos, 
                                        const Vector<Double>& locPos,
                                        UInt dir )
  {

    Double factor, maxPos, delta, diffCoord;
    
    // vector with local direction
    Vector<Double> locDir(locPos.GetSize());
    locDir.Init();
    locDir[dir] = 1.0;
    
    factor = 1.0;
    if ( dampingTypePML_ == "constant" ) {
      factor = dampingFactor_;
    }

    else if ( dampingTypePML_ == "quadDist" ) {
	      
      //get correct layer thickness
      if ( locPos[dir] < minLoc_[dir] ) {
        delta = layerThickness_[0][dir];
        diffCoord = abs(locPos[dir] - minLoc_[dir]);
      }
      else {
        delta = layerThickness_[1][dir];
        diffCoord = abs(locPos[dir] - maxLoc_[dir]);
      }
      factor = dampingFactor_ * ( diffCoord*diffCoord )/ ( delta*delta );
    } 
    
    else if ( dampingTypePML_ == "inverseDist" ) {

      //get correct maximal PML local coordinate
      if ( locPos[dir] < minLoc_[dir] ) {
        maxPos = minLoc_[dir] - layerThickness_[0][dir];
      }
      else {
        maxPos = maxLoc_[dir] + layerThickness_[1][dir];
      }

      //! This is just for Lobatto integration
      Double  divisor = abs(maxPos - locPos[dir]);
      if ( abs (maxPos - locPos[dir]) < 1e-12 ) {
        factor = 0;
      }else{
        factor = abs (dampingFactor_ / divisor );
      }
    }
    
    else if ( dampingTypePML_ == "smoothDamp" ) {
      
      //get correct layer thickness
      if ( locPos[dir] < minLoc_[dir] ) {
        delta = layerThickness_[0][dir];
        diffCoord = abs(locPos[dir]) - minLoc_[dir];
      } else {
        delta = layerThickness_[1][dir];
        diffCoord = abs(locPos[dir]) - maxLoc_[dir];
      }
      factor = dampingFactor_ * 
        ( diffCoord/delta - std::sin(2*PI*diffCoord/delta) / ( 2*PI ) );
    }
    
    // factors vector is just scalar factor in "local" direction
    factors = locDir * factor;
  }


  void PMLBasics::SetPosPML( Matrix<Double> & inner, 
                             Matrix<Double> & outer,
                             const std::string& coordSysId)
  {

    // inner/outer:   xmin  ymin  zmin
    //                xmax  ymax  zmax

    // store maximum / minimum coordinates (in local coordinates)
    // of pml layer (we assume axis-parallel boundaries)
    UInt dim = outer.GetNumCols();
    minLoc_.Resize(dim);
    maxLoc_.Resize(dim);
    for( UInt i = 0; i < dim; ++i ) {
      minLoc_[i] = inner[0][i];
      maxLoc_[i] = inner[1][i];
    }

    // compute thickness of pml layer
    layerThickness_.Resize(2,inner.GetNumCols());
    for (UInt i=0; i<inner.GetNumCols(); i++) {
      layerThickness_[0][i] = abs(outer[0][i] - inner[0][i]);
      layerThickness_[1][i] = abs(outer[1][i] - inner[1][i]);
    }

    // store local coordinate system
    coordSysId_ = coordSysId;
    coordSys_ = domain->GetCoordSystem(coordSysId_);

  }


} // end namespace CoupledField
