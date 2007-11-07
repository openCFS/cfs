// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "pmlBasics.hh"

namespace CoupledField
{

  PMLBasics::PMLBasics( std::string dampingTypePML, Double damp, std::string type)
  {

    //check for correct type
    if ( type == "laplaceInt" ) {
      formsType_ = type;
    }
    else if ( type == "massInt" ) {
      formsType_ = type;
    }
    else {
      Error("PMLInt: type must be laplaceInt or massInt", __FILE__, __LINE__);
    }

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
    else {
      Error("Damping type for PML not known", __FILE__, __LINE__);
    }

    dampingFactor_ = damp;
  }


  PMLBasics::~PMLBasics()
  {
  }


  void PMLBasics::ComputeFactorPML(Vector<Complex>& factorsPML, Complex& pmlDet, 
                                   Vector<Double> & pos, Double omega)
  {

    UInt numVals = pos.GetSize();

    Vector<Complex> factors(numVals);
    Double omegaInv = 1.0 / omega;
    Double imagVal, factor;

    if ( pos[0] < minX_ || pos[0] > maxX_ ) {
      factor = ComputeDampingFactor(pos, X);
      imagVal = factor * omegaInv;
      factors[0] = Complex(1.0,-imagVal);

      if ( pos[1] < minY_ || pos[1] > maxY_ ) {
	factor = ComputeDampingFactor(pos, Y);
	imagVal = factor * omegaInv;
	factors[1] = Complex(1.0,-imagVal);

	//check for 3D
	if (numVals == 3) {
	  if  (pos[2] < minZ_ || pos[2] > maxZ_ ) {
	    //compute z-value
	    factor = ComputeDampingFactor(pos, Z);
	    imagVal = factor * omegaInv;
	    factors[2] = Complex(1.0,-imagVal);
	  }
	  else {
	    factors[2] = Complex(1.0,0);
	  }
	}
      }
      
      else {
	factors[1] = Complex(1.0,0);

	//check for 3D
	if (numVals == 3) {
	  if  (pos[2] < minZ_ || pos[2] > maxZ_ ) {
	    //compute z-value
	    factor  = ComputeDampingFactor(pos, Z);
	    imagVal = factor * omegaInv;
	    factors[2] = Complex(1.0,-imagVal);
	  }
	  else {
	    factors[2] = Complex(1.0,0);
	  }
	}
      }
    }

    else {
      factors[0] = Complex(1.0,0.0); 

      if ( pos[1] < minY_ || pos[1] > maxY_ ) {
	//compute y-value
	factor  = ComputeDampingFactor(pos, Y);
	imagVal = factor * omegaInv;
	factors[1] = Complex(1.0,-imagVal);

	//check for 3D
	if (numVals == 3) {
	  if  (pos[2] < minZ_ || pos[2] > maxZ_ ) {
	    //compute z-value
	    factor  = ComputeDampingFactor(pos, Z);
	    imagVal = factor * omegaInv;
	    factors[2] = Complex(1.0,-imagVal);
	  }
	  else {
	    factors[2] = Complex(1.0,0);
	  }
	}
      }

      else {
	factors[1] = Complex(1.0,0); 

	//check for 3D
	if (numVals == 3) {
	  if  (pos[2] < minZ_ || pos[2] > maxZ_ ) {
	    factor  = ComputeDampingFactor(pos, Z);
	    imagVal = factor * omegaInv;
	    factors[2] = Complex(1,-imagVal);
	  }
	  else {
	    factors[2] = Complex(1.0,0);
	  }
	}
      }
    }


//      std::cout << "pos:\n" << pos << std::endl;
//      std::cout << "factors:\n" << factors << std::endl;

    factorsPML.Resize(numVals);
    factorsPML.Init();

    if ( formsType_ == "laplaceInt") {
      if (numVals == 3) {
        factorsPML[0] = Complex(1.0,0) / factors[0];  
        //y-part
        factorsPML[1] = Complex(1.0,0) / factors[1];   
        //z-part
        factorsPML[2] = Complex(1.0,0) / factors[2];  

        // determinant
        pmlDet = factors[0]*factors[1]*factors[2];
      }
      
      else {
        factorsPML[0] = Complex(1.0,0) / factors[0];  
        //y-part
        factorsPML[1] = Complex(1.0,0) / factors[1]; 
        //determinant
        pmlDet = factors[0]*factors[1];
      }
    }
    else {
      //PML factor for mass integrator, we just need the determinant
      if (numVals == 3) {
	pmlDet = factors[0]*factors[1]*factors[2];
      }     
      else {
	pmlDet = factors[0]*factors[1];
      }

    }

  }


  Double PMLBasics::ComputeDampingFactor(Vector<Double>& pos, Directions dir)
  {

    Double factor, maxPos, delta, diffCoord;
    UInt idx;

    factor = 1.0;
    if ( dampingTypePML_ == "constant" ) {
      factor = dampingFactor_;
    }

    else if ( dampingTypePML_ == "quadDist" ) {
      if ( dir == X ) {
	//get correct layer thickness
	if ( pos[0] < minX_ ) {
	  delta = layerThickness_[0][0];
	  diffCoord = abs(pos[0]) - minX_;
	}
	else {
	  delta = layerThickness_[1][0];
	  diffCoord = abs(pos[0]) - maxX_;
	}
      }
      else if ( dir == Y ) {
	//get correct layer thickness
	if ( pos[1] < minY_ ) {
	  delta = layerThickness_[0][1];
	  diffCoord = abs(pos[1]) - minY_;
	}
	else {
	  delta = layerThickness_[1][1];
	  diffCoord = abs(pos[1]) - maxY_;
	}
      }
      else {
	//get correct layer thickness
	if ( pos[2] < minZ_ ) {
	  delta = layerThickness_[0][2];
	  diffCoord = abs(pos[2]) - minZ_;
	}
	else {
	  delta = layerThickness_[1][2];
	  diffCoord = abs(pos[2]) - maxZ_;
	}
      }

      factor = dampingFactor_ * ( diffCoord*diffCoord )/ ( delta*delta );

    }

    else if ( dampingTypePML_ == "inverseDist" ) {
      if ( dir == X ) {
	//get correct maximal PML y-coordinate
	if ( pos[0] < minX_ ) {
	  maxPos = minX_ - layerThickness_[0][0];
	}
	else {
	  maxPos = maxX_ + layerThickness_[1][0];
	}
	idx = 0;
      }

      else if ( dir == Y ) {
	//get correct maximal PML y-coordinate
	if ( pos[1] < minY_ ) {
	  maxPos = minY_ - layerThickness_[0][1];
	}
	else {
	  maxPos = maxY_ + layerThickness_[1][1];
	}
	idx = 1;
      }

      else {
	//get correct maximal PML z-coordinate
	if ( pos[2] < minZ_ ) {
	  maxPos = minZ_ - layerThickness_[0][2];
	}
	else {
	  maxPos = maxZ_ + layerThickness_[1][2];
	}
	idx = 2;
      }

      if ( abs (maxPos - pos[idx]) < 1e-12 ) {
	Error("PML damping inverseDist divides by factor smaller 1E-12",
	      __FILE__,__LINE__);
      }

      //      std::cout << "maxPos =" << maxPos << std::endl;

      factor = abs (dampingFactor_ / ( maxPos  - pos[idx] ) );
    }

    return factor;
  }


  void PMLBasics:: SetPosPML(Matrix<Double> & inner, Matrix<Double> & outer)
  {

    // inner/outer:   xmin  ymin  zmin
    //                xmax  ymax  zmax

    minX_ = inner[0][0];
    maxX_ = inner[1][0];
    minY_ = inner[0][1];
    maxY_ = inner[1][1];

    if (inner.GetSizeCol() > 2 ) {
      minZ_ = inner[0][2];
      maxZ_ = inner[1][2];
    }

    //get layer thickness
    layerThickness_.Resize(2,inner.GetSizeCol());
    for (UInt i=0; i<inner.GetSizeCol(); i++) {
      layerThickness_[0][i] = abs(outer[0][i] - inner[0][i]);
      layerThickness_[1][i] = abs(outer[1][i] - inner[1][i]);
    }

    std::cout << "LayerThickness:\n" << layerThickness_ << std::endl;

  }





} // end namespace CoupledField
