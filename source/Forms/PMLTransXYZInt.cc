// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "PMLTransXYZInt.hh"

namespace CoupledField
{

  PMLTransXYZInt::PMLTransXYZInt(std::string type, Double factor, std::string dampingTypePML, Double damp, 
				 const UInt nrDofsPerNode, bool axi)
    : BaseForm( NULL ), 
      nrDofsPerNode_(nrDofsPerNode), formsType_(type),formsFactor_(factor),  dampingFactor_(damp)
       
  {
    ENTER_FCN( "PMLTransXYZInt::PMLTransXYZInt" );

    name_ = "PMLTransXYZInt";
    
    //check for correct type
    if ( type == "stiffness" ) {
      formsType_ = type;
    }
    else if ( type == "damping" ) {
      formsType_ = type;
    }
    else {
      Error("PMLTransXYZInt: type must be stiffness or damping", __FILE__, __LINE__);
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

    isaxi_ = axi;
  }


 
  PMLTransXYZInt::~PMLTransXYZInt()
  {
    ENTER_FCN( "PMLTransXYZInt::~PMLTransXYZInt" );
  }



  void PMLTransXYZInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                          EntityIterator& ent1, 
                                          EntityIterator& ent2 )
  {
    ENTER_FCN( "PMLTransXYZInt::CalcElementMatrix" );
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );


    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat, singleMat;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs*nrDofsPerNode_);
    elemMat.Init();
    singleMat.Resize(numFncs); 
    singleMat.Init();

    Vector<Double> factorsPML;
    ComputeFactorPML( factorsPML, ptCoord_);

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem() );
        
      ptelem-> GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );
        
      partElemMat.DyadicMult(shapeFncAtIp);
        
      if (isaxi_) {
        CoordAtIP = ptCoord_ * shapeFncAtIp;
        partElemMat *= 2 * PI * intWeights[actIntPt-1] * jacDet * CoordAtIP[0];
      }
      else 
        partElemMat *= intWeights[actIntPt-1] * jacDet;
        
      singleMat += partElemMat;
    }

    //    std::cout << "PML-SingleMatrix:\n" << singleMat << std::endl;

    Double factor;
    UInt   startCol, startRow;

    // put x-part
    factor   = factorsPML[0] * formsFactor_;
    startCol = 0;
    startRow = 0;
    for (UInt i=0; i < numFncs; i++) {
      for (UInt j=0; j < numFncs; j++) {
	elemMat[i*nrDofsPerNode_ + startRow][j*nrDofsPerNode_ + startCol] = factor*singleMat[i][j]; 
      }
    }

    //    std::cout << "PML-ElemMatrix_X:" << factorsPML[0] << "\n" << elemMat << std::endl;

    // put y-part
    factor   = factorsPML[1] * formsFactor_;
    startCol = 1;
    startRow = 1;
    for (UInt i=0; i < numFncs; i++) {
      for (UInt j=0; j < numFncs; j++) {
	elemMat[i*nrDofsPerNode_ + startRow][j*nrDofsPerNode_ + startCol] = factor*singleMat[i][j]; 
      }
    }

//     std::cout << "Type:" << formsType_ << std::endl;
//     std::cout << "formFactor=" << formsFactor_ << std::endl;
//     std::cout << "PML-ElemMatrix_XY:" << "fx=" << factorsPML[0] << " fy=" 
//  	      << factorsPML[1] << "\n" << elemMat << std::endl;

    if (nrDofsPerNode_ == 3) {
      // put z-part
      factor   = factorsPML[2] * formsFactor_;
      startCol = 2;
      startRow = 2;
      for (UInt i=0; i < numFncs; i++) {
	for (UInt j=0; j < numFncs; j++) {
	  elemMat[i*nrDofsPerNode_ + startRow][j*nrDofsPerNode_ + startCol] = factor*singleMat[i][j]; 
	}
      }
    }

    //    std::cout << "PML-ElemMattrix:\n" << elemMat << std::endl;

  }


  void PMLTransXYZInt::ComputeFactorPML(Vector<Double>& factorsPML, Matrix<Double> & coordMat)
  {
    ENTER_FCN( "PMLTransXYZInt ::ComputefactorPML"); 

    UInt numVals = coordMat.GetSizeRow();
    Vector<Double> pos;
    pos.Resize(numVals);
    pos.Init(0);
    for (UInt i=0; i<coordMat.GetSizeRow(); i++) {
      for (UInt j=0; j<coordMat.GetSizeCol(); j++) {
	pos[i] += coordMat[i][j];
      }
      pos[i] /= (Double) coordMat.GetSizeCol();
    }

    Vector<Double> factors(numVals);

    if ( pos[0] < minX_ || pos[0] > maxX_ ) {
      factors[0] = ComputeDampingFactor(pos, X);

      if ( pos[1] < minY_ || pos[1] > maxY_ ) {
	factors[1] = ComputeDampingFactor(pos, Y);

	//check for 3D
	if (numVals == 3) {
	  if  (pos[2] < minZ_ || pos[2] > maxZ_ ) {
	    //compute z-value
	    factors[2] = ComputeDampingFactor(pos, Z);
	  }
	  else {
	    factors[2] = 0.0;
	  }
	}
      }
      
      else {
	factors[1] = 0.0;

	//check for 3D
	if (numVals == 3) {
	  if  (pos[2] < minZ_ || pos[2] > maxZ_ ) {
	    //compute z-value
	    factors[2] = ComputeDampingFactor(pos, Z);
	  }
	  else {
	    factors[2] = 0.0;
	  }
	}
      }
    }

    else {
      factors[0] = 0.0;

      if ( pos[1] < minY_ || pos[1] > maxY_ ) {
	//compute y-value
	factors[1] = ComputeDampingFactor(pos, Y);

	//check for 3D
	if (numVals == 3) {
	  if  (pos[2] < minZ_ || pos[2] > maxZ_ ) {
	    //compute z-value
	    factors[2] = ComputeDampingFactor(pos, Z);
	  }
	  else {
	    factors[2] = 0.0;
	  }
	}
      }

      else {
	factors[1] = 0.0;

	//check for 3D
	if (numVals == 3) {
	  if  (pos[2] < minZ_ || pos[2] > maxZ_ ) {
	    factors[2] = ComputeDampingFactor(pos, Z);
	  }
	  else {
	    factors[2] = 0.0;
	  }
	}
      }
    }


//     std::cout << "pos:\n" << pos << std::endl;
//     std::cout << "factors:\n" << factors << std::endl;

    factorsPML.Resize(numVals);
    factorsPML.Init();

    if ( formsType_ == "stiffness") {
      //x-part
      factorsPML[0] = factors[0]*factors[0];
      //y-part
      factorsPML[1] = factors[1]*factors[1];
      if (numVals == 3) {
	//z-part
	factorsPML[2] = factors[2]*factors[2];
      }
    }
    else {
      //PML factor for mass integrator
      //x-part
      factorsPML[0] = 2.0*factors[0];
      //y-part
      factorsPML[1] = 2.0*factors[1];
      if (numVals == 3) {
	//z-part
	factorsPML[2] = 2.0*factors[2];
      }
    }

  }





  Double PMLTransXYZInt:: ComputeDampingFactor(Vector<Double>& pos, Directions dir)
  {
    ENTER_FCN( "PMLTransXYZInt :: ComputeDampingFactor"); 

    Double factor, maxPos, delta, diffCoord;
    UInt idx;

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


  void PMLTransXYZInt:: SetPosPML(Matrix<Double> & inner, Matrix<Double> & outer)
  {
    ENTER_FCN( "PMLTransXYZInt ::SetPosXML"); 

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

    //    std::cout << "LayerThickness:\n" << layerThickness_ << std::endl;

  }


} // end namespace CoupledField
