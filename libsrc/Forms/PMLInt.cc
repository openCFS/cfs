#include <iostream>
#include <fstream>

#include "PMLInt.hh"

namespace CoupledField
{

  PMLInt::PMLInt(std::string type, Double factor, std::string dampingTypePML, Double damp, 
		 Double thickness, Boolean axi)
    : BaseForm(), formsFactor_(factor), formsType_(type), dampingFactor_(damp), 
      layerThickness_(thickness)
  {
    ENTER_FCN( "PMLInt::PMLInt" );

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
    else {
      Error("Damping type for PML not known", __FILE__, __LINE__);
    }

    isaxi_ = axi;

    //    std::cout << "damp=" << damp << std::endl;
  }


 
  PMLInt::~PMLInt()
  {
    ENTER_FCN( "PMLInt::~PMLInt" );
  }



  void PMLInt::CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
  {
    ENTER_FCN( "PMLInt::CalcElementMatrix" );
  
    if ( formsType_ =="laplaceInt" ) {
      CalcElementMatrixStiff(ptCoord, elemMat);
    }
    else {
      CalcElementMatrixMass(ptCoord, elemMat);
    }
    
  }


  void PMLInt::CalcElementMatrixStiff(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
  {
    ENTER_FCN( "PMLInt::CalcElementMatrixStiff" );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Double> xiDxTransp;
    Matrix<Double> partElemMat;
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes); elemMat.Init();

    Vector<Double> factorsPML;
    ComputeFactorPML( factorsPML, ptCoord);

    //    std::cout << "in Forms:\n" << factorsPML << std::endl;

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, jacDet);

        xiDx.Transpose(xiDxTransp);

	//multiply the derivatives with the x-,y- and z-factors
	for (UInt i=0; i<xiDx.GetSizeCol(); i++) {
	  for (UInt j=0; j<xiDx.GetSizeRow(); j++) {
	    xiDx[j][i] *= factorsPML[i];
	  }
	}

        partElemMat = xiDx * xiDxTransp;
        
        if (isaxi_)
          {
            ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
            CoordAtIP = ptCoord * ShpFncAtIp;
            partElemMat *= 2 * PI * intWeights[actIntPt-1] * jacDet * formsFactor_ * CoordAtIP[0];

          }
        else 
          partElemMat *= intWeights[actIntPt-1] * jacDet * formsFactor_;

        elemMat += partElemMat;
      }

    //    std::cout << "PML-ElemMatSiff:\n" << elemMat << std::endl;
  }


  void PMLInt::CalcElementMatrixMass(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
  {
    ENTER_FCN( "PMLInt::CalcElementMatrixMass" );
  
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes);
    
    Vector<Double> factorsPML;
    ComputeFactorPML( factorsPML, ptCoord);

    Double factorPML = factorsPML[0] * formsFactor_;

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
        
      ptelem-> GetShFncAtIp(shapeFncAtIp, actIntPt);
        
      partElemMat.DyadicMult(shapeFncAtIp);
        
      if (isaxi_) {
        CoordAtIP = ptCoord * shapeFncAtIp;
        partElemMat *= 2 * PI * intWeights[actIntPt-1] * factorPML * jacDet * CoordAtIP[0];
      }
      else 
        partElemMat *= intWeights[actIntPt-1] * factorPML * jacDet;
        
      elemMat += partElemMat;
    }

    //    std::cout << "PML-ElemMatMass:\n" << elemMat << std::endl;

  }


  void PMLInt::ComputeFactorPML(Vector<Double>& factorsPML, Matrix<Double> & coordMat)
  {
    ENTER_FCN( "PMLInt ::ComputefactorPML"); 

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

    Vector<Complex> factors(numVals);
    Double omegaInv = 1.0 / frequency_;
    Double imagVal;

    if ( pos[0] < minX_ || pos[0] > maxX_ ) {
      imagVal = dampingFactor_ * omegaInv;
      factors[0] = Complex(1.0,-imagVal);

      if ( pos[1] < minY_ || pos[1] > maxY_ ) {
	imagVal = dampingFactor_ * omegaInv;
	factors[1] = Complex(1.0,-imagVal);

	//check for 3D
	if (numVals == 3) {
	  if  (pos[2] < minZ_ || pos[2] > maxZ_ ) {
	    //compute z-value
	    imagVal = dampingFactor_ * omegaInv;
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
	  if  (pos[2] < minZ_ || pos[1] > maxZ_ ) {
	    //compute z-value
	    imagVal = dampingFactor_ * omegaInv;
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
	imagVal = dampingFactor_ * omegaInv;
	factors[1] = Complex(1.0,-imagVal);

	//check for 3D
	if (numVals == 3) {
	  if  (pos[2] < minZ_ || pos[2] > maxZ_ ) {
	    //compute z-value
	    imagVal = dampingFactor_ * omegaInv;
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
	    imagVal = dampingFactor_ * omegaInv;
	    factors[2] = Complex(1,-imagVal);
	  }
	  else {
	    factors[2] = Complex(1.0,0);
	  }
	}
      }
    }


    factorsPML.Resize(numVals);
    Complex val;

    if ( formsType_ == "laplaceInt") {
      if (numVals == 3) {

	if (piezoMatType_ == REALMATERIALPARAMETER) {
	  //x-part
	  val = factors[1]*factors[2] / factors[0];
	  factorsPML[0] = val.real();
	  //y-part
	  val = factors[0]*factors[2] / factors[1];
	  factorsPML[1] = val.real();
	  //z-part
	  val = factors[0]*factors[1] / factors[2];
	  factorsPML[2] = val.real();
	}
	else if (piezoMatType_ == IMAGMATERIALPARAMETER) {
	  //x-part
	  val = factors[1]*factors[2] / factors[0];
	  factorsPML[0] = val.imag();
	  //y-part
	  val = factors[0]*factors[2] / factors[1];
	  factorsPML[1] = val.imag();
	  //z-part
	  val = factors[0]*factors[1] / factors[2];
	  factorsPML[2] = val.imag();
	}
      }
      
      else {
	if (piezoMatType_ == REALMATERIALPARAMETER) {
	  //x-part
	  val = factors[1] / factors[0];
	  factorsPML[0] = val.real();
	  //y-part
	  val = factors[0] / factors[1];
	  factorsPML[1] = val.real();
	}
	else if (piezoMatType_ == IMAGMATERIALPARAMETER) {
	  //x-part
	  val = factors[1] / factors[0];
	  factorsPML[0] = val.imag();
	  //y-part
	  val = factors[0] / factors[1];
	  factorsPML[1] = val.imag();
	}
      }
    }
    else {
      //PML factor for mass integrator
      if (numVals == 3) {
	val = factors[0]*factors[1]*factors[2];
	if (piezoMatType_ == REALMATERIALPARAMETER) {
	  factorsPML[0] = val.real();
	}
	else if (piezoMatType_ == IMAGMATERIALPARAMETER) {
	  factorsPML[0] = val.imag();
	}
      }
      
      else {
	val = factors[0]*factors[1];
	if (piezoMatType_ == REALMATERIALPARAMETER) {
	  factorsPML[0] = val.real();
	}
	else if (piezoMatType_ == IMAGMATERIALPARAMETER) {
	  factorsPML[0] = val.imag();
	}
      }

    }

  }


//   void PMLInt::ComputeFactorPML(Vector<Double>& factorsPML, Matrix<Double> & coordMat)
//   {
//     ENTER_FCN( "PMLInt ::ComputefactorPML"); 

//     //compute average position of element;
//     UInt numVals = coordMat.GetSizeRow();
//     Vector<Double> pos;
//     pos.Resize(numVals);
//     pos.Init(0);
//     for (UInt i=0; i<coordMat.GetSizeRow(); i++) {
//       for (UInt j=0; j<coordMat.GetSizeCol(); j++) {
// 	pos[i] += coordMat[i][j];
//       }
//       pos[i] /= (Double) coordMat.GetSizeCol();
//     }

//     Vector<Complex> factors(numVals);
//     Double omegaInv = 1.0 / frequency_;
//     Double imagVal;

//     if ( pos[0] < minX_ || pos[0] > maxX_ ) {
//       //compute x-factor
//       if (  pos[0] < minX_ ) {
// 	imagVal = dampingFactor_ * omegaInv * (pos[0]-minX_) * (pos[0]-minX_);
//       }
//       else {
// 	imagVal = dampingFactor_ * omegaInv * (pos[0]-maxX_) * (pos[0]-maxX_);
//       }
//       factors[0] = Complex(1.0,-imagVal);

//       if ( pos[1] < minY_ || pos[1] > maxY_ ) {
// 	//compute y-value
// 	if (  pos[1] < minY_ ) {
// 	  imagVal = dampingFactor_ * omegaInv * (pos[1]-minY_) * (pos[1]-minY_);
// 	}
// 	else {
// 	  imagVal = dampingFactor_ * omegaInv * (pos[1]-maxY_) * (pos[1]-maxY_);
// 	}
// 	factors[1] = Complex(1.0,-imagVal);

// 	//check for 3D
// 	if (numVals == 3) {
// 	  if  (pos[2] < minZ_ || pos[2] > maxZ_ ) {
// 	    //compute z-value
// 	    if (  pos[2] < minZ_ ) {
// 	      imagVal = dampingFactor_ * omegaInv * (pos[2]-minZ_) * (pos[2]-minZ_);
// 	    }
// 	    else {
// 	      imagVal = dampingFactor_ * omegaInv * (pos[2]-maxZ_) * (pos[2]-maxZ_);
// 	    }
// 	    factors[2] = Complex(1.0,-imagVal);
// 	  }
// 	  else {
// 	    factors[2] = Complex(1.0,0);
// 	  }
// 	}
//       }
      
//       else {
// 	factors[1] = Complex(1.0,0);

// 	//check for 3D
// 	if (numVals == 3) {
// 	  if  (pos[2] < minZ_ || pos[1] > maxZ_ ) {
// 	    //compute z-value
// 	    if (  pos[2] < minZ_ ) {
// 	      imagVal = dampingFactor_ * omegaInv * (pos[2]-minZ_) * (pos[2]-minZ_);
// 	    }
// 	    else {
// 	      imagVal = dampingFactor_ * omegaInv * (pos[2]-maxZ_) * (pos[2]-maxZ_);
// 	    }
// 	    factors[2] = Complex(1.0,-imagVal);
// 	  }
// 	  else {
// 	    factors[2] = Complex(1.0,0);
// 	  }
// 	}
//       }
//     }

//     else {
//       factors[0] = Complex(1.0,0.0); 

//       if ( pos[1] < minY_ || pos[1] > maxY_ ) {
// 	//compute y-value
// 	if (  pos[1] < minY_ ) {
// 	  imagVal = dampingFactor_ * omegaInv * (pos[1]-minY_) * (pos[1]-minY_);
// 	}
// 	else {
// 	  imagVal = dampingFactor_ * omegaInv * (pos[1]-maxY_) * (pos[1]-maxY_);
// 	}
// 	factors[1] = Complex(1.0,-imagVal);

// 	//check for 3D
// 	if (numVals == 3) {
// 	  if  (pos[2] < minZ_ || pos[2] > maxZ_ ) {
// 	    //compute z-value
// 	    if (  pos[2] < minZ_ ) {
// 	      imagVal = dampingFactor_ * omegaInv * (pos[2]-minZ_) * (pos[2]-minZ_);
// 	    }
// 	    else {
// 	      imagVal = dampingFactor_ * omegaInv * (pos[2]-maxZ_) * (pos[2]-maxZ_);
// 	    }
// 	    factors[2] = Complex(1.0,-imagVal);
// 	  }
// 	  else {
// 	    factors[2] = Complex(1.0,0);
// 	  }
// 	}
//       }

//       else {
// 	factors[1] = Complex(1.0,0); 

// 	//check for 3D
// 	if (numVals == 3) {
// 	  if  (pos[2] < minZ_ || pos[2] > maxZ_ ) {
// 	    //compute z-value
// 	    if (  pos[2] < minZ_ ) {
// 	      imagVal = dampingFactor_ * omegaInv * (pos[2]-minZ_) * (pos[2]-minZ_);
// 	    }
// 	    else {
// 	      imagVal = dampingFactor_ * omegaInv * (pos[2]-maxZ_) * (pos[2]-maxZ_);
// 	    }
// 	    factors[2] = Complex(1,-imagVal);
// 	  }
// 	  else {
// 	    factors[2] = Complex(1.0,0);
// 	  }
// 	}
//       }
//     }


//     factorsPML.Resize(numVals);
//     Complex val;

//     if ( formsType_ == "laplaceInt") {
//       if (numVals == 3) {

// 	if (piezoMatType_ == REALMATERIALPARAMETER) {
// 	  //x-part
// 	  val = factors[1]*factors[2] / factors[0];
// 	  factorsPML[0] = val.real();
// 	  //y-part
// 	  val = factors[0]*factors[2] / factors[1];
// 	  factorsPML[1] = val.real();
// 	  //z-part
// 	  val = factors[0]*factors[1] / factors[2];
// 	  factorsPML[2] = val.real();
// 	}
// 	else if (piezoMatType_ == IMAGMATERIALPARAMETER) {
// 	  //x-part
// 	  val = factors[1]*factors[2] / factors[0];
// 	  factorsPML[0] = val.imag();
// 	  //y-part
// 	  val = factors[0]*factors[2] / factors[1];
// 	  factorsPML[1] = val.imag();
// 	  //z-part
// 	  val = factors[0]*factors[1] / factors[2];
// 	  factorsPML[2] = val.imag();
// 	}
//       }
      
//       else {
// 	if (piezoMatType_ == REALMATERIALPARAMETER) {
// 	  //x-part
// 	  val = factors[1] / factors[0];
// 	  factorsPML[0] = val.real();
// 	  //y-part
// 	  val = factors[0] / factors[1];
// 	  factorsPML[1] = val.real();
// 	}
// 	else if (piezoMatType_ == IMAGMATERIALPARAMETER) {
// 	  //x-part
// 	  val = factors[1] / factors[0];
// 	  factorsPML[0] = val.imag();
// 	  //y-part
// 	  val = factors[0] / factors[1];
// 	  factorsPML[1] = val.imag();
// 	}
//       }

//       std::cout << "x-Pos: " << pos[0] << " val= " << factorsPML[0] << std::endl;
//       std::cout << "y-Pos: " << pos[1] << " val= " << factorsPML[1] << std::endl;
//       if (numVals == 3) 
// 	std::cout << "z-Pos: " << pos[0] << " val= " << factorsPML[2] << std::endl;
//     }

//     else {
//       //PML factor for mass integrator
//       if (numVals == 3) {
// 	val = factors[0]*factors[1]*factors[2];
// 	if (piezoMatType_ == REALMATERIALPARAMETER) {
// 	  factorsPML[0] = val.real();
// 	}
// 	else if (piezoMatType_ == IMAGMATERIALPARAMETER) {
// 	  factorsPML[0] = val.imag();
// 	}
//       }
      
//       else {
// 	val = factors[0]*factors[1];
// 	if (piezoMatType_ == REALMATERIALPARAMETER) {
// 	  factorsPML[0] = val.real();
// 	}
// 	else if (piezoMatType_ == IMAGMATERIALPARAMETER) {
// 	  factorsPML[0] = val.imag();
// 	}
//       }

//       std::cout << "x-Pos: " << pos[0] << " y-Pos=" << pos[1] << std::endl;
//       std::cout << "val= " << factorsPML[0] << std::endl;
//       std::cout << "Factors:\n " << factors << std::endl;
//     }

//   }

  void PMLInt:: SetPosPML(Matrix<Double> & pos)
  {
    ENTER_FCN( "PMLInt ::SetPosXML"); 

    minX_ = pos[0][0];
    maxX_ = pos[1][0];
    minY_ = pos[0][1];
    maxY_ = pos[1][1];

    if (pos.GetSizeCol() > 2 ) {
      minZ_ = pos[0][2];
      maxZ_ = pos[1][2];
    }

  }

  void PMLInt::Print(std::ostream * out, const Matrix<Double> Result) const
  {
    ENTER_FCN( "PMLInt ::Print"); 
    (*out)<< "Laplace stiffness matrix:" << std::endl << Result;
  }


} // end namespace CoupledField
