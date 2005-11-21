#include <iostream>
#include <fstream>

#include "laplacePMLInt.hh"

namespace CoupledField
{

//   LaplacePMLInt::LaplacePMLInt(BaseFE * aptelem, Double aVal, Boolean axi)
//     : BaseForm(aptelem),laplVal_ (aVal)
//   {
//     ENTER_FCN( "LaplacePMLInt::LaplacePMLInt");

//     isaxi_ = axi;
//   }


  LaplacePMLInt::LaplacePMLInt(Double damp, Double thickness, Boolean axi)
    : BaseForm(),damping_(damp),layerThickness_(thickness)
  {
    ENTER_FCN( "LaplacePMLInt::LaplacePMLInt" );

    isaxi_ = axi;
  }


 
  LaplacePMLInt::~LaplacePMLInt()
  {
    ENTER_FCN( "LaplacePMLInt::~LaplacePMLInt" );
  }



  void LaplacePMLInt::CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
  {
    ENTER_FCN( "LaplacePMLInt::CalcElementMatrix" );
  
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

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, jacDet);

        xiDx.Transpose(xiDxTransp);

	//multiply the derivatives with the x-,y- and z-factors
	for (UInt i=0; i<xiDx.GetSizeRow(); i++) {
	  for (UInt j=0; j<xiDx.GetSizeCol(); j++) {
	    xiDx[i][j] += factorsPML[i];
	  }
	}

        partElemMat = xiDx * xiDxTransp;
        
        if (isaxi_)
          {
            ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
            CoordAtIP = ptCoord * ShpFncAtIp;
            partElemMat *= 2 * PI * intWeights[actIntPt-1] * jacDet * damping_ * CoordAtIP[0];

          }
        else 
          partElemMat *= intWeights[actIntPt-1] * jacDet * damping_;

        elemMat += partElemMat;
      }

    //    std::cout << "ElemMat:\n" << elemMat << std::endl;
  }


  void LaplacePMLInt::ComputeFactorPML(Vector<Double>& factorsPML, Matrix<Double> & coordMat)
  {
    ENTER_FCN( "LaplacePMLInt ::ComputefactorPML"); 

    //compute average position of element;
    UInt numVals = coordMat.GetSizeRow();
    Vector<Double> pos;
    pos.Resize(numVals);
    pos.Init(0);
    for (UInt i=0; i<coordMat.GetSizeRow(); i++) {
      for (UInt j=0; j<coordMat.GetSizeCol(); j++) {
	pos[i] += coordMat[i][j];
      }
      pos[i] /= (Double) coordMat.GetSizeRow();
    }

    Vector<Complex> factors(numVals);
    Double factor = damping_ / (2*PI*frequency_);
    Double imagVal;

    if ( pos[0] < minMax_[0][0] || pos[0] > minMax_[1][0] ) {
      imagVal = factor*pos[0];
      factors[0] = (1,-imagVal);

      if ( pos[1] < minMax_[0][1] || pos[1] > minMax_[1][1] ) {
	imagVal = factor*pos[1];
	factors[1] = (1,-imagVal);

	//check for 3D
	if (numVals == 3) {
	  if  (pos[1] < minMax_[0][1] || pos[1] > minMax_[1][1] ) {
	    imagVal = factor*pos[2];
	    factors[2] = (1,-imagVal);
	  }
	  else {
	    factors[2] = (1,0);
	  }
	}
      }
      
      else {
	factors[1] = (1,0);

	//check for 3D
	if (numVals == 3) {
	  if  (pos[1] < minMax_[0][1] || pos[1] > minMax_[1][1] ) {
	    imagVal = factor*pos[2];
	    factors[2] = (1,-imagVal);
	  }
	  else {
	    factors[2] = (1,0);
	  }
	}
      }
    }

    else {
      factors[0] = (1,0); 

      if ( pos[1] < minMax_[0][1] || pos[1] > minMax_[1][1] ) {
	imagVal = factor*pos[1];
	factors[1] = (1,-imagVal);

	//check for 3D
	if (numVals == 3) {
	  if  (pos[1] < minMax_[0][1] || pos[1] > minMax_[1][1] ) {
	    imagVal = factor*pos[2];
	    factors[2] = (1,-imagVal);
	  }
	  else {
	    factors[2] = (1,0);
	  }
	}
      }

      else {
	factors[1] = (1,0); 

	//check for 3D
	if (numVals == 3) {
	  if  (pos[1] < minMax_[0][1] || pos[1] > minMax_[1][1] ) {
	    imagVal = factor*pos[2];
	    factors[2] = (1,-imagVal);
	  }
	  else {
	    factors[2] = (1,0);
	  }
	}
      }
    }


    factorsPML.Resize(numVals);
    Complex val;

    if (numVals == 3) {
      if (piezoMatType_ == REALMATERIALPARAMETER) {
	//x-part
	val = factors[1] * factors[2] / factors[0];
	factorsPML[0] = val.real();
	//y-part
	val = factors[0] * factors[2] / factors[1];
	factorsPML[1] = val.real();
	//z-part
	val = factors[0] * factors[1] / factors[2];
	factorsPML[2] = val.real();
      }
      else if (piezoMatType_ == IMAGMATERIALPARAMETER) {
	//x-part
	val = factors[1] * factors[2] / factors[0];
	factorsPML[0] = val.imag();
	//y-part
	val = factors[0] * factors[2] / factors[1];
	factorsPML[1] = val.imag();
	//z-part
	val = factors[0] * factors[1] / factors[2];
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

    std::cout << "x-Pos: " << pos[0] << " val= " << factorsPML[0] << std::endl;
    std::cout << "y-Pos: " << pos[0] << " val= " << factorsPML[1] << std::endl;
    if (numVals == 3) 
      std::cout << "z-Pos: " << pos[0] << " val= " << factorsPML[2] << std::endl;

  }

  void LaplacePMLInt::Print(std::ostream * out, const Matrix<Double> Result) const
  {
    ENTER_FCN( "LaplacePMLInt ::Print"); 
    (*out)<< "Laplace stiffness matrix:" << std::endl << Result;
  }

} // end namespace CoupledField
