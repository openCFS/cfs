// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "FlatShellElecInt.hh"
#define PI 3.141592654

namespace CoupledField {

 void FlatShellElecInt::CalcElementMatrix( Matrix<Double>& elemMat,
					   EntityIterator& ent1, 
					   EntityIterator& ent2 ) {

    
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    
    //std::cout << "FlatShellStiffInt::CalcElementMatrix\n" << std::endl;

    const UInt nrLayers  = composite_->thickness.GetSize();
    //std::cout << "The Number of Layers is\n" << nrLayers << std::endl;
        
    Matrix<Double> eMat; //E matrix
    Matrix<Double> TransMat; //Transformation matrix
    Matrix<Double> ShellCoord; //Shell Coordinates
    
    double Omega = 0.0;

    //Element matrix in the case of Flat Shells with linear elements is 24x24
    elemMat.Resize(nrLayers);   
    elemMat.Init();
       
    FlatShellInt::CoordTrans(ptCoord_, TransMat, ShellCoord );
    
    FlatShellElecInt::calcEMat(eMat);
    
    Omega = ptelem->CalcVolume(ShellCoord,FALSE);    
    //See equation 4.65 Pietfort
    elemMat = eMat*(- Omega);

    // Do the back-Rotation and return
    // Back Rotation is not needed, as this matrix does not depend on 
    // shell coordinates! (Andreas)
    //FlatShellInt::LocaltoGlob(elemMat,TransMat );

    //std::cerr << "elemMat =\n" << elemMat << std::endl << std::endl;
        
  }
  


  
  void FlatShellElecInt::calcEMat( Matrix<Double> &eMat ){


    const UInt nrLayers  = composite_->thickness.GetSize();
    Double thickness_all = 0.0;
    //Piezoelectric coefficients e31 and e32
    Double epsillon_k = 0.0; 
    
    eMat.Resize(nrLayers);
    eMat.Init();
    	  
    //Initialisation of the composite structure variables
    z_.Resize(nrLayers + 1);
    z_.Init();
    orAngle_.Resize(nrLayers + 1);
    orAngle_.Init();

    for (UInt k=0; k < nrLayers; k++) {
        thickness_all += composite_->thickness[k];
        }
    z_[0]= - (thickness_all/2.0);
    orAngle_[0]= 0.0;
    
    for(UInt i=1; i <= nrLayers; i++) {
    	z_[i] = z_[i-1] + composite_->thickness[i-1];
        orAngle_[i] = composite_->orientation[i-1]*PI/180;
        }
     
     //std::cout << "z coordinates are \n" << z_ << std::endl;
     //std::cout << "orientations are are \n" << orAngle_ << std::endl;
                       
    for (UInt k=1 ; k <= nrLayers ; k++) {
        
        composite_->materials[k-1]->GetScalar(epsillon_k,ELEC_PERMITTIVITY,REAL);
        eMat[k-1][k-1] = epsillon_k/composite_->thickness[k-1]; 
	}
        //std::cout << "The E matrix is\n" << eMat << std::endl;
    
  }

  // *************************************
  //   Constructor for Composite Material
  // *************************************
  FlatShellElecInt::FlatShellElecInt( Composite * composite ) 
    : FlatShellInt(composite) {

    name_ = "FlatShellElecInt";
    
    baseType_ = STIFFNESS;             
    }

  // **************
  //   Destructor
  // **************
  FlatShellElecInt::~FlatShellElecInt() {
  }
  
} // end namespace CoupledField
