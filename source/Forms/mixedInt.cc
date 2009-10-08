#include <iostream>
#include <fstream>

#include "mixedInt.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"


namespace CoupledField {

//=========================== Mass-MPP ==========================//

  MassMixedInt_PP::MassMixedInt_PP(const Double aFactor, bool axi, bool coordUpdate )
    : BaseForm(NULL), 
      factor_(aFactor)     
  {
    name_ = "MassMixedPressureInt";

    isaxi_ = axi;
    coordUpdate_ = coordUpdate;
    baseType_ = MASS;
  }

 
  MassMixedInt_PP::~MassMixedInt_PP()
  {
  }

 


  void MassMixedInt_PP::CalcElementMatrix( Matrix<Double>& elemMat,
                                   EntityIterator& ent1,
                                   EntityIterator& ent2  )  {

     // Extract pointer to reference element and get coordinates
     Vector<Double> CoordAtIP;

     ExtractElemInfo( ent1 );
     ptelem->SetAnsatzFct( ansatzFct1_ );
     UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
     const UInt nrIntPts= ptelem->GetNumIntPoints();    
     const Vector<Double> & intWeights = ptelem->GetIntWeights();
     const Vector<Double> * intPoints = ptelem->GetIntPoints();

     Double jacDet;
     Vector<Double> shapeFncAtIp;
     Matrix<Double> partElemMat;
     
     elemMat.Resize(numFncs);
     elemMat.Init();
     partElemMat.Resize(numFncs);
     partElemMat.Init();
     
     for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
       
       jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem() );
       
       ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );
       
       partElemMat.DyadicMult(shapeFncAtIp);
       if (isaxi_) {
	        ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
				    ptCoord_, it1_.GetElem() );
	        partElemMat *= 2 * PI * intWeights[actIntPt-1] * factor_
	                       * jacDet * CoordAtIP[0];
       }
       else {
	        partElemMat *= intWeights[actIntPt-1] * factor_ * jacDet;
       }
       elemMat += partElemMat;
     }
     
     //std::cout << "Part Element Mass Matrix, MPP:\n" << elemMat << std::endl;
  }



//=========================== Mass-MVV ==========================//

  MassMixedInt_VV::MassMixedInt_VV(const Double aFactor, UInt dim, bool axi, bool coordUpdate )
    : BaseForm(NULL), 
      factor_(aFactor)
  {
    name_ = "MassMixedVelocityInt";

    isaxi_ = axi;
    dim_   = dim;
    coordUpdate_ = coordUpdate;
    baseType_    = MASS;
  }

 
  MassMixedInt_VV::~MassMixedInt_VV()
  {
  }

 


  void MassMixedInt_VV::CalcElementMatrix( Matrix<Double>& elemMat,
                                   EntityIterator& ent1,
                                   EntityIterator& ent2  )  {
     // Extract pointer to reference element and get coordinates
     Vector<Double> CoordAtIP;
     Matrix<Double> tempElemMat;

     ExtractElemInfo( ent1 );
     ptelem->SetAnsatzFct( ansatzFct1_ );
     UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
     const UInt nrIntPts= ptelem->GetNumIntPoints();    
     const Vector<Double> & intWeights = ptelem->GetIntWeights();
     const Vector<Double> * intPoints = ptelem->GetIntPoints();
     Double jacDet;

     elemMat.Resize(dim_*numFncs);
     elemMat.Init();

     tempElemMat.Resize(numFncs);
     tempElemMat.Init();

     Matrix<Double>  partElemMat;
     partElemMat.Resize(numFncs);
     partElemMat.Init();

     Vector<Double> shapeFncAtIp;

     for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
       
       jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem() );       
       ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );
       partElemMat.DyadicMult(shapeFncAtIp);

       if (isaxi_) {
	 ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
				    ptCoord_, it1_.GetElem() );
	 partElemMat *= 2 * PI * intWeights[actIntPt-1] * factor_
	   * jacDet * CoordAtIP[0];
       }
       else {
         partElemMat *= intWeights[actIntPt-1] * factor_ * jacDet;
       }

       tempElemMat += partElemMat;
     }


     //Blowing up the element matrix to number of velocity components 
     const UInt singleDofSize = tempElemMat.GetNumRows();
     for (UInt i=0; i < singleDofSize; i++)
       for (UInt j=0; j < singleDofSize; j++)
         for (UInt actDof = 0; actDof < dim_ ; actDof++)
           elemMat[i*dim_ + actDof][j*dim_ + actDof] = tempElemMat[i][j]; 

     //std::cout << "ElemMatMass:\n" << elemMat << std::endl;

   }


//=========================== Mass-Piola-MVV ==========================//

  MassPiolaMixedInt_VV::MassPiolaMixedInt_VV(const Double aFactor, bool axi, bool coordUpdate )
    : BaseForm(NULL), 
      factor_(aFactor)
  {
    name_ = "MassMixedVelocityInt with PIOLA";

    isaxi_ = axi;
    coordUpdate_ = coordUpdate;
    baseType_ = MASS;
  }

 
  MassPiolaMixedInt_VV::~MassPiolaMixedInt_VV()
  {
  }

 


  void MassPiolaMixedInt_VV::CalcElementMatrix( Matrix<Double>& elemMat,
                                   EntityIterator& ent1,
                                   EntityIterator& ent2  )  {

     // Extract pointer to reference element and get coordinates
     Vector<Double> CoordAtIP;
     Matrix<Double> tempElemMat;

     // set matrix to desired size and set all elements to zero
    //    partElemMat.Resize(numFncs);


     //New additions
     ExtractElemInfo( ent1 );
     ptelem->SetAnsatzFct( ansatzFct1_ );
     UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );

     //check for spectral element approximation if not, throw error
     assert(ansatzFct1_->GetType() == AnsatzFct::SPECTRAL );

     const UInt nrIntPts= ptelem->GetNumIntPoints();    
     const UInt spaceDim = ptelem->GetDim();  
     const Vector<Double> & intWeights = ptelem->GetIntWeights();
     Double jacDet;

     elemMat.Resize(spaceDim*numFncs);
     elemMat.Init();

     Matrix<Double>  partElemMat;
     partElemMat.Resize(spaceDim*numFncs);
     partElemMat.Init();

     Vector<Double> shapeFncAtIp;
     Matrix<Double> JacMat;
     Matrix<Double> JacMatT;
     Matrix<Double> subMat; 
     JacMat.Resize(spaceDim,spaceDim);
     JacMatT.Resize(spaceDim,spaceDim);
     subMat.Resize(spaceDim,spaceDim);
     JacMat.Init();
     JacMatT.Init();
     subMat.Init();


     for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
       
       jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem() );       
       ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );
       ptelem->CalcJacobianAtIp(JacMat,actIntPt,ptCoord_,ent1.GetElem() );

       JacMat.Transpose(JacMatT);
       //subMat = JacMat;
       JacMatT.Mult(JacMat,subMat);
       //subMat *= 1.0/jacDet;
       //subMat *= factor_;

       for( UInt intPti = 0 ; intPti< shapeFncAtIp.GetSize();intPti++){
         for( UInt intPtj = 0 ; intPtj< shapeFncAtIp.GetSize();intPtj++){
           for ( UInt i = 0;i<spaceDim ; i++ ) {
             for ( UInt j=0; j<spaceDim;j++ ) {
               partElemMat[(intPti*spaceDim)+i][(intPtj*spaceDim)+j] =  intWeights[actIntPt-1] * 
                                                                        subMat[i][j]  * 
                                                                        (1.0 / jacDet) *
                                                                        factor_ *
                                                                        shapeFncAtIp[intPti] * shapeFncAtIp[intPtj];
             }
           }
         }
       }

       elemMat += partElemMat;
     }
     //std::cout << "MassPiolaMixedInt_VV Matrix:\n" << elemMat << std::endl;
   }



 //=========================== Stiff-KPV==========================//
  StiffMixedInt_KPV::StiffMixedInt_KPV(Double aVal, UInt dim, bool axi, bool coordUpdate )
    : BaseForm(NULL), factor_ (aVal)  {
    
    name_ = "StiffMixed_KPV_Int";
    isaxi_ = axi;
    dim_   = dim;
    coordUpdate_ = coordUpdate;
  }


 
  StiffMixedInt_KPV::~StiffMixedInt_KPV()
  {
  }


  void StiffMixedInt_KPV::CalcElementMatrix( Matrix<Double>& elemMat,
                                             EntityIterator& ent1,
                                             EntityIterator& ent2 ) {
    //Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    
    // ansatzFct1_ corresponds to functions of pressure P (Taylor Hood: 1st order)
    // ansatzFct2_ corresponds to functions of veloctity V (TaylorHood: 2nd order )
    UInt numFncs1 = ptelem->GetNumFncs( ansatzFct1_ );
    UInt numFncs2 = ptelem->GetNumFncs( ansatzFct2_ );

    ptelem->SetAnsatzFct( ansatzFct2_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    Double jacDet;

    Vector<Double> shapeFncAtIp, ShpFncAtIp;
    Vector<Double> CoordAtIP;
    Matrix<Double> xiDx;
    Matrix<Double> K_x, K_y, K_z;

    //set matrix to desired size and set all elements to zero
    xiDx.Init();
    shapeFncAtIp.Init();
    K_x.Resize(numFncs1, numFncs2);
    K_y.Resize(numFncs1, numFncs2);
    K_x.Init();
    K_y.Init();
    if ( dim_ ==3 ) {
      K_z.Resize(numFncs1, numFncs2);
      K_z.Init();
    }

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      //get derivatives
      ptelem->SetAnsatzFct( ansatzFct1_ , false);
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet, ent1.GetElem());

      if (isaxi_) {
        ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                   ptCoord_, it1_.GetElem() );
        jacDet *= 2 * PI * CoordAtIP[0];
      }

      //Get the shape functions vector
      ptelem->SetAnsatzFct( ansatzFct2_ , false);
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );

      Double val = jacDet * intWeights[actIntPt-1] * factor_;
      for(UInt i=0; i<numFncs1; i++ ) {
        for(UInt j=0; j<numFncs2; j++ ) {
          K_x[i][j] += xiDx[i][0]*shapeFncAtIp[j] * val;
          K_y[i][j] += xiDx[i][1]*shapeFncAtIp[j] * val;
	        if ( dim_ == 3 ) {
	          K_z[i][j] += xiDx[i][2]*shapeFncAtIp[j] * val;
          }
        }
      }
    }

    elemMat.Resize(numFncs1, numFncs2*dim_);
    elemMat.Init();

    //Here the solution vector is P-vector, the element matrix is of dimension (numFncs1 X dim*numFncs2) 
    for (UInt i=0; i<numFncs1; i++) {
      for (UInt j=0; j<numFncs2; j++) {
        elemMat[i][j*dim_] = -K_x[i][j] ;
        elemMat[i][j*dim_+1] = -K_y[i][j];
	if ( dim_ == 3 )
	  elemMat[i][j*dim_+2] = -K_z[i][j];
      }
    }

    //std::cout << "ElemeMatStiff KPV:\n" << elemMat << std::endl;

}



//=========================== Stiff-Piola-KPV==========================//
  StiffPiolaMixedInt_KPV::StiffPiolaMixedInt_KPV(Double aVal, bool axi, bool coordUpdate )
   : BaseForm(NULL), factor_ (aVal)  {
    
    name_ = "StiffMixed_KPV_Int with PIOLA";
    isaxi_ = axi;
    coordUpdate_ = coordUpdate;
  }


 
  StiffPiolaMixedInt_KPV::~StiffPiolaMixedInt_KPV()
  {
  }


  void StiffPiolaMixedInt_KPV::CalcElementMatrix( Matrix<Double>& elemMat,
                                             EntityIterator& ent1,
                                             EntityIterator& ent2 ) {
    //Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

     //check for spectral element approximation if not, throw error
     assert(ansatzFct1_->GetType() == AnsatzFct::SPECTRAL );
    
    // ansatzFct1_ corresponds to functions of pressure P (Taylor Hood: 1st order)
    // ansatzFct2_ corresponds to functions of veloctity V (TaylorHood: 2nd order )
    UInt numFncs1 = ptelem->GetNumFncs( ansatzFct1_ );
    UInt numFncs2 = ptelem->GetNumFncs( ansatzFct2_ );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const UInt spaceDim = ptelem->GetDim();  
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Vector<Double> CoordAtIP;
    Matrix<Double> xiDx;

    //set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs1, numFncs2*spaceDim);
    elemMat.Init();
    xiDx.Init();
    shapeFncAtIp.Init();

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      //get derivatives
      ptelem->SetAnsatzFct( ansatzFct1_ , false);
      ptelem->GetLocDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet, ent1.GetElem());

      //Get the shape functions vector
      ptelem->SetAnsatzFct( ansatzFct2_ , false);
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );

      Double val = intWeights[actIntPt-1] * factor_;
      for(UInt i=0; i<numFncs1; i++ ) {
        for(UInt j=0; j<numFncs2; j++ ) {
          for(UInt d = 0; d<spaceDim;d++) {
            elemMat[i][(j*spaceDim)+d] -= shapeFncAtIp[j] * val * xiDx[i][d];
          }
        }
      }
    }
    //std::cout << "ElemeMatStiff KPV PIOLA:\n" << elemMat << std::endl;

}

  //=========================== Stiff-KVP==========================//
  StiffMixedInt_KVP::StiffMixedInt_KVP(Double aVal, UInt dim, bool axi, bool coordUpdate )
    : BaseForm(NULL), factor_ (aVal)
  {


    name_ = "StiffMixed_KVP_Int";
    isaxi_ = axi;
    dim_   = dim;
    coordUpdate_ = coordUpdate;
  }


 
  StiffMixedInt_KVP::~StiffMixedInt_KVP()
  {


  }


  void StiffMixedInt_KVP::CalcElementMatrix( Matrix<Double>& elemMat,
                                             EntityIterator& ent1,
                                             EntityIterator& ent2 ) {

    
    
    //Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    // ansatzFct1_ corresponds to functions of velocity V (Taylor Hood: 2nd order)
    // ansatzFct2_ corresponds to functions of pressure P (TaylorHood: 1st order )
    UInt numFncs1 = ptelem->GetNumFncs( ansatzFct1_ );
    UInt numFncs2 = ptelem->GetNumFncs( ansatzFct2_ );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    Double jacDet;

    Vector<Double> shapeFncAtIp, ShpFncAtIp;
    Vector<Double> CoordAtIP;
    Matrix<Double> xiDx;
    Matrix<Double> K_x, K_y, K_z;

    //set matrix to desired size and set all elements to zero
    xiDx.Init();
    shapeFncAtIp.Init();
    K_x.Resize(numFncs1, numFncs2);
    K_y.Resize(numFncs1, numFncs2);
    K_x.Init();
    K_y.Init();
    if ( dim_ == 3 ) {
      K_z.Resize(numFncs1, numFncs2);
      K_z.Init();
    }


    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      ptelem->SetAnsatzFct( ansatzFct2_ , false);
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet, ent1.GetElem());

      if (isaxi_) {
        ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                   ptCoord_, it1_.GetElem() );
        jacDet *= 2 * PI * CoordAtIP[0];
      }

       //Get the shape functions vector
      ptelem->SetAnsatzFct( ansatzFct1_ , false);
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );

      Double val = jacDet * intWeights[actIntPt-1] * factor_;
      for(UInt i = 0; i < numFncs1; i++ ) {
        for(UInt j = 0; j < numFncs2; j++ ) {
	  K_x[i][j] += xiDx[j][0]*shapeFncAtIp[i] * val;
          K_y[i][j] += xiDx[j][1]*shapeFncAtIp[i] * val;
	  if ( dim_ == 3 ) 
	    K_z[i][j] += xiDx[j][2]*shapeFncAtIp[i] * val;
        }
      }
    }

    elemMat.Resize(numFncs1*dim_, numFncs2);
    elemMat.Init();

    //Here the solution vector is (U,V)-vector, the element matrix is of dimension (numFncs1*2 X numFncs2) 
    for (UInt i = 0; i < numFncs1; i++) {
      for (UInt j = 0; j < numFncs2; j++) {
        elemMat[i*dim_][j] = K_x[i][j] ;
        elemMat[i*dim_+1][j] = K_y[i][j] ;
	if ( dim_ == 3 ) 
	  elemMat[i*dim_+2][j] = K_z[i][j] ;
      }
    }
    //    std::cout << "ElemeMatStiff KVP:\n" << elemMat << std::endl;
}


//=========================== Stiff-Piola-KVP==========================//
  StiffPiolaMixedInt_KVP::StiffPiolaMixedInt_KVP(Double aVal, bool axi, bool coordUpdate )
    : BaseForm(NULL), factor_ (aVal)
  {
    name_ = "StiffMixed_KVP_Int with PIOLA";
    isaxi_ = axi;
    coordUpdate_ = coordUpdate;
  }


 
  StiffPiolaMixedInt_KVP::~StiffPiolaMixedInt_KVP()
  {


  }


  void StiffPiolaMixedInt_KVP::CalcElementMatrix( Matrix<Double>& elemMat,
                                             EntityIterator& ent1,
                                             EntityIterator& ent2 ) {

    
    
    //Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

     //check for spectral element approximation if not, throw error
     assert(ansatzFct1_->GetType() == AnsatzFct::SPECTRAL );

    // ansatzFct1_ corresponds to functions of velocity V (Taylor Hood: 2nd order)
    // ansatzFct2_ corresponds to functions of pressure P (TaylorHood: 1st order )
    UInt numFncs1 = ptelem->GetNumFncs( ansatzFct1_ );
    UInt numFncs2 = ptelem->GetNumFncs( ansatzFct2_ );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const UInt spaceDim = ptelem->GetDim();  
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Vector<Double> CoordAtIP;
    Matrix<Double> xiDx;
    
    //set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs1*spaceDim, numFncs2);
    elemMat.Init();
    xiDx.Init();
    shapeFncAtIp.Init();

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      //get derivatives
      ptelem->SetAnsatzFct( ansatzFct1_ , false);
      ptelem->GetLocDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet, ent1.GetElem());

      //Get the shape functions vector
      ptelem->SetAnsatzFct( ansatzFct2_ , false);
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );

      //std::cout << "LocalFnc intPT# " << actIntPt << ":\n" << shapeFncAtIp << std::endl;

      Double val = intWeights[actIntPt-1] * factor_;
      for(UInt i=0; i<numFncs1; i++ ) {
        for(UInt j=0; j<numFncs2; j++ ) {
          for(UInt d = 0; d<spaceDim;d++) {
            elemMat[(j*spaceDim)+d][i] += shapeFncAtIp[j] * val * xiDx[i][d];
          }
        }
      }
    }
    //std::cout << "ElemeMatStiff KVP PIOLA:\n" << elemMat << std::endl;
}


//=========================== Absorbing Boundaries 1st order==========================//


  ABC_MixedInt:: ABC_MixedInt(Double aVal, bool axi, bool coordUpdate )
    : BaseForm(NULL), factor_ (aVal) {
    
    name_ = "ABC_MixedInt";
    isaxi_ = axi;
    coordUpdate_ = coordUpdate;
  }


 
  ABC_MixedInt::~ABC_MixedInt()
  {

  }


  void ABC_MixedInt::CalcElementMatrix( Matrix<Double>& elemMat,
          EntityIterator& ent1,
          EntityIterator& ent2 ) {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    Double jacDet;
    Vector<Double> shapeFncAtIp, CoordAtIP;
    Matrix<Double> partElemMat;

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    const Vector<Double> * intPoints = ptelem->GetIntPoints();

    elemMat.Resize(numFncs);
    elemMat.Init();

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, 
                                           ent1.GetElem() );
        
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );
        
      partElemMat.DyadicMult(shapeFncAtIp);
        
      if ( isaxi_ ) {
        ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                   ptCoord_, it1_.GetElem() );
        partElemMat *= 2 * PI * intWeights[actIntPt-1] 
          * factor_ * jacDet * CoordAtIP[0];
      }
      else 
        partElemMat *= intWeights[actIntPt-1] * factor_ * jacDet;
        
      elemMat += partElemMat;
    }

    //    std::cout << "ABC Mat:\n" << elemMat << std::endl;
  }


  //=========================== surface bilinear form: normal particle velocity  ==============//


  SurfVel_MixedInt:: SurfVel_MixedInt(Double aVal, UInt dim, bool axi, bool coordUpdate ) 
  {    
    name_ = "SurfVel_MixedInt";
    isaxi_ = axi;
    dim_   = dim;
    coordUpdate_ = coordUpdate;
    factor_ = aVal;
  }


 
  SurfVel_MixedInt::~SurfVel_MixedInt()
  {

  }


  void SurfVel_MixedInt::CalcElementMatrix( Matrix<Double>& elemMat,
              EntityIterator& ent1,
              EntityIterator& ent2 ) {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    Double jacDet;
    Vector<Double> shapeFncAtIp1, shapeFncAtIp2, CoordAtIP;
    Matrix<Double> partMat, helpElemMat;

    // ansatzFct1_ corresponds to functions of pressure P (Taylor Hood: 1st order)
    // ansatzFct2_ corresponds to functions of veloctity V (TaylorHood: 2nd order )
    UInt numFncs1 = ptelem->GetNumFncs( ansatzFct1_ );
    UInt numFncs2 = ptelem->GetNumFncs( ansatzFct2_ );


    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    const Vector<Double> * intPoints = ptelem->GetIntPoints();

    elemMat.Resize(numFncs1, dim_*numFncs2);
    elemMat.Init();
    helpElemMat.Resize(numFncs1, numFncs2);
    helpElemMat.Init();
    partMat.Resize(numFncs1, numFncs2);

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      ptelem->SetAnsatzFct( ansatzFct1_ );
      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, 
                                           ent1.GetElem() );
      ptelem->GetShFncAtIp(shapeFncAtIp1, actIntPt, ent1.GetElem() );

      ptelem->SetAnsatzFct( ansatzFct2_ );
      ptelem->GetShFncAtIp(shapeFncAtIp2, actIntPt, ent1.GetElem() );
        
      for ( UInt i=0; i<numFncs1; i++ ){
        for ( UInt j=0; j<numFncs2; j++ ){
          partMat[i][j] = shapeFncAtIp1[i] * shapeFncAtIp2[j];
        }
      }
      if (isaxi_) {
        ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                  ptCoord_, it1_.GetElem() );
        partMat *= 2 * PI * intWeights[actIntPt-1] 
                    * factor_ * jacDet * CoordAtIP[0];
      }
      else { 
        partMat *= intWeights[actIntPt-1] * factor_ * jacDet;
      }

      helpElemMat += partMat;
    }

    // Create a multi-dof matrix, multiplied by normal vector
    for ( UInt iRow = 0; iRow<numFncs1; iRow++ ) {
      for ( UInt iCol = 0; iCol<numFncs2; iCol++ ) {
        for ( UInt iDof = 0; iDof<dim_; iDof++ ) {
          elemMat[iRow][iCol*dim_+iDof] = 
            normal_[iDof] * helpElemMat[iRow][iCol];
        }
      }
    }
   
    //    std::cout << "Surf Mat:\n" << elemMat << std::endl;
  }


  //================== Neumann integrator for surface velocity =====================!

  LinSurfVelocity::LinSurfVelocity( std::string amplitudeStr,
            std::string phaseStr,
            Double materialParam,
            bool isaxi ) 
    : LinearSurfForm() {

    name_ = "LinSurfVelocity";
    isaxi_ = isaxi;

    // store value and phase string
    amplitude_ = amplitudeStr;
    phase_ = phaseStr;
    matFactor_ = materialParam;
  }


  LinSurfVelocity::~LinSurfVelocity() {

  }

  void LinSurfVelocity::PrepareElemVector( Vector<Double> & elemVec,
             EntityIterator& ent) {
    
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const Vector<Double> * intPoints = ptelem->GetIntPoints();

    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> shapeFnc, CoordAtIP;


    // Calculate element vector  
    elemVec.Resize(numFncs);
    elemVec.Init(0.0);
    Double value = 0.0;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      
      ptelem->GetShFncAtIp(shapeFnc,actIntPt,ent.GetElem() );
      value = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent.GetElem()) * 
        intWeights[actIntPt-1] * matFactor_;
      
      if (isaxi_) {
        ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                   ptCoord_, ent.GetElem() );
        value *=  2 * PI * CoordAtIP[0];
      }
      
      shapeFnc *= value;
      elemVec += shapeFnc;
    }
  }



  void LinSurfVelocity::CalcElemVector( Vector<Double> & elemVec,
          EntityIterator& ent) {
    
    // compute element vector
    PrepareElemVector( elemVec, ent );

    // evaluate value for current element
    mParser_->SetExpr( mHandle_, amplitude_ );
    Double factor = mParser_->Eval( mHandle_ );

    // multiply element vector with factor
    elemVec *= factor;

  }

  void LinSurfVelocity::CalcElemVector( Vector<Complex> & elemVec,
          EntityIterator& ent) {

    // compute element vector
    Vector<Double> helpVec;
    PrepareElemVector( helpVec, ent );

    // evaluate value and phase for current element
    mParser_->SetExpr( mHandle_, amplitude_ );
    Double value = mParser_->Eval( mHandle_ );

    mParser_->SetExpr( mHandle_, phase_ );
    Double phase = mParser_->Eval( mHandle_ );

    // Note: Since phase is in (grad), we have to transform it into
    //        rad-value
    Double realPart = value * cos( phase / 180 * PI );
    Double imagPart = value * sin( phase / 180 * PI );

    // multiply element vector with complex factor
    Complex factor(realPart, imagPart);
    elemVec = helpVec * factor;
  }




 //=========================== Convective-KPP==========================//
  ConvectiveMixedInt_KPP::ConvectiveMixedInt_KPP( Vector<Double> valVec, UInt dim, 
						  bool axi, bool coordUpdate )
    : BaseForm(NULL)  {
    
    name_ = "ConvectiveMixed_KPV_Int";
    isaxi_ = axi;
    dim_   = dim;
    coordUpdate_ = coordUpdate;
    velVec_ = valVec;
  }


 
  ConvectiveMixedInt_KPP::~ConvectiveMixedInt_KPP()
  {
  }


  void ConvectiveMixedInt_KPP::CalcElementMatrix( Matrix<Double>& elemMat,
						  EntityIterator& ent1,
						  EntityIterator& ent2 ) {
    //Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    
    // ansatzFct1_ corresponds to functions of pressure P (Taylor Hood: 1st order)
    UInt numFncs1 = ptelem->GetNumFncs( ansatzFct1_ );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Vector<Double> CoordAtIP;
    Matrix<Double> xiDx;

    //set matrix to desired size and set all elements to zero
    xiDx.Init();
    shapeFncAtIp.Init();
    elemMat.Resize(numFncs1);
    elemMat.Init();

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      //get derivatives
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet, ent1.GetElem());

      if (isaxi_) {
        ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                   ptCoord_, it1_.GetElem() );
        jacDet *= 2 * PI * CoordAtIP[0];
      }

      //Get the shape functions vector
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );

      Double val = jacDet * intWeights[actIntPt-1];
      for(UInt i=0; i<numFncs1; i++ ) {
        for(UInt j=0; j<numFncs1; j++ ) {
	  Double velPart;
	  velPart =  xiDx[j][0] * velVec_[0] +  xiDx[j][1] * velVec_[1];
	  if ( dim_==3 ) 
	    velPart += xiDx[j][2] * velVec_[2];

	  elemMat[i][j] += shapeFncAtIp[i] * velPart * val;
	}
      }
    }

    //   std::cout << "ConvectiveElemMat KPV:\n" << elemMat << std::endl;

}


 //=========================== Convective-KVP==========================//
  ConvectiveMixedInt_KVV::ConvectiveMixedInt_KVV( Vector<Double> valVec, UInt dim, 
						  bool axi, bool coordUpdate )
    : BaseForm(NULL)  {
    
    name_ = "ConvectiveMixed_KVP_Int";
    isaxi_ = axi;
    dim_   = dim;
    coordUpdate_ = coordUpdate;
    velVec_ = valVec;
  }


 
  ConvectiveMixedInt_KVV::~ConvectiveMixedInt_KVV()
  {
  }


  void ConvectiveMixedInt_KVV::CalcElementMatrix( Matrix<Double>& elemMat,
						  EntityIterator& ent1,
						  EntityIterator& ent2 ) {
    //Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    
    // ansatzFct1_ corresponds to functions of velocity V (Taylor Hood: 2nd order)
    UInt numFncs1 = ptelem->GetNumFncs( ansatzFct1_ );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Vector<Double> CoordAtIP;
    Matrix<Double> xiDx;

    //set matrix to desired size and set all elements to zero
    xiDx.Init();
    shapeFncAtIp.Init();
    Matrix<Double>  partElemMat;
    partElemMat.Resize(numFncs1);
    partElemMat.Init();

     for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      //get derivatives
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet, ent1.GetElem());

      if (isaxi_) {
        ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                   ptCoord_, it1_.GetElem() );
        jacDet *= 2 * PI * CoordAtIP[0];
      }

      //Get the shape functions vector
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );
      
      Double val = jacDet * intWeights[actIntPt-1];
      for(UInt i=0; i<numFncs1; i++ ) {
        for(UInt j=0; j<numFncs1; j++ ) {
	  Double velPart;
	  velPart =  xiDx[j][0] * velVec_[0] +  xiDx[j][1] * velVec_[1];
	  if ( dim_==3 ) 
	    velPart += xiDx[j][2] * velVec_[2];

	  partElemMat[i][j] += val *  shapeFncAtIp[i] * velPart; 
	}
      }
     }

     elemMat.Resize(numFncs1*dim_);
     elemMat.Init();

     //Blowing up the element matrix to number of velocity components 
     for (UInt i=0; i < numFncs1; i++)
       for (UInt j=0; j < numFncs1; j++)
	 for (UInt actDof = 0; actDof < dim_ ; actDof++)
	   elemMat[i*dim_ + actDof][j*dim_ + actDof] = partElemMat[i][j]; 
     

     //    std::cout << "ConvectiveElemMat KVP:\n" << elemMat << std::endl;

}


} // end namespace CoupledField

