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
	 CoordAtIP = ptCoord_ * shapeFncAtIp;
	 partElemMat *= 2 * PI * intWeights[actIntPt-1] * factor_
	   * jacDet * CoordAtIP[0];
       }
       else
	 partElemMat *= intWeights[actIntPt-1] * factor_ * jacDet;
       
       elemMat += partElemMat;
     }
     
     //std::cout << "Part Element Mass Matrix, MPP:\n" << elemMat << std::endl;
  }




//=========================== Mass-MVV ==========================//

  MassMixedInt_VV::MassMixedInt_VV(const Double aFactor, bool axi, bool coordUpdate )
    : BaseForm(NULL), 
      factor_(aFactor)
  {
    name_ = "MassMixedVelocityInt";

    isaxi_ = axi;
    coordUpdate_ = coordUpdate;
    baseType_ = MASS;
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

     // set matrix to desired size and set all elements to zero
    //    partElemMat.Resize(numFncs);


     //New additions
     ExtractElemInfo( ent1 );
     ptelem->SetAnsatzFct( ansatzFct1_ );
     UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
     const UInt nrIntPts= ptelem->GetNumIntPoints();    
     const Vector<Double> & intWeights = ptelem->GetIntWeights();
     Double jacDet;

     elemMat.Resize(2*numFncs);
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
	 CoordAtIP = ptCoord_ * shapeFncAtIp;
	 partElemMat *= 2 * PI * intWeights[actIntPt-1] * factor_
	   * jacDet * CoordAtIP[0];
       }
       else
	 partElemMat *= intWeights[actIntPt-1] * factor_ * jacDet;
       
       tempElemMat += partElemMat;
     }
     


     const UInt singleDofSize = tempElemMat.GetSizeRow();
     UInt i, j, actDof;
     //    std::cout << "Part Element Mass Matrix, MPP:\n" << tempElemMat << std::endl;
     //Blowing up the element matrix to number of velocity components 
     for (i=0; i < singleDofSize; i++)
       for ( j=0; j < singleDofSize; j++)
	 for (actDof=0; actDof < 2 ; actDof++)
	   elemMat[i*2 + actDof][j*2 + actDof] = tempElemMat[i][j]; 
     
     //std::cout << "ElemMatMass:\n" << elemMat << std::endl;

   }



  //=========================== Stiff-KPV==========================//
  StiffMixedInt_KPV::StiffMixedInt_KPV(Double aVal, bool axi, bool coordUpdate )
    : BaseForm(NULL), factor_ (aVal)
  {


    name_ = "StiffMixed_KPV_Int";
    isaxi_ = axi;
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

    // ansatzFct1_ corresponds to functions of pressure V (Taylor Hood: 2nd order)
    // ansatzFct2_ corresponds to functions of veloctity P (TaylorHood: 1st order )
    UInt numFncs1 = ptelem->GetNumFncs( ansatzFct1_ );
    UInt numFncs2 = ptelem->GetNumFncs( ansatzFct2_ );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Vector<Double> CoordAtIP;
    Matrix<Double> xiDx;
    Matrix<Double> xiDxTransp;

    UInt  dim = ptelem->GetDim()+1;
    int i, j, actDof;

    //set matrix to desired size and set all elements to zero
    xiDx.Init();
    xiDxTransp.Init();
    shapeFncAtIp.Init();
    Matrix<Double> K_x, K_y, K_x_Temp, K_y_Temp;
    K_x.Resize(numFncs1, numFncs2);
    K_y.Resize(numFncs1, numFncs2);
    K_x_Temp.Resize(numFncs1, numFncs2);
    K_y_Temp.Resize(numFncs1, numFncs2);
    K_x.Init();
    K_y.Init();
    K_x_Temp.Init();
    K_y_Temp.Init();



    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) 
      {
	//	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem() );
	ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet, ent1.GetElem());

	//Get the traspose	matrix
        ptelem->SetAnsatzFct( ansatzFct2_ , false);
	xiDx.Transpose(xiDxTransp);
	//Get the shape functions vector

	ptelem->SetAnsatzFct( ansatzFct1_ , false);
        ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );
	for(i = 0; i < numFncs1; i++ )
	  for(j = 0; j < numFncs2; j++ )
	    {
	      // 						K_x[i][j] += xiDxTransp[0][j]*shapeFncAtIp[i];
	      // 						K_y[i][j] += xiDxTransp[1][j]*shapeFncAtIp[i];
	      K_x_Temp[i][j] = xiDxTransp[0][j]*shapeFncAtIp[i] * jacDet * intWeights[actIntPt-1] * factor_;
	      K_y_Temp[i][j] = xiDxTransp[1][j]*shapeFncAtIp[i] * jacDet * intWeights[actIntPt-1] * factor_;;
	    }
	K_x += K_x_Temp;
	K_y += K_y_Temp;
      }
		
//     std::cout << "K_x = "<< std::endl << K_x << std::endl;
//     std::cout << "K_y = "<< std::endl << K_y << std::endl;

    const UInt singleDofSize = K_x.GetSizeRow();
    elemMat.Resize(numFncs1*2, numFncs2);
    elemMat.Init();

    //Here the solution vector is (U,V)-vector, the element matrix is of dimension (numFncs X numFncs*2) 
    for (i = 0; i < numFncs1; i++) 
      for (j = 0; j < numFncs2; j++)
	  {
		elemMat[i*2][j] = K_x[i][j] ;
		elemMat[i*2+1][j] = K_y[i][j] ;
	  }

    //  std::cout << "ElemeMatStiff KPV:\n" << elemMat << std::endl;

}




  //=========================== Stiff-KVP==========================//
  StiffMixedInt_KVP::StiffMixedInt_KVP(Double aVal, bool axi, bool coordUpdate )
    : BaseForm(NULL), factor_ (aVal)
  {


    name_ = "StiffMixed_KVP_Int";
    isaxi_ = axi;
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
    
    // ansatzFct1_ corresponds to functions of pressure P (Taylor Hood: 1st order)
    // ansatzFct2_ corresponds to functions of veloctity V (TaylorHood: 2nd order )
    UInt numFncs1 = ptelem->GetNumFncs( ansatzFct1_ );
    UInt numFncs2 = ptelem->GetNumFncs( ansatzFct2_ );

    ptelem->SetAnsatzFct( ansatzFct2_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Vector<Double> CoordAtIP;
    Matrix<Double> xiDx;
    Matrix<Double> xiDxTransp;

    UInt  dim = ptelem->GetDim()+1;
    int i, j, actDof;

    //set matrix to desired size and set all elements to zero
    xiDx.Init();
    xiDxTransp.Init();
    shapeFncAtIp.Init();
    Matrix<Double> K_x, K_y, K_x_Temp, K_y_Temp;
    K_x.Resize(numFncs1, numFncs2);
    K_y.Resize(numFncs1, numFncs2);
    K_x_Temp.Resize(numFncs1, numFncs2);
    K_y_Temp.Resize(numFncs1, numFncs2);
    K_x.Init();
    K_y.Init();
    K_x_Temp.Init();
    K_y_Temp.Init();



    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) 
      {
	//	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem() );
	ptelem->SetAnsatzFct( ansatzFct2_ , false);
	ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet, ent1.GetElem());

	//Get the traspose	matrix
	xiDx.Transpose(xiDxTransp);
	//Get the shape functions vector

        ptelem->SetAnsatzFct( ansatzFct1_ , false);
	ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent1.GetElem() );
	for(i = 0; i < numFncs1; i++ )
	  for(j = 0; j < numFncs2; j++ )
	    {
	      // 						K_x[i][j] += xiDxTransp[0][j]*shapeFncAtIp[i];
	      // 						K_y[i][j] += xiDxTransp[1][j]*shapeFncAtIp[i];
	      K_x_Temp[i][j] = xiDxTransp[0][j]*shapeFncAtIp[i] * jacDet * intWeights[actIntPt-1] * factor_;
	      K_y_Temp[i][j] = xiDxTransp[1][j]*shapeFncAtIp[i] * jacDet * intWeights[actIntPt-1] * factor_;;
	    }
	K_x += K_x_Temp;
	K_y += K_y_Temp;
      }
		
//     std::cout << "K_x = "<< std::endl << K_x << std::endl;
//     std::cout << "K_y = "<< std::endl << K_y << std::endl;

    const UInt singleDofSize = K_x.GetSizeRow();
    elemMat.Resize(numFncs1, numFncs2*2);
    elemMat.Init();

    //Here the solution vector is P-vector, the element matrix is of dimension (numFncs1 X numFncs2) 
    for (i = 0; i < numFncs1; i++) 
      for (j = 0; j < numFncs2; j++)
	  {
		elemMat[i][j*2] = K_x[i][j] ;
		elemMat[i][j*2+1] = K_y[i][j] ;
	  }

    //  std::cout << "ElemeMatStiff KVP:\n" << elemMat << std::endl;
}



} // end namespace CoupledField









































