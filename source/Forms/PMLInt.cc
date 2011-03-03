// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include "Domain/domain.hh"
#include "Forms/pmlBasics.hh"

#include "PMLInt.hh"

namespace CoupledField
{

  PMLInt::PMLInt(std::string type, Double factor, std::string dampingTypePML, Double damp, 
		 bool axi)
    : BaseForm(NULL), formsFactor_(factor)
  {
    name_ = "PMLInt";

    isComplex_ = true;
    isaxi_     = axi;

    nrDofsPerNode_ = 1;
    pmlFnc_    = new PMLBasics( dampingTypePML, damp, type);

    // Set Expression for parser
    mParser_->SetExpr( mHandle_, "f" );
  }


 
  PMLInt::~PMLInt()
  {
  }



  void PMLInt::CalcElementMatrix( Matrix<Complex>& elemMat,
                                  EntityIterator& ent1, 
                                  EntityIterator& ent2 )
  {
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    //correct PML
    if ( pmlFnc_->GetFormsType() =="laplaceInt" ) {
      CalcElementMatrixStiff(ptCoord_, elemMat);
    }
    else if ( pmlFnc_->GetFormsType() =="massInt" ) {
      CalcElementMatrixMass(ptCoord_, elemMat);
    }
    //almost PML
    else if ( pmlFnc_->GetFormsType() =="laplaceIntAPML" ) {
      CalcElementMatrixStiff4APML(ptCoord_, elemMat);
    }
    else if ( pmlFnc_->GetFormsType() =="massIntAPML" ) {
      CalcElementMatrixMass4APML(ptCoord_, elemMat);
    } 
  }


  void PMLInt::CalcElementMatrixStiff(Matrix<Double> & ptCoord, Matrix<Complex> & elemMat)
  {

    
    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    Double jacDet;  

    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx, xiDxTransp;
    Matrix<Complex> xiDxC, xiDxTranspC;
    Matrix<Complex> partElemMat;
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs); 
    elemMat.Init();
   
    //set correct size for complex xiDx
    const UInt spaceDim = ptelem->GetDim();  
    xiDxC.Resize(numFncs, spaceDim);

    Complex jacDetC;
    Vector<Complex> factorsPML;
    Double omega = 2 * PI * mParser_->Eval( mHandle_ );

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, 
                                      jacDet, it1_.GetElem() );
 	// compute PML factor 
        ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                   ptCoord, it1_.GetElem() );
	pmlFnc_->ComputeFactorPML( factorsPML, jacDetC, CoordAtIP, omega);        

	    //multiply the derivatives with the x-,y- and z-factors
	    for (UInt i=0; i<xiDx.GetNumCols(); i++) {
	      for (UInt j=0; j<xiDx.GetNumRows(); j++) {
	        xiDxC[j][i] = xiDx[j][i] * factorsPML[i];
	      }
	    }

        xiDxC.Transpose(xiDxTranspC);
        partElemMat = xiDxC * xiDxTranspC;

        if (isaxi_) {
	  partElemMat *= 2 * PI * intWeights[actIntPt-1] * jacDet * 
            formsFactor_ * CoordAtIP[0] * jacDetC;
	}
        else 
          partElemMat *= intWeights[actIntPt-1] * jacDet * 
            formsFactor_ * jacDetC;

        elemMat += partElemMat;
      }

    //        std::cout << "PML-ElemMatStiff:\n" << elemMat << std::endl;
  }


  void PMLInt::CalcElementMatrixMass(Matrix<Double> & ptCoord, Matrix<Complex> & elemMat)
  {
    
    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs);
    elemMat.Init();
    
    Complex jacDetC;
    Vector<Complex> factorsPML;
    Complex factor, factorPML;
    Double omega = 2 * PI * mParser_->Eval( mHandle_ );

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      ptelem->SetAnsatzFct( ansatzFct1_ );
      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord,
                                           it1_.GetElem() );
        
      ptelem-> GetShFncAtIp(shapeFncAtIp, actIntPt, it1_.GetElem() );
        
      partElemMat.DyadicMult(shapeFncAtIp);
        
      // compute PML factor 
      ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                 ptCoord, it1_.GetElem() );

      //std::cout << "Mass CoordAtIP:\n" << CoordAtIP << std::endl;
      
      pmlFnc_->ComputeFactorPML( factorsPML, jacDetC, CoordAtIP, omega); 

      factorPML = jacDetC * formsFactor_;

      if (isaxi_) 
        factor = 2 * PI * intWeights[actIntPt-1] * jacDet * 
          CoordAtIP[0] * factorPML;
      else 
        factor = intWeights[actIntPt-1] * jacDet * factorPML;
        
      for ( UInt i=0; i<numFncs; i++ ) {
        for ( UInt j=0; j<numFncs; j++ ) {
          elemMat[i][j] += factor * partElemMat[i][j];
        }
      }
    }
    
    if (nrDofsPerNode_ > 1 ) {
     Matrix <Complex> singleDofMass = elemMat;
     UInt singleDofSize = singleDofMass.GetNumRows();

     elemMat.Resize( nrDofsPerNode_* singleDofSize );

     for (UInt i=0; i < singleDofSize; i++)
       for (UInt j=0; j < singleDofSize; j++)
         for (UInt actDof=0; actDof < nrDofsPerNode_; actDof++)
           elemMat[i*nrDofsPerNode_ + actDof][j*nrDofsPerNode_ + actDof] = singleDofMass[i][j]; 
   }

  }


  void PMLInt::SetPosPML(Matrix<Double> & inner, Matrix<Double> & outer)
  {

    pmlFnc_-> SetPosPML( inner, outer );
  }


  //------------------------harmonic case for almost PML---------------------------------

  void PMLInt::CalcElementMatrixStiff4APML(Matrix<Double> & ptCoord, 
                                           Matrix<Complex> & elemMat) {
    
    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    Double jacDet;  

    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Complex> xiDxC,xiDxTranspC;
    Matrix<Complex> partElemMat;
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs); 
    elemMat.Init();
   
    //set correct size for complex xiDx
    const UInt spaceDim = ptelem->GetDim();  
    xiDxC.Resize(numFncs, spaceDim);

    Vector<Complex> factorsPML;
    Double omega = 2 * PI * mParser_->Eval( mHandle_ );

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, 
                                      jacDet, it1_.GetElem() );
 	// compute PML factor 
        ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                   ptCoord, it1_.GetElem() );
	pmlFnc_->ComputeFactorAPML( factorsPML, CoordAtIP, omega);        

	//multiply the derivatives with the x-,y- and z-factors
	for (UInt i=0; i<xiDx.GetNumCols(); i++) {
	  for (UInt j=0; j<xiDx.GetNumRows(); j++) {
	    xiDxC[j][i] = xiDx[j][i] * factorsPML[i];
	  }
	}

        xiDxC.Transpose(xiDxTranspC);
        partElemMat = xiDx * xiDxTranspC;

        if (isaxi_) {
	  partElemMat *= 2 * PI * intWeights[actIntPt-1] * jacDet * 
            formsFactor_ * CoordAtIP[0];
	}
        else 
          partElemMat *= intWeights[actIntPt-1] * jacDet * 
            formsFactor_;

        elemMat += partElemMat;
      }

    //std::cout << "PML-ElemMatStiff:\n" << elemMat << std::endl;
  }


  void PMLInt::CalcElementMatrixMass4APML(Matrix<Double> & ptCoord, 
                                          Matrix<Complex> & elemMat) {
    
    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs);
    elemMat.Init();
    
    Vector<Complex> factorsPML;
    Complex factor, factorPML;
    Double omega = 2 * PI * mParser_->Eval( mHandle_ );

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      ptelem->SetAnsatzFct( ansatzFct1_ );
      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord,
                                           it1_.GetElem() );
        
      ptelem-> GetShFncAtIp(shapeFncAtIp, actIntPt, it1_.GetElem() );
        
      partElemMat.DyadicMult(shapeFncAtIp);
        
      // compute PML factor 
      ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                 ptCoord, it1_.GetElem() );

      pmlFnc_->ComputeFactorAPML( factorsPML, CoordAtIP, omega); 

      factorPML = factorsPML[0] * formsFactor_;

      if (isaxi_) 
        factor = 2 * PI * intWeights[actIntPt-1] * jacDet * 
          CoordAtIP[0] * factorPML;
      else 
        factor = intWeights[actIntPt-1] * jacDet * factorPML;
        
      for ( UInt i=0; i<numFncs; i++ ) {
        for ( UInt j=0; j<numFncs; j++ ) {
          elemMat[i][j] += factor * partElemMat[i][j];
        }
      }
    }
    
    if (nrDofsPerNode_ > 1 ) {
     Matrix <Complex> singleDofMass = elemMat;
     UInt singleDofSize = singleDofMass.GetNumRows();

     elemMat.Resize( nrDofsPerNode_* singleDofSize );

     for (UInt i=0; i < singleDofSize; i++)
       for (UInt j=0; j < singleDofSize; j++)
         for (UInt actDof=0; actDof < nrDofsPerNode_; actDof++)
           elemMat[i*nrDofsPerNode_ + actDof][j*nrDofsPerNode_ + actDof] = singleDofMass[i][j]; 
   }


    // std::cout << "PML-ElemMatMass:\n" << elemMat << std::endl;
  }


  //========================================= MECHANICAL ===============================

  MechPMLInt::MechPMLInt( std::string type,  BaseMaterial* matData,  
                          std::string dampingTypePML, Double damp, 
                          SubTensorType tensorType)
    : linElastInt( matData, tensorType)
  {
    name_ = "MechPMLInt";

    isComplex_ = true;
    pmlFnc_    = new PMLBasics( dampingTypePML, damp, type);

    // Set Expression for parser
    mParser_->SetExpr( mHandle_, "f" );
  }


 
  MechPMLInt::~MechPMLInt()
  {
  }



  void MechPMLInt::CalcElementMatrix( Matrix<Complex>& elemMat,
                                      EntityIterator& ent1, 
                                      EntityIterator& ent2 ) {


    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );


    // First of all, set ansatz function to element
    ptelem->SetAnsatzFct( ansatzFct1_ );

    UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrDofs   = getNrDofs();  

    double jacDet;

    Matrix<Double> bMat; 
    Matrix<Double> dMat; 
    Matrix<Complex> bMatC, dbMatC, dMatC; 
    Vector<Double> CoordAtIP;
    Complex aux, *ptr1, *ptr2, jacDetC, fac;

    elemMat.Resize( nrFncs * nrDofs);
    elemMat.Init();

    //get integration points    
    const UInt nrIntPts = ptelem->GetNumIntPoints(); 
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    const Vector<Double> * intPoints = ptelem->GetIntPoints();

    // get material tensor
    calcDMat( dMat, ent1.GetElem());
    dMatC =  dMat * Complex(1.0,0);

    // Loop over all integration points
    for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {
      
      // Setup the B matrix for current integration point
      calcBMat( bMat, actIntPt, ptCoord_ );
      
      // compute PML factor 
      ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                 ptCoord_, it1_.GetElem() );
 
      calcBMatPML( bMat, CoordAtIP, bMatC, jacDetC );

      // Compute Jacobian for integration point
      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, ent1.GetElem() );
      
      // Perform a safety check
      if ( jacDet < 0.0 ) {
        EXCEPTION( "BDBInt::CalcElementMatrix: Encountered "
                 << "negative Jacobian determinant!" );
      }
      
      // Special things must be done in the axi-symmetric case
      if ( isaxi_ ) {
        Vector<Double> ShpFncAtIp;
        ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );
        
        CoordAtIP = ptCoord_ * ShpFncAtIp;
        jacDet *= 2 * PI * CoordAtIP[0];
      }
      
      // Compute the matrix product D * B and store as intermediate matrix
      //      dbMatC.Resize( dMat.GetNumRows(), bMat.GetNumCols() );
      dbMatC = dMatC * bMatC;
      
      // We now compute B^T * D * B and scale it by the determinant
      // of the Jacobian and the weight of the current integration
      // point. The result is added to the element matrix.
      fac = jacDet * intWeights[actIntPt-1] * jacDetC;
      for ( UInt k = 0; k < bMat.GetNumRows(); k++ ) {
        ptr1 =  bMatC[k];
        ptr2 = dbMatC[k];
        for ( UInt i = 0; i < bMatC.GetNumCols(); i++ ) {
          aux = fac * ptr1[i];
          for ( UInt j = 0; j < dbMatC.GetNumCols(); j++ ) {
            elemMat[i][j] += aux * ptr2[j];
          }
        }
      }
    }
  }
  

  void MechPMLInt::calcBMatPML( Matrix<Double>& bMat, Vector<Double>& CoordAtIP, 
                                Matrix<Complex>& bMatC, Complex& jacDetC)
  {

    bMatC.Resize( bMat.GetNumRows(),  bMat.GetNumCols() );
    bMatC.Init();

    Double omega = 2 * PI * mParser_->Eval( mHandle_ );
    Vector<Complex> factorsPML;

    pmlFnc_->ComputeFactorPML( factorsPML, jacDetC, CoordAtIP, omega);
    //    std::cout << "FactorsPML:\n" << factorsPML << std::endl;

    UInt idx;
    const UInt spaceDim = ptelem->GetDim();  
    const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );

    for( UInt actNode=0; actNode < nrFncs; actNode++) {
      if ( spaceDim < 3 ) {
        //2D problem; here plane and axi is the same!!!
        idx = actNode * spaceDim;
        bMatC[0][idx]   = bMat[0][idx]   * factorsPML[0]; // d/dx
        bMatC[1][idx+1] = bMat[1][idx+1] * factorsPML[1]; // d/dy
        bMatC[2][idx]   = bMat[2][idx]   * factorsPML[1]; // d/dy
        bMatC[2][idx+1] = bMat[2][idx+1] * factorsPML[0]; // d/dx
      }
      else {
        //3D problem
        idx = actNode * spaceDim;
        bMatC[0][idx]   = bMat[0][idx]   * factorsPML[0]; // d/dx
        bMatC[1][idx+1] = bMat[1][idx+1] * factorsPML[1]; // d/dy
        bMatC[2][idx+2] = bMat[2][idx+2] * factorsPML[2]; // d/dz
        bMatC[3][idx+1] = bMat[3][idx+1] * factorsPML[2]; // d/dz
        bMatC[3][idx+2] = bMat[3][idx+2] * factorsPML[1]; // d/dy
        bMatC[4][idx]   = bMat[4][idx]   * factorsPML[2]; // d/dz
        bMatC[4][idx+2] = bMat[4][idx+2] * factorsPML[0]; // d/dx
        bMatC[5][idx]   = bMat[5][idx]   * factorsPML[1]; // d/dy
        bMatC[5][idx+1] = bMat[5][idx+1] * factorsPML[0]; // d/dx
      }
    }
//     std::cout << "Bop:\n" << bMat << std::endl;
//     std::cout << "BopC:\n" << bMatC << std::endl;
  }


  void MechPMLInt::SetPosPML(Matrix<Double> & inner, Matrix<Double> & outer)
  {

    pmlFnc_-> SetPosPML( inner, outer );
  }



  //---------------------------Time domain PML -----------------------------------//


  PMLTimeInt::PMLTimeInt(std::string type, Double factor, std::string dampingTypePML, 
                         Double damp, bool axi)
    : BaseForm(NULL), formsFactor_(factor)
  {
    name_ = "PMLTimeInt";

    isaxi_     = axi;
    formsType_ = type;
    pmlFnc_    = new PMLBasics( dampingTypePML, damp, type);

  }


 
  PMLTimeInt::~PMLTimeInt()
  {
  }



  void PMLTimeInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                      EntityIterator& ent1, 
                                      EntityIterator& ent2 )
  {
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    if ( formsType_ == "pressureStiff" || formsType_ =="pressureDamp" 
         || formsType_ == "scalarAuxStiff" ) {
      CalcElementMatrixPressureOrAux(ptCoord_, elemMat);
    }
    else if ( formsType_ == "pressureGrad" || formsType_ == "auxGrad") {
      CalcElementMatrixPressureOrAuxGrad(ptCoord_, elemMat);
    }
    else if ( formsType_ == "vecAuxillaryDiv" ) {
      CalcElementMatrixVecAuxillaryDiv(ptCoord_, elemMat);
    }
    else if ( formsType_ == "vecAuxillaryStiff" ) {
      CalcElementMatrixVecAuxillaryStiff(ptCoord_, elemMat);
    }
  }


  void PMLTimeInt::CalcElementMatrixPressureOrAux(Matrix<Double>& ptCoord, 
                                                  Matrix<Double>& elemMat)
  {
    
    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    const UInt dim = ptCoord.GetNumRows();
    Double jacDet;  

    // some variables
    Matrix<Double> partElemMat;
    Vector<Double> shapeFncAtIp;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs); 
    elemMat.Init();
   

    Vector<Double> factorsPML;
    Double pmlVal;

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      
      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, it1_.GetElem() );    
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, it1_.GetElem() );
      partElemMat.DyadicMult(shapeFncAtIp);

      // compute PML factor 
      ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                 ptCoord, it1_.GetElem() );
      pmlFnc_->ComputeTimeFactorPML( factorsPML, CoordAtIP );  

      if ( formsType_ == "pressureStiff" ) {
        pmlVal = factorsPML[0] * factorsPML[1];
        if ( dim == 3 ) {
          pmlVal += factorsPML[0] * factorsPML[2];
          pmlVal += factorsPML[1] * factorsPML[2];
        }
      }
      else if ( formsType_ == "pressureDamp" ) {
        pmlVal = 0.0;
        for ( UInt k=0; k<dim; k++)
          pmlVal += factorsPML[k];
      }
      else {
        // stiff matrix for scalar auxillary, just in 3D!!
        pmlVal = 1.0;
        for ( UInt k=0; k<dim; k++)
          pmlVal *= factorsPML[k];
      }

      //std::cout << "Type: " << formsType_ << "  pmlVal: " << pmlVal << std::endl;

      if (isaxi_) {
        partElemMat *= 2 * PI * intWeights[actIntPt-1] * formsFactor_ 
                       * pmlVal * jacDet * CoordAtIP[0];
      }
      else 
        partElemMat *= intWeights[actIntPt-1] * formsFactor_ * pmlVal * jacDet;
      
      elemMat += partElemMat;
    }

    //std::cout << "Type: " << formsType_ << "\n" << elemMat << std::endl;

  }

  void PMLTimeInt::CalcElementMatrixPressureOrAuxGrad(Matrix<Double>& ptCoord, Matrix<Double>& elemMat)
  {

    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    const UInt dim = ptCoord.GetNumRows();

    Vector<Double> shapeFncAtIp, ShpFncAtIp;
    Vector<Double> CoordAtIP;
    Matrix<Double> xiDx;
    Matrix<Double> K_x, K_y, K_z;
    Double jacDet;

    //set matrix to desired size and set all elements to zero
    xiDx.Init();
    shapeFncAtIp.Init();
    K_x.Resize(numFncs);
    K_y.Resize(numFncs);
    K_x.Init();
    K_y.Init();
    if ( dim == 3 ) {
      K_z.Resize(numFncs);
      K_z.Init();
    }

    Vector<Double> factorsPML;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet, it1_.GetElem());

      // compute PML factor 
      ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                 ptCoord, it1_.GetElem() );
      pmlFnc_->ComputeTimeFactorPML( factorsPML, CoordAtIP );  

      if (isaxi_) {
        jacDet *= 2 * PI * CoordAtIP[0];
      }

      //Get the shape functions vector
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, it1_.GetElem() );
    
      Double val, pmlX, pmlY, pmlZ;
      val = jacDet * intWeights[actIntPt-1] * formsFactor_;
      pmlX = pmlY = pmlZ = 0.0;
      if ( formsType_ == "pressureGrad" ) {
        pmlX = factorsPML[0] - factorsPML[1];
        pmlY = factorsPML[1] - factorsPML[0];
        if ( dim == 3 ) {
          pmlX -= factorsPML[2];
          pmlY -= factorsPML[2];
          pmlZ = factorsPML[2] - factorsPML[0] - factorsPML[1];
        }
      }
      else {
        // aux grad, just in 3D
        pmlX = factorsPML[1] * factorsPML[2];
        pmlY = factorsPML[0] * factorsPML[2];
        pmlZ = factorsPML[0] * factorsPML[1];
      }

//       std::cout << "Type: " << formsType_ << "  pmlValX: " << pmlX << std::endl;
//       std::cout << "Type: " << formsType_ << "  pmlValY: " << pmlY << std::endl;
//       std::cout << "Type: " << formsType_ << "  pmlValZ: " << pmlZ << std::endl;

      for(UInt i = 0; i < numFncs; i++ ) {
        for(UInt j = 0; j < numFncs; j++ ) {
	  K_x[i][j] += xiDx[j][0]*shapeFncAtIp[i] * val * pmlX;
          K_y[i][j] += xiDx[j][1]*shapeFncAtIp[i] * val * pmlY;
	  if ( dim == 3 ) 
	    K_z[i][j] += xiDx[j][2]*shapeFncAtIp[i] * val * pmlZ;
        }
      }
    }

    elemMat.Resize(numFncs*dim, numFncs);
    elemMat.Init();

    // Here the solution vector is (U,V,W)-vector, the element matrix is of 
    // dimension (numFncs*dim X numFncs) 
    for (UInt i = 0; i < numFncs; i++) {
      for (UInt j = 0; j < numFncs; j++) {
        elemMat[i*dim][j] = K_x[i][j] ;
        elemMat[i*dim+1][j] = K_y[i][j] ;
	if ( dim == 3 ) 
	  elemMat[i*dim+2][j] = K_z[i][j] ;
      }
    }

    //    std::cout << "Type: " << formsType_ << "\n" << elemMat << std::endl;

  }

  void PMLTimeInt::CalcElementMatrixVecAuxillaryStiff(Matrix<Double>& ptCoord, Matrix<Double>& elemMat)
  {
    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const UInt dim = ptCoord.GetNumRows();

    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat, elemMatX, elemMatY, elemMatZ;
    Vector<Double> CoordAtIP;
    Double jacDet;

    // set matrix to desired size and set all elements to zero
    //    partElemMat.Resize(numFncs);
    elemMatX.Resize(numFncs);
    elemMatX.Init();
    elemMatY.Resize(numFncs);
    elemMatY.Init();
    if ( dim == 3 ) {
      elemMatZ.Resize(numFncs);
      elemMatZ.Init();
    }

    Vector<Double> factorsPML(dim);
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      
      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, it1_.GetElem() );
      
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, it1_.GetElem() );
      
      partElemMat.DyadicMult(shapeFncAtIp);
      
      // compute PML factor 
      ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                 ptCoord, it1_.GetElem() );
      pmlFnc_->ComputeTimeFactorPML( factorsPML, CoordAtIP );  

//       std::cout << "Type: " << formsType_ << "  pmlValX: " << factorsPML[0] << std::endl;
//       std::cout << "Type: " << formsType_ << "  pmlValY: " << factorsPML[1] << std::endl;
//       std::cout << "Type: " << formsType_ << "  pmlValZ: " << factorsPML[2] << std::endl;

      if (isaxi_) {
        partElemMat *= 2 * PI * intWeights[actIntPt-1] * formsFactor_ * jacDet * CoordAtIP[0];
      }
      else 
        partElemMat *= intWeights[actIntPt-1] * formsFactor_ * jacDet;
      
      elemMatX += partElemMat  * factorsPML[0];
      elemMatY += partElemMat  * factorsPML[1];
      if ( dim==3 )
        elemMatZ += partElemMat  * factorsPML[2];
    }

    elemMat.Resize(numFncs*dim);
    elemMat.Init();

    for ( UInt i=0; i < numFncs; i++ )
      for ( UInt j=0; j < numFncs; j++ ) {
        elemMat[i*dim][j*dim]         = elemMatX[i][j]; 
        elemMat[i*dim + 1][j*dim + 1] = elemMatY[i][j]; 
        if ( dim == 3 ) 
          elemMat[i*dim + 2][j*dim + 2] = elemMatZ[i][j]; 
      }

    //std::cout << "Type: " << formsType_ << "\n" << elemMat << std::endl;

  }


  void PMLTimeInt::CalcElementMatrixVecAuxillaryDiv(Matrix<Double>& ptCoord, Matrix<Double>& elemMat)
  {
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    const UInt dim = ptCoord.GetNumRows();

    Vector<Double> shapeFncAtIp, ShpFncAtIp;
    Vector<Double> CoordAtIP;
    Matrix<Double> xiDx;
    Matrix<Double> K_x, K_y, K_z;
    Double jacDet;

    //set matrix to desired size and set all elements to zero
    xiDx.Init();
    shapeFncAtIp.Init();
    K_x.Resize(numFncs);
    K_y.Resize(numFncs);
    K_x.Init();
    K_y.Init();
    if ( dim ==3 ) {
      K_z.Resize(numFncs);
      K_z.Init();
    }

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      //get derivatives
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet, it1_.GetElem());

      if (isaxi_) {
        ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                   ptCoord, it1_.GetElem() );
        jacDet *= 2 * PI * CoordAtIP[0];
      }

      //Get the shape functions vector
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, it1_.GetElem() );

      Double val = jacDet * intWeights[actIntPt-1] * formsFactor_;
      for(UInt i=0; i<numFncs; i++ ) {
        for(UInt j=0; j<numFncs; j++ ) {
          K_x[i][j] += xiDx[i][0]*shapeFncAtIp[j] * val;
          K_y[i][j] += xiDx[i][1]*shapeFncAtIp[j] * val;
          if ( dim == 3 ) {
            K_z[i][j] += xiDx[i][2]*shapeFncAtIp[j] * val;
          }
        }
      }
    }

    elemMat.Resize(numFncs, numFncs*dim);
    elemMat.Init();

    //Here the solution vector is P-vector, the element matrix is of dimension (numFncs1 X dim*numFncs2) 
    for (UInt i=0; i<numFncs; i++) {
      for (UInt j=0; j<numFncs; j++) {
        elemMat[i][j*dim] = K_x[i][j] ;
        elemMat[i][j*dim+1] = K_y[i][j];
	if ( dim == 3 )
	  elemMat[i][j*dim+2] = K_z[i][j];
      }
    }
    //std::cout << "Type: " << formsType_ << "\n" << elemMat << std::endl;

  }

  void PMLTimeInt::SetPosPML(Matrix<Double> & inner, Matrix<Double> & outer)
  {

    pmlFnc_-> SetPosPML( inner, outer );
  }

  PMLMixedInt::PMLMixedInt(std::string type, Double factor, std::string dampingTypePML, Double damp, 
		 bool axi)
    : BaseForm(NULL), formsFactor_(factor)
  {
    name_ = "PMLMixedInt";

    isComplex_ = true;
    isaxi_     = axi;
    //std::cout << type << std::endl;

    nrDofsPerNode_ = 1;
    pmlFnc_    = new PMLBasics( dampingTypePML, damp, type);

    // Set Expression for parser
    mParser_->SetExpr( mHandle_, "f" );
  }


 
  PMLMixedInt::~PMLMixedInt()
  {
  }



  void PMLMixedInt::CalcElementMatrix( Matrix<Complex>& elemMat,
                                  EntityIterator& ent1, 
                                  EntityIterator& ent2 )
  {
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    if(pmlFnc_->GetFormsType() =="mixedPMLMassPPInt"){
      CalcElementMatrixMassPP(ptCoord_, elemMat);
    }else if (pmlFnc_->GetFormsType() =="mixedPMLMassVVInt"){
      CalcElementMatrixMassVV(ptCoord_, elemMat);
    }
  }

  void PMLMixedInt::CalcElementMatrixMassPP(Matrix<Double> & ptCoord, Matrix<Complex> & elemMat) {
    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Vector<Double> CoordAtIP;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs);
    elemMat.Init();
    
    Complex jacDetC;
    Vector<Complex> factorsPML;
    Complex factor, factorPML;
    Double omega = 2 * PI * mParser_->Eval( mHandle_ );

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      ptelem->SetAnsatzFct( ansatzFct1_ );
      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord,
                                           it1_.GetElem() );
        
      ptelem-> GetShFncAtIp(shapeFncAtIp, actIntPt, it1_.GetElem() );
        
      partElemMat.DyadicMult(shapeFncAtIp);
        
      // compute PML factor 
      ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                 ptCoord, it1_.GetElem() );

      //std::cout << "Mass CoordAtIP:\n" << CoordAtIP << std::endl;
      
      pmlFnc_->ComputeFactorPML( factorsPML, jacDetC, CoordAtIP, omega); 

      factorPML = jacDetC * formsFactor_;

      if (isaxi_) 
        factor = 2 * PI * intWeights[actIntPt-1] * jacDet * 
          CoordAtIP[0] * factorPML;
      else 
        factor = intWeights[actIntPt-1] * jacDet * factorPML;
        
      for ( UInt i=0; i<numFncs; i++ ) {
        for ( UInt j=0; j<numFncs; j++ ) {
          elemMat[i][j] += factor * partElemMat[i][j];
        }
      }
    }
    
    if (nrDofsPerNode_ > 1 ) {
     Matrix <Complex> singleDofMass = elemMat;
     UInt singleDofSize = singleDofMass.GetNumRows();

     elemMat.Resize( nrDofsPerNode_* singleDofSize );

     for (UInt i=0; i < singleDofSize; i++)
       for (UInt j=0; j < singleDofSize; j++)
         for (UInt actDof=0; actDof < nrDofsPerNode_; actDof++)
           elemMat[i*nrDofsPerNode_ + actDof][j*nrDofsPerNode_ + actDof] = singleDofMass[i][j]; 
   }

  }

  void PMLMixedInt::CalcElementMatrixMassVV(Matrix<Double> & ptCoord, Matrix<Complex> & elemMat)
  {
     // Extract pointer to reference element and get coordinates
     Vector<Double> CoordAtIP;
     Matrix<Double> tempElemMat;

     ptelem->SetAnsatzFct( ansatzFct1_ );
     UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );

     //check for spectral element approximation if not, throw error
     assert(ansatzFct1_->GetType() == AnsatzFct::SPECTRAL );

     const UInt nrIntPts= ptelem->GetNumIntPoints();    
     const UInt spaceDim = ptelem->GetDim();  
     const Vector<Double> & intWeights = ptelem->GetIntWeights();
     const Vector<Double> * intPoints = ptelem->GetIntPoints();
     Double jacDet;
     Complex jacDetC;
     Vector<Complex> factorsPML;

     elemMat.Resize(spaceDim*numFncs);
     elemMat.Init();

     Matrix<Complex>  partElemMat;
     partElemMat.Resize(spaceDim*numFncs);
     partElemMat.Init();

     Vector<Double> shapeFncAtIp;
     Matrix<Double> JacMat;
     Matrix<Complex> JacMatC;
     Matrix<Complex> JacMatT;
     Matrix<Complex> subMat; 
     JacMat.Resize(spaceDim,spaceDim);
     JacMatC.Resize(spaceDim,spaceDim);
     JacMatT.Resize(spaceDim,spaceDim);
     subMat.Resize(spaceDim,spaceDim);
     JacMat.Init();
     JacMatT.Init();
     subMat.Init();

     Double omega = 2 * PI * mParser_->Eval( mHandle_ );
     for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
       
       jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, it1_.GetElem() );       
       ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, it1_.GetElem() );
       ptelem->CalcJacobianAtIp(JacMat,actIntPt,ptCoord_,it1_.GetElem() );

 	     // compute PML factor 
       ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                  ptCoord, it1_.GetElem() );
	     pmlFnc_->ComputeFactorPML( factorsPML, jacDetC, CoordAtIP, omega);        

       for ( UInt i = 0;i<spaceDim ; i++ ) {
         for ( UInt j=0; j<spaceDim;j++ ) {
           JacMatT[j][i] = JacMat[i][j] * (Complex(1.0,0.0) / factorsPML[j]); 
           JacMatC[i][j] = JacMat[i][j] * (Complex(1.0,0.0) / factorsPML[j]); 
         }
       }
       JacMatC.Mult(JacMatT,subMat);

       for( UInt intPti = 0 ; intPti< shapeFncAtIp.GetSize();intPti++){
         for( UInt intPtj = 0 ; intPtj< shapeFncAtIp.GetSize();intPtj++){
           for ( UInt i = 0;i<spaceDim ; i++ ) {
             for ( UInt j=0; j<spaceDim;j++ ) {
               partElemMat[(intPti*spaceDim)+i][(intPtj*spaceDim)+j] =  intWeights[actIntPt-1] * 
                                                                        subMat[i][j]  * 
                                                                        (Complex(1.0,0.0) / jacDet) *
                                                                        (Complex(1.0,0.0) / jacDetC) *
                                                                        formsFactor_ *
                                                                        shapeFncAtIp[intPti] * shapeFncAtIp[intPtj];
             }
           }
         }
       }

       elemMat += partElemMat;
     }
     //std::cout << "MassPiolaMixedInt_VV Matrix:\n" << elemMat << std::endl;

    //std::cout << elemMat << std::endl << std::endl;
  }

  void PMLMixedInt::SetPosPML(Matrix<Double> & inner, Matrix<Double> & outer)
  {

    pmlFnc_-> SetPosPML( inner, outer );
  }

    //! Constructor
  PMLMixedTimeInt::PMLMixedTimeInt(std::string type, Double factor, std::string dampingTypePML, 
                    Double damp, bool axi)
                  : BaseForm(NULL), factor_(factor)
  {    

    name_ = "PMLMixedTimeInt, " + type;

    isaxi_     = axi;
    formsType_ = type;
    pmlFnc_    = new PMLBasics( dampingTypePML, damp, type);
  }
    
  PMLMixedTimeInt::~PMLMixedTimeInt(){
  }
    
  //! Calculation of stiffmess matrix
  void PMLMixedTimeInt::CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 ){
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    ent_ = ent1;
    
    if ( formsType_ == "PMLVelStiff" ) { 
      CalcElementMatrixVelStiff(ptCoord_, elemMat);
    }   
    else if ( formsType_ == "PMLAccelStiff" ) { 
      CalcElementMatrixAccelStiff(ptCoord_, elemMat);
    }   
    else if ( formsType_ == "PMLAuxVecMass" ) { 
      CalcElementMatrixAuxVecMass(ptCoord_, elemMat);
    }   
    else if ( formsType_ == "PMLStiffPhi" ) { 
      CalcElementMatrixStiffPhi(ptCoord_, elemMat);
    }   
    else if ( formsType_ == "PMLGradR_PhiSigma" ) { 
      CalcElementMatrixGradRV(ptCoord_, elemMat);
    }  
    else
      EXCEPTION("Requested a Type which is not supported by PMLMixedTimeInt " << formsType_);
  }
    
    //! set min/max of x,y,z coordinates form where PML starts and ends
  void PMLMixedTimeInt::SetPosPML(Matrix<Double> & inner, Matrix<Double> & outer)
  {

      pmlFnc_-> SetPosPML( inner, outer );
  }
    
    //! Calculation of stiffmess matrix
  void PMLMixedTimeInt::CalcElementMatrixGradRV(Matrix<Double> & ptCoord, Matrix<Double> & elemMat){

     //check for spectral element approximation if not, throw error
     assert(ansatzFct1_->GetType() == AnsatzFct::SPECTRAL );
    
    // ansatzFct1_ corresponds to functions of pressure P (Taylor Hood: 1st order)
    // ansatzFct2_ corresponds to functions of veloctity V (TaylorHood: 2nd order )
    //UInt numFncs1 = ptelem->GetNumFncs( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct2_ );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
     const Vector<Double> * intPoints = ptelem->GetIntPoints();
    const UInt spaceDim = ptelem->GetDim();  
    Double jacDet;

    Vector<Double> shapeFncAtIp;
    Vector<Double> CoordAtIP;
    Matrix<Double> xiDx;
    Matrix<Double> JacMat;
    Matrix<Double> partElemMat;
    Vector<Double> factorsPML;
    Vector<Double> CoordCenter;
    Vector<Double> LocalCenter;

    LocalCenter.Resize(spaceDim);
    LocalCenter.Init();

    //set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs*spaceDim, numFncs*spaceDim);
    partElemMat.Resize(numFncs*spaceDim, numFncs*spaceDim);
    JacMat.Resize(spaceDim);
    factorsPML.Resize(spaceDim);
    factorsPML.Init();
    JacMat.Init();

    elemMat.Init();
    xiDx.Init();
    shapeFncAtIp.Init();

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      // compute PML factor 
      ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                 ptCoord, it1_.GetElem() );
      // compute Element Center  
      ptelem->Local2GlobalCoord( CoordCenter, LocalCenter,
                                 ptCoord, it1_.GetElem() );
      //pmlFnc_->ComputeTimeFactorPML( factorsPML, CoordAtIP,CoordCenter );  
      pmlFnc_->ComputeTimeFactorPML( factorsPML, CoordAtIP );  

      //get derivatives
      ptelem->SetAnsatzFct( ansatzFct1_ , false);
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet, ent_.GetElem());
      ptelem->CalcJacobianAtIp(JacMat,actIntPt,ptCoord_,ent_.GetElem() );

      //Get the shape functions vector
      ptelem->SetAnsatzFct( ansatzFct2_ , false);
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, ent_.GetElem() );
      
      partElemMat.Init();
      for( UInt intPti = 0 ; intPti< shapeFncAtIp.GetSize();intPti++){
        for( UInt intPtj = 0 ; intPtj< shapeFncAtIp.GetSize();intPtj++){
          for ( UInt i = 0;i<spaceDim ; i++ ) {
            for ( UInt j=0; j<spaceDim;j++ ) {
              partElemMat[(intPti*spaceDim)+i][(intPtj*spaceDim)+j] = //factorsPML[i] * 
                                                                      xiDx[intPti][i] * 
                                                                      //shapeFncAtIp[intPti] * 
                                                                      shapeFncAtIp[intPtj] * 
                                                                      JacMat[i][j] *
                                                                      intWeights[actIntPt-1]*
                                                                      factor_;
            }
          }
        }
      }
      elemMat +=partElemMat;
    }

    //std::cerr << "CalcElementMatrixGradRV Matrix for element #" << it1_.GetElem()->elemNum << ":\n" << elemMat << std::endl;
  }
    

    //! Calculation of mass matrix
  void PMLMixedTimeInt::CalcElementMatrixStiffPhi(Matrix<Double> & ptCoord, Matrix<Double> & elemMat){
    // Extract pointer to reference element and get coordinates
    Vector<Double> CoordAtIP;

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    const UInt spaceDim = ptelem->GetDim();  

    Double jacDet;
    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Matrix<Double> singleMatrix;
    Vector<Double> factorsPML;
    Vector<Double> CoordCenter;
    Vector<Double> LocalCenter;

    LocalCenter.Resize(spaceDim);
    LocalCenter.Init();
    
    elemMat.Resize(numFncs,numFncs*spaceDim);
    partElemMat.Resize(numFncs,numFncs);
    singleMatrix.Resize(numFncs,numFncs);
    elemMat.Init();
    singleMatrix.Init();
    partElemMat.Init();
    

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      // compute PML factor 
      ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                 ptCoord, it1_.GetElem() );
      // compute Element Center  
      ptelem->Local2GlobalCoord( CoordCenter, LocalCenter,
                                 ptCoord, it1_.GetElem() );
      //pmlFnc_->ComputeTimeFactorPML( factorsPML, CoordAtIP,CoordCenter );  
      pmlFnc_->ComputeTimeFactorPML( factorsPML, CoordAtIP );  
       
      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, it1_.GetElem() );
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, it1_.GetElem() );
      
      partElemMat.DyadicMult(shapeFncAtIp);

      if (isaxi_) {
	       ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
			     ptCoord_, it1_.GetElem() );
	       partElemMat *= 2 * PI * intWeights[actIntPt-1] * factor_
	                      * jacDet * CoordAtIP[0];
      }
      else {
	       partElemMat *= intWeights[actIntPt-1]  * jacDet * factor_;
      }
      //now blow the matrix up
      for(UInt i=0; i < numFncs; i++){
        for(UInt j=0; j < numFncs; j++){
          for(UInt d=0; d < spaceDim; d++){
            elemMat[i][(j*spaceDim)+d] += partElemMat[i][j] * factorsPML[d];
          }
        }
      }
    }
  }
    
  void PMLMixedTimeInt::CalcElementMatrixVelStiff(Matrix<Double> & ptCoord, Matrix<Double> & elemMat){
    // Extract pointer to reference element and get coordinates
    Vector<Double> CoordAtIP;
    Matrix<Double> tempElemMat;

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );

    //check for spectral element approximation if not, throw error
    assert(ansatzFct1_->GetType() == AnsatzFct::SPECTRAL );

    const UInt nrIntPts= ptelem->GetNumIntPoints();    
    const UInt spaceDim = ptelem->GetDim();  
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    Double jacDet;
    Vector<Double> CoordCenter;
    Vector<Double> LocalCenter;

    LocalCenter.Resize(spaceDim);
    LocalCenter.Init();

    elemMat.Resize(spaceDim*numFncs);
    elemMat.Init();

    Matrix<Double>  partElemMat;
    partElemMat.Resize(spaceDim*numFncs);
    partElemMat.Init();

    Vector<Double> shapeFncAtIp;
    Matrix<Double> JacMat;
    Matrix<Double> JacMatT;
    Matrix<Double> subMat; 
    Vector<Double> factorsPML(spaceDim);
    JacMat.Resize(spaceDim,spaceDim);
    JacMatT.Resize(spaceDim,spaceDim);
    subMat.Resize(spaceDim,spaceDim);
    JacMat.Init();
    JacMatT.Init();
    subMat.Init();
    factorsPML.Init();


    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      // compute PML factor 
      ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                 ptCoord, it1_.GetElem() );
      // compute Element Center  
      ptelem->Local2GlobalCoord( CoordCenter, LocalCenter,
                                 ptCoord, it1_.GetElem() );
      //pmlFnc_->ComputeTimeFactorPML( factorsPML, CoordAtIP,CoordCenter );  
      pmlFnc_->ComputeTimeFactorPML( factorsPML, CoordAtIP );  
      
      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, it1_.GetElem() );       
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, it1_.GetElem() );
      ptelem->CalcJacobianAtIp(JacMat,actIntPt,ptCoord_,it1_.GetElem() );

      JacMat.Transpose(JacMatT);
      for ( UInt i = 0;i<spaceDim ; i++ ) {
        for ( UInt j=0; j<spaceDim;j++ ) {
         JacMat[i][j] = JacMat[i][j] * factorsPML[i];
        }
      }
      JacMatT.Mult(JacMat,subMat);
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

   // std::cout << "CalcElementMatrixVelStiff Matrix #" << it1_.GetElem()->elemNum << "\n" << elemMat << std::endl;
  }
    
    //! Calculation of mass matrix
  void PMLMixedTimeInt::CalcElementMatrixAuxVecMass(Matrix<Double> & ptCoord, Matrix<Double> & elemMat){
    // Extract pointer to reference element and get coordinates
    //Vector<Double> CoordAtIP;

    //ptelem->SetAnsatzFct( ansatzFct1_ );
    //UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    //const UInt nrIntPts= ptelem->GetNumIntPoints();    
    //const Vector<Double> & intWeights = ptelem->GetIntWeights();
    //const Vector<Double> * intPoints = ptelem->GetIntPoints();
    //const UInt spaceDim = ptelem->GetDim();  

    //Double jacDet;
    //Vector<Double> shapeFncAtIp;
    //Matrix<Double> partElemMat;
    //Matrix<Double> singleMatrix;
    //Matrix<Double> JacMat;

    //JacMat.Resize(spaceDim,spaceDim);
    //singleMatrix.Resize(numFncs);
    //elemMat.Resize(numFncs*spaceDim,numFncs*spaceDim);
    //partElemMat.Resize(numFncs,numFncs);
    //elemMat.Init();
    //singleMatrix.Init();
    //partElemMat.Init();
    //
    //Vector<Double> factorsPML(spaceDim);
    //Vector<Double> JacVec(spaceDim);

    //for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
    //   
    //  partElemMat.Init();


    //  jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, it1_.GetElem() );
    //  ptelem->CalcJacobianAtIp(JacMat,actIntPt,ptCoord_,it1_.GetElem() );
    //  
    //  ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, it1_.GetElem() );
    //  
    //  partElemMat.DyadicMult(shapeFncAtIp);
    //  if (isaxi_) {
	  //     ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
		//	     ptCoord_, it1_.GetElem() );
	  //     partElemMat *= 2 * PI * intWeights[actIntPt-1] * factor_
	  //                    * jacDet * CoordAtIP[0];
    //  }
    //  else {
	  //     partElemMat *= intWeights[actIntPt-1] * factor_;
    //  }
    //  singleMatrix += partElemMat;
    //}
    ////now blow the matrix up
    //for(UInt i=0; i < numFncs; i++){
    //  for(UInt j=0; j < numFncs; j++){
    //    for(UInt d=0; d < spaceDim; d++){
    //      elemMat[(i*spaceDim)+d][(j*spaceDim)+d] = singleMatrix[i][j];
    //    }
    //  }
    //}
    
  }

  //! Calculation of mass matrix
  void PMLMixedTimeInt::CalcElementMatrixAccelStiff(Matrix<Double> & ptCoord, Matrix<Double> & elemMat){
    //check for spectral element approximation if not, throw error
    assert(ansatzFct1_->GetType() == AnsatzFct::SPECTRAL );

    // Extract pointer to reference element and get coordinates
    Vector<Double> CoordAtIP;
    Vector<Double> factorsPML;
    StdVector< Matrix<Double> > tempElemMat;

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();    
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    const UInt spaceDim = ptelem->GetDim();  
    Double jacDet;
    Vector<Double> CoordCenter;
    Vector<Double> LocalCenter;

    LocalCenter.Resize(spaceDim);
    LocalCenter.Init();

    elemMat.Resize(spaceDim*numFncs);
    elemMat.Init();

    tempElemMat.Resize(spaceDim);
    tempElemMat.Init();

    for (UInt i = 0; i < spaceDim; i += 1 ){
      tempElemMat[i].Resize(numFncs);
      tempElemMat[i].Init();
    }

    factorsPML.Resize(spaceDim);
    factorsPML.Init();

    Matrix<Double>  partElemMat;
    partElemMat.Resize(numFncs);
    partElemMat.Init();

    Vector<Double> shapeFncAtIp;

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
      // compute PML factor 
      ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                 ptCoord, it1_.GetElem() );
      // compute Element Center  
      ptelem->Local2GlobalCoord( CoordCenter, LocalCenter,
                                 ptCoord, it1_.GetElem() );
      //pmlFnc_->ComputeTimeFactorPML( factorsPML, CoordAtIP,CoordCenter );  
      pmlFnc_->ComputeTimeFactorPML( factorsPML, CoordAtIP );  

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, it1_.GetElem() );       
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt, it1_.GetElem() );
      partElemMat.DyadicMult(shapeFncAtIp);

      if (isaxi_) {
	       ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
	        		    ptCoord_, it1_.GetElem() );
	       partElemMat *= 2 * PI * intWeights[actIntPt-1] * factor_
	         * jacDet * CoordAtIP[0];
      }
      else {
        partElemMat *= intWeights[actIntPt-1] * jacDet * factor_;
      }

      for (UInt i = 0; i < spaceDim; i += 1 ){
        tempElemMat[i] += (partElemMat * factorsPML[i]);
      }
    }
    //Blowing up the element matrix to number of velocity components 
    for (UInt i=0; i < numFncs; i++){
      for (UInt j=0; j < numFncs; j++){
        for (UInt actDof = 0; actDof < spaceDim ; actDof++){
          elemMat[i*spaceDim + actDof][j*spaceDim + actDof] = tempElemMat[actDof][i][j]; 
        }
      }
    }

    //std::cout << "Damped ElemMatStiff for substitution variable:\n" << elemMat << std::endl;

   }
} // end namespace CoupledField
