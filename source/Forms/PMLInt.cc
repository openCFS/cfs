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

    if ( pmlFnc_->GetFormsType() =="laplaceInt" ) {
      CalcElementMatrixStiff(ptCoord_, elemMat);
    }
    else {
      CalcElementMatrixMass(ptCoord_, elemMat);
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
	for (UInt i=0; i<xiDx.GetSizeCol(); i++) {
	  for (UInt j=0; j<xiDx.GetSizeRow(); j++) {
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
     UInt singleDofSize = singleDofMass.GetSizeRow();

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
        (*error) << "BDBInt::CalcElementMatrix: Encountered "
                 << "negative Jacobian determinant!";
        Error( __FILE__, __LINE__ );
      }
      
      // Special things must be done in the axi-symmetric case
      if ( isaxi_ ) {
        Vector<Double> ShpFncAtIp;
        ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );
        
        CoordAtIP = ptCoord_ * ShpFncAtIp;
        jacDet *= 2 * PI * CoordAtIP[0];
      }
      
      // Compute the matrix product D * B and store as intermediate matrix
      //      dbMatC.Resize( dMat.GetSizeRow(), bMat.GetSizeCol() );
      dbMatC = dMatC * bMatC;
      
      // We now compute B^T * D * B and scale it by the determinant
      // of the Jacobian and the weight of the current integration
      // point. The result is added to the element matrix.
      fac = jacDet * intWeights[actIntPt-1] * jacDetC;
      for ( UInt k = 0; k < bMat.GetSizeRow(); k++ ) {
        ptr1 =  bMatC[k];
        ptr2 = dbMatC[k];
        for ( UInt i = 0; i < bMatC.GetSizeCol(); i++ ) {
          aux = fac * ptr1[i];
          for ( UInt j = 0; j < dbMatC.GetSizeCol(); j++ ) {
            elemMat[i][j] += aux * ptr2[j];
          }
        }
      }
    }
  }
  

  void MechPMLInt::calcBMatPML( Matrix<Double>& bMat, Vector<Double>& CoordAtIP, 
                                Matrix<Complex>& bMatC, Complex& jacDetC)
  {

    bMatC.Resize( bMat.GetSizeRow(),  bMat.GetSizeCol() );
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


} // end namespace CoupledField
