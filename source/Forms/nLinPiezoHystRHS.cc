// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <math.h>

#include "nLinPiezoHystRHS.hh"
#include "mechStressStrain.hh"
#include "Utils/coordSystem.hh"
#include "Domain/domain.hh"
#include "Utils/SmoothSpline.hh"
#include "CoupledPDE/PiezoCoupling.hh"

namespace CoupledField {



  // ==================================================================
  // piezoelectric polarization: 
  // ==================================================================

  // base class
  PiezoPolarizationRhsInt::PiezoPolarizationRhsInt( BaseMaterial* matDataCouple, 
                                                    BaseMaterial* matDataMech, 
                                                    BaseMaterial* matDataElec, 
                                                    SubTensorType type)
    : LinearForm() {

    matDataCouple_ = matDataCouple;
    matDataMech_   = matDataMech;
    matDataElec_   = matDataElec;
    
    subTensorType_ = type;

    isaxi_ = false;
    if ( subTensorType_ == AXI )
      isaxi_ = true;

    piezoBilinearForm_ = new linPiezoCoupling(matDataCouple, type); 
  }

  PiezoPolarizationRhsInt::~PiezoPolarizationRhsInt() {

    delete EfieldOp_;
    delete piezoBilinearForm_;
    delete EfieldPrevOp_;
    delete mechStrainOp_;
    //delete matDataCouple_;
    //delete matDataMech_;
    //delete matDataElec_;
  }


  void PiezoPolarizationRhsInt::Set4NonLinMaterial( Grid* ptGrid, 
                             StdPDE* ptPDE2Elec,
                             shared_ptr<EqnMap> eqnMapElec,
                             shared_ptr<ResultInfo> resultElec,
                                                    PiezoCoupling* ptPiezoCoupling) {

    // pointer to piezo-coupling object
    piezoCoupling_ =  ptPiezoCoupling;

    // electric field operator
    EfieldOp_ =  new GradientFieldOp<Double>(ptGrid, ptPDE2Elec, 
                                             eqnMapElec, *solElec_, resultElec->fctType, isaxi_,
                                             coordUpdate_);

    EfieldPrevOp_ =  new GradientFieldOp<Double>(ptGrid, ptPDE2Elec, 
                                                 eqnMapElec, *solPrevElec_, resultElec->fctType, isaxi_,
                                                 coordUpdate_);

    // mechanical strain operator
    mechStrainOp_ = 
      new MechStressStrain<Double>( matDataMech_ , subTensorType_);
  }


  // ==================================================================
  // piezoelectric polarization: RHS of electric equation
  // ==================================================================

  PiezoPolarizationElecRhsInt::PiezoPolarizationElecRhsInt( BaseMaterial* matDataCouple, 
                                                            BaseMaterial* matDataMech, 
                                                            BaseMaterial* matDataElec, 
                                                            SubTensorType type)
    : PiezoPolarizationRhsInt( matDataCouple, matDataMech, matDataElec, type ) {

    name_  = "PiezoPolarisationElecRhsInt";

  }


  PiezoPolarizationElecRhsInt::~PiezoPolarizationElecRhsInt() {
  }


  void PiezoPolarizationElecRhsInt::CalcElemVector( Vector<Double> & elemVec,
						    EntityIterator& ent ) 
    
  {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
   
     // compute the correct material tensor
    Vector<Double> LCoord, Efield, mechStrain;
    Vector<Double> EfieldPrev, mechStrainPrev, elecD;
    Matrix<Double> dMat1, dMat2;

    // compute electric field in the midpoint of element
    ptelem->GetCoordMidPoint(LCoord);
    EfieldOp_->CalcElemGradField( Efield, ent, LCoord, 1);
    EfieldPrevOp_->CalcElemGradField( EfieldPrev, ent, LCoord, 1);

    // calc mechanical strain
    Matrix<Double> elemDispl, elemDisplPrev;
    solMech_->GetElemSolutionAsMatrix( elemDispl, ent);
    mechStrainOp_->SetActElemSol(elemDispl);
    mechStrainOp_->SetIntPoint(LCoord);
    mechStrainOp_->CalcStrainVec( mechStrain, 1, ent);

    solPrevMech_->GetElemSolutionAsMatrix( elemDisplPrev, ent);
    mechStrainOp_->SetActElemSol(elemDisplPrev);
    mechStrainOp_->SetIntPoint(LCoord);
    mechStrainOp_->CalcStrainVec( mechStrainPrev, 1, ent);

    //    std::cout << "EfieldPrev:  \n" << EfieldPrev << std::endl;
    // compute the material tensors

    std::string matType = "ElecRHS1";
    piezoCoupling_->GetNonlinMaterialTensor( dMat1, elecD,
                                             matType, 
                                             matDataMech_, 
                                             matDataElec_,
                                             matDataCouple_, 
                                             subTensorType_, 
                                             Efield, mechStrain, 
                                             EfieldPrev, mechStrainPrev, 
                                             ent );

    //    std::cout << "ELEC-RHS1: \n " << dMat1 << std::endl;

    matType = "ElecRHS2";
    piezoCoupling_->GetNonlinMaterialTensor( dMat2, elecD,
                                             matType, 
                                             matDataMech_, 
                                             matDataElec_,
                                             matDataCouple_, 
                                             subTensorType_, 
                                             Efield, mechStrain, 
                                             EfieldPrev, mechStrainPrev, 
                                             ent );

    //    std::cout << "ELEC-RHS2: \n " << elecD << std::endl;

    Double jacDet, aux, *ptr1, *ptr2;
    Matrix<Double> bMat, bMatTrans, dbMat;
    Matrix<Double> elemMat;
    Vector<Double> helpVec;

    elemMat.Resize( numFncs );
    elemMat.Init();
    dbMat.Resize( dMat1.GetNumRows(), numFncs );
    elemVec.Resize ( numFncs );
    elemVec.Init();

    piezoBilinearForm_->ExtractElemInfo( ent );

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

      // Setup the B matrix for current integration point
      piezoBilinearForm_->CalcBMat( bMat, actIntPt, ptCoord_ );

      // Compute Jacobian for integration point
      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_,ent.GetElem() );

      if ( isaxi_ ) {
        Vector<Double>  ShpFncAtIp;
        Double rad = 0.0;

        ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent.GetElem() );
        for ( UInt i = 0; i < numFncs; i++ ) {
          rad += ptCoord_[0][i] * ShpFncAtIp[i];
        }
        jacDet *= 2.0 * PI * rad;
      }
      jacDet *= intWeights[actIntPt-1];

      //part 2 of RHS
      bMat.Transpose( bMatTrans);
      helpVec = bMatTrans * elecD;
      helpVec *= jacDet;
      elemVec -= helpVec;
  
      // part 1 of RHS
      // Compute the matrix product D * B and store as intermediate matrix
      dMat1.Mult( bMat, dbMat );

      // We now compute B^T * D * B and scale it by the determinant
      // of the Jacobian and the weight of the current integration
      // point. The result is added to the element matrix.
      for ( UInt k = 0; k < bMat.GetNumRows(); k++ ) {
        ptr1 = bMat[k];
        ptr2 = dbMat[k];
        for ( UInt i = 0; i < bMat.GetNumCols(); i++ ) {
          aux = jacDet * ptr1[i];
          for ( UInt j = 0; j < dbMat.GetNumCols(); j++ ) {
            elemMat[i][j] += aux * ptr2[j];
          }
        }
      }
    }
    
    Matrix<Double> tmp;
    Vector<Double> pot;
    solPrevElec_->GetElemSolutionAsMatrix(  tmp, ent);
    tmp.ConvertToVec_AppendRows( pot );

    //add part 1 and part 2!
    elemVec -= elemMat * pot; 

    //    std::cout << "elecRHS:\n" << elemVec << std::endl;
  }


  // ==================================================================
  // piezoelectric polarization: RHS of mechanic equation
  // ==================================================================

  PiezoPolarizationMechRhsInt::PiezoPolarizationMechRhsInt( BaseMaterial* matDataCouple, 
                                                            BaseMaterial* matDataMech, 
                                                            BaseMaterial* matDataElec, 
                                                            SubTensorType type)
    : PiezoPolarizationRhsInt( matDataCouple, matDataMech, matDataElec, type ) {

    name_  = "PiezoPolarisationMechRhsInt";

  }


  PiezoPolarizationMechRhsInt::~PiezoPolarizationMechRhsInt()
  {
    //delete matDataMech_; 
  }


  void PiezoPolarizationMechRhsInt::CalcElemVector( Vector<Double> & elemVec,
						    EntityIterator& ent ) 
    
  {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  

    Double jacDet, aux;
    Matrix<Double> bMat, aMat, dMat;
    Matrix<Double> dbMat, elemMat;

    elemMat.Resize( numFncs * piezoBilinearForm_->getNumDofsA(), 
                    numFncs * piezoBilinearForm_->getNumDofsB() ); 
    elemMat.Init();

    // compute the correct material tensor
    Vector<Double> LCoord, Efield, mechStrain;
    Vector<Double> EfieldPrev, mechStrainPrev, elecD;

    // compute electric field in the midpoint of element
    ptelem->GetCoordMidPoint(LCoord);
    EfieldOp_->CalcElemGradField( Efield, ent, LCoord, 1);
    EfieldPrevOp_->CalcElemGradField( EfieldPrev, ent, LCoord, 1);

    // calc mechanical strain
    Matrix<Double> elemDispl, elemDisplPrev;
    solMech_->GetElemSolutionAsMatrix( elemDispl, ent);
    mechStrainOp_->SetActElemSol(elemDispl);
    mechStrainOp_->SetIntPoint(LCoord);
    mechStrainOp_->CalcStrainVec( mechStrain, 1, ent);

    solPrevMech_->GetElemSolutionAsMatrix( elemDisplPrev, ent);
    mechStrainOp_->SetActElemSol(elemDisplPrev);
    mechStrainOp_->SetIntPoint(LCoord);
    mechStrainOp_->CalcStrainVec( mechStrainPrev, 1, ent);

    // compute the nonlinear coupling tensor
    std::string matType = "MechRHS";
    piezoCoupling_->GetNonlinMaterialTensor( dMat, elecD,
                                             matType, 
                                             matDataMech_, 
                                             matDataElec_,
                                             matDataCouple_, 
                                             subTensorType_, 
                                             Efield, mechStrain, 
                                             EfieldPrev, mechStrainPrev, 
                                             ent );

    //    std::cout << "MECH-RHS: \n " << dMat << std::endl;


    piezoBilinearForm_->ExtractElemInfo( ent );

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

      // Setup the A matrix for current integration point
      piezoBilinearForm_->calcAMat( aMat, actIntPt, ptCoord_ );

      // Setup the B matrix for current integration point
      piezoBilinearForm_->CalcBMat( bMat, actIntPt, ptCoord_ );

      // Compute Jacobian for integration point
      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_,ent.GetElem() );

      if ( isaxi_ ) {
        Vector<Double>  ShpFncAtIp;
        Double rad = 0.0;

        ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent.GetElem() );
        for ( UInt i = 0; i < numFncs; i++ ) {
          rad += ptCoord_[0][i] * ShpFncAtIp[i];
        }
        jacDet *= 2.0 * PI * rad;
      }

      jacDet *= intWeights[actIntPt-1];

      // Compute the matrix product D * B and store as intermediate matrix
      dbMat.Resize( dMat.GetNumRows(), bMat.GetNumCols() );
      dMat.Mult( bMat, dbMat );

      // We now compute A * D * B and scale it by the determinant
      // of the Jacobian 
      for ( UInt i = 0; i < aMat.GetNumRows(); i++ ) {
        for ( UInt j = 0; j < dbMat.GetNumCols(); j++ ) {

          // Compute entry (i,j) of A * D * B
          aux = 0.0;
          for ( UInt k = 0; k < aMat.GetNumCols(); k++ ) {
            aux += aMat[i][k] * dbMat[k][j];
          }

          // Scale result and add to corresponding entry
          // of element matrix
          elemMat[i][j] -= aux * jacDet;
        }
      }
    }
    Matrix<Double> tmp;
    Vector<Double> pot;
    solPrevElec_->GetElemSolutionAsMatrix( tmp, ent);
    tmp.ConvertToVec_AppendRows( pot );
    elemVec = elemMat * pot; 
    
    // std::cout << "mechRHS:\n" << elemVec << std::endl;
  }


} // end of namespace
