// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>

//#include "Forms/forms_header.hh"
#include "CoupledPDE/PiezoCoupling.hh"
#include "Domain/resultInfo.hh"
#include "Elements/basefe.hh"
#include "Forms/adbInt.hh"
#include "Forms/bdbInt.hh"
#include "Forms/gradfieldop.hh"
#include "Forms/linGradBDBInt.hh"
#include "Forms/mechStressStrain.hh"
#include "Forms/nLinPiezoMicro.hh"
#include "General/environment.hh"
#include "Utils/nodestoresol.hh"


namespace CoupledField {


  // ============
  //   Constructor of coupling bilinear form in case of micro-piezo-model
  // ============
class BaseMaterial;
class EqnMap;
class Grid;
class StdPDE;
struct Elem;

  nLinMicroPiezoCouple::nLinMicroPiezoCouple( BaseMaterial* matDataCouple,
                                              BaseMaterial* matDataMech,
                                              BaseMaterial* matDataElec, 
                                              SubTensorType type,
                                              bool isMech ) 
    : linPiezoCoupling(matDataCouple, type){

    name_ = "nLinMicroPiezoCouple";
    isSolDependent_ = true;

    // true, if coupling part belongs to mechanical equation    
    isMech_ = isMech;

    matDataCouple_ = matDataCouple;
    matDataMech_   = matDataMech;
    matDataElec_   = matDataElec;
    
    subTensorType_ = type;

    if ( type == FULL ) {
      numDofsA_  = 3;
      numDofsB_  = 1;
      matDimRow_ = 6;
      matDimCol_ = 3;
    }
    else if ( type == PLANE_STRAIN || type == PLANE_STRESS ) {
      numDofsA_  = 2;
      numDofsB_  = 1;
      matDimRow_ = 3;
      matDimCol_ = 2;
    }
    else if ( type == AXI ) {
      numDofsA_  = 2;
      numDofsB_  = 1;
      matDimRow_ = 4;
      matDimCol_ = 2;
      isaxi_     = true;
    }
  }

    //! Destructor
  nLinMicroPiezoCouple::~nLinMicroPiezoCouple() {
    delete EfieldOp_;
    delete EfieldPrevOp_;
    delete mechStrainOp_;
  }



  void nLinMicroPiezoCouple::Set4NonLinMaterial(Grid* ptGrid, 
                                               StdPDE* ptPDE2Elec,
                                               shared_ptr<EqnMap> eqnMapElec,
                                               shared_ptr<ResultInfo> resultElec,
                                               PiezoCoupling* ptPiezoCoupling) 
  {
    // pointer to piezo-coupling object
    piezoCoupling_ =  ptPiezoCoupling;

    // electric field operator
    EfieldOp_ =  new GradientFieldOp<Double>(ptGrid, ptPDE2Elec, 
                                             eqnMapElec, *solElec_, resultElec->fctType, isaxi_,
                                             coordUpdate_);

    // electric field operator with previous solution
    EfieldPrevOp_ =  new GradientFieldOp<Double>(ptGrid, ptPDE2Elec, 
                                                 eqnMapElec, *solPrevElec_, resultElec->fctType, isaxi_,
                                                 coordUpdate_);

    // mechanical strain operator
    mechStrainOp_ = 
      new MechStressStrain<Double>( matDataMech_ , subTensorType_);
  }

  // ====================
  //   calcElementMatrix
  // ====================
  void nLinMicroPiezoCouple::CalcElementMatrix( Matrix<Double>& elemMat,
                                               EntityIterator& ent1, 
                                               EntityIterator& ent2 ) {
    
    ent1_ = ent1;
      
    // get displacements of element
    solMech_->GetElemSolutionAsMatrix( elemDispl_, ent1);
    solPrevMech_->GetElemSolutionAsMatrix( elemDisplPrev_, ent1);

    // get elecPot of element
    solElec_->GetElemSolution( elemPot_, ent2);
    solPrevElec_->GetElemSolution( elemPotPrev_, ent2);

    Vector<Double> LCoord, Efield;

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 ); 
    ExtractElemInfo( ent2 ); 

    Matrix<double> matCouple;
    ADBInt::CalcElementMatrix( matCouple, ent1, ent2 );

    if ( isMech_ ) {
      elemMat = matCouple;
    }
    else {
      matCouple.Transpose( elemMat ); 
    }
    //std::cout << "CoupleMatrix: \n" <<  matCouple << std::endl;
  }


  void nLinMicroPiezoCouple::calcDMat(Matrix<Double> & dMat, const Elem* elem ) {

    Vector<Double> LCoord, Efield, mechStrain;
    Vector<Double> EfieldPrev, mechStrainPrev, elecD;

    //get mid-point of element
    ptelem->GetCoordMidPoint(LCoord);

    // calc E-Field:
    EfieldOp_->CalcElemGradField( Efield, ent1_, LCoord, 1);
    EfieldPrevOp_->CalcElemGradField( EfieldPrev, ent1_, LCoord, 1);

    // calc mechanical strain
    mechStrainOp_->SetActElemSol(elemDispl_);
    mechStrainOp_->SetIntPoint(LCoord);
    mechStrainOp_->CalcStrainVec( mechStrain, 1, ent1_ );

    mechStrainOp_->SetActElemSol(elemDisplPrev_);
    mechStrainOp_->SetIntPoint(LCoord);
    mechStrainOp_->CalcStrainVec( mechStrainPrev, 1, ent1_ );

    // compute the nonlinear coupling tensor
    std::string matType;
    if ( isMech_ ) 
      matType = "CouplingTensorMechEQ";
    else 
      matType = "CouplingTensorElecEQ";

    piezoCoupling_->GetMaterialTensorMicroPiezo( dMat, elecD,
                                                 matType, 
                                                 matDataMech_, 
                                                 matDataElec_,
                                                 matDataCouple_, 
                                                 subTensorType_, 
                                                 Efield, mechStrain, 
                                                 EfieldPrev, mechStrainPrev, 
                                                 ent1_ );
    //    std::cout << matType << "\n" << dMat << std::endl;
  }
 

  // ============
  //   Constructor of mechanical bilinear form in case of micro-piezo-model
  // ============
  nLinMicroPiezoMech::nLinMicroPiezoMech( BaseMaterial* matDataCouple,
                                          BaseMaterial* matDataMech,
                                          BaseMaterial* matDataElec,                                       
                                          SubTensorType type) 
    : linElastInt(matDataElec, type){
    
    name_ = "nLinMicroPiezoMech";
    isSolDependent_ = true;

    matDataCouple_ = matDataCouple;
    matDataMech_   = matDataMech;
    matDataElec_   = matDataElec;
    
    subTensorType_ = type;

    if ( type == AXI )
      isaxi_ = true;
  }

    // Destructor
  nLinMicroPiezoMech::~nLinMicroPiezoMech() {
    delete EfieldOp_;
    delete EfieldPrevOp_;
    delete mechStrainOp_;
  }



  void nLinMicroPiezoMech::Set4NonLinMaterial(Grid* ptGrid, 
                                              StdPDE* ptPDE2Elec,
                                              shared_ptr<EqnMap> eqnMapElec,
                                              shared_ptr<ResultInfo> resultElec,
                                              PiezoCoupling* ptPiezoCoupling) 
  {

    // pointer to piezo-coupling object
    piezoCoupling_ =  ptPiezoCoupling;

    // electric field operator
    EfieldOp_ =  new GradientFieldOp<Double>(ptGrid, ptPDE2Elec, 
                                             eqnMapElec, *solElec_, resultElec->fctType, isaxi_,
                                             coordUpdate_);
    // electric field operator
    EfieldPrevOp_ =  new GradientFieldOp<Double>(ptGrid, ptPDE2Elec, 
                                                 eqnMapElec, *solPrevElec_, resultElec->fctType, isaxi_,
                                                 coordUpdate_);

    // mechanical strain operator
    mechStrainOp_ = 
      new MechStressStrain<Double>( matDataMech_ , subTensorType_);
  }

  // ====================
  //   calcElementMatrix
  // ====================
  void nLinMicroPiezoMech::CalcElementMatrix( Matrix<Double>& elemMat,
                                             EntityIterator& ent1, 
                                             EntityIterator& ent2 ) {

    ent1_ = ent1;
      
    // get displacements of element
    solMech_->GetElemSolutionAsMatrix( elemDispl_, ent1);
    solPrevMech_->GetElemSolutionAsMatrix( elemDisplPrev_, ent1);

    // get elecPot of element
    solElec_->GetElemSolution( elemPot_, ent2);

    Vector<Double> LCoord, Efield;

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 ); 
    ExtractElemInfo( ent2 ); 

    BDBInt::CalcElementMatrix( elemMat, ent1, ent2 );
    //    std::cout << "Mech Matrix:\n " << elemMat << std::endl;

  }


  void nLinMicroPiezoMech::calcDMat(Matrix<Double> & dMat, const Elem* elem) {

    Vector<Double> LCoord, Efield, mechStrain;
    Vector<Double> EfieldPrev, mechStrainPrev, elecD;

    //get mid-point of element
    ptelem->GetCoordMidPoint(LCoord);

    // calc E-Field:
    EfieldOp_->CalcElemGradField( Efield, ent1_, LCoord, 1);
    EfieldPrevOp_->CalcElemGradField( EfieldPrev, ent1_, LCoord, 1);

    // calc mechanical strain
    mechStrainOp_->SetActElemSol(elemDispl_);
    mechStrainOp_->SetIntPoint(LCoord);
    mechStrainOp_->CalcStrainVec( mechStrain, 1, ent1_ );
    mechStrainOp_->SetActElemSol(elemDisplPrev_);
    mechStrainOp_->SetIntPoint(LCoord);
    mechStrainOp_->CalcStrainVec( mechStrainPrev, 1, ent1_ );

    // compute the nonlinear coupling tensor
    std::string matType = "MechTensor";
    piezoCoupling_->GetMaterialTensorMicroPiezo( dMat, elecD,
                                                 matType, 
                                                 matDataMech_, 
                                                 matDataElec_,
                                                 matDataCouple_, 
                                                 subTensorType_, 
                                                 Efield, mechStrain, 
                                                 EfieldPrev, mechStrainPrev, 
                                                 ent1_ );

    //std::cout << matType << "\n " << dMat << std::endl;
  }

  // ============
  //   Constructor of electric form in case of micro-piezo-model
  // ============
  nLinMicroPiezoElec::nLinMicroPiezoElec( BaseMaterial* matDataCouple,
                                          BaseMaterial* matDataMech,
                                          BaseMaterial* matDataElec,                                       
                                          SubTensorType type) 
    : linGradBDBInt(matDataElec, ELEC_PERMITTIVITY, type){
    
    name_ = "nLinMicroPiezoElec";
    isSolDependent_ = true;

    matDataCouple_ = matDataCouple;
    matDataMech_   = matDataMech;
    matDataElec_   = matDataElec;
    
    subTensorType_ = type;

    if ( type == AXI )
      isaxi_ = true;
  }

    // Destructor
  nLinMicroPiezoElec::~nLinMicroPiezoElec() {
    delete EfieldOp_;
    delete mechStrainOp_;
  }



  void nLinMicroPiezoElec::Set4NonLinMaterial(Grid* ptGrid, 
                                             StdPDE* ptPDE2Elec,
                                             shared_ptr<EqnMap> eqnMapElec,
                                             shared_ptr<ResultInfo> resultElec,
                                             PiezoCoupling* ptPiezoCoupling) 
  {

    // pointer to piezo-coupling object
    piezoCoupling_ =  ptPiezoCoupling;

    // electric field operator
    EfieldOp_ =  new GradientFieldOp<Double>(ptGrid, ptPDE2Elec, 
                                             eqnMapElec, *solElec_, resultElec->fctType, isaxi_,
                                             coordUpdate_);
    // electric field operator
    EfieldPrevOp_ =  new GradientFieldOp<Double>(ptGrid, ptPDE2Elec, 
                                                 eqnMapElec, *solPrevElec_, resultElec->fctType, isaxi_,
                                                 coordUpdate_);

    // mechanical strain operator
    mechStrainOp_ = 
      new MechStressStrain<Double>( matDataMech_ , subTensorType_);
  }

  // ====================
  //   calcElementMatrix
  // ====================
  void nLinMicroPiezoElec::CalcElementMatrix( Matrix<Double>& elemMat,
                                             EntityIterator& ent1, 
                                             EntityIterator& ent2 ) {

    ent1_ = ent1;
      
    // get displacements of element
    solMech_->GetElemSolutionAsMatrix( elemDispl_, ent1);
    solPrevMech_->GetElemSolutionAsMatrix( elemDisplPrev_, ent1);

    // get elecPot of element
    solElec_->GetElemSolution( elemPot_, ent2);

    Vector<Double> LCoord, Efield;

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 ); 
    ExtractElemInfo( ent2 ); 

    BDBInt::CalcElementMatrix( elemMat, ent1, ent2 );
    elemMat *= -1.0;
    //std::cout << "Elec Matrix:\n " << elemMat << std::endl;

  }


  void nLinMicroPiezoElec::calcDMat(Matrix<Double> & dMat, const Elem* elem ) {

    Vector<Double> LCoord, Efield, mechStrain;
    Vector<Double> EfieldPrev, mechStrainPrev, elecD;

    //get mid-point of element
    ptelem->GetCoordMidPoint(LCoord);

    // calc E-Field:
    EfieldOp_->CalcElemGradField( Efield, ent1_, LCoord, 1);
    EfieldPrevOp_->CalcElemGradField( EfieldPrev, ent1_, LCoord, 1);

    // calc mechanical strain
    mechStrainOp_->SetActElemSol(elemDispl_);
    mechStrainOp_->SetIntPoint(LCoord);
    mechStrainOp_->CalcStrainVec( mechStrain, 1, ent1_ );
    mechStrainOp_->SetActElemSol(elemDisplPrev_);
    mechStrainOp_->SetIntPoint(LCoord);
    mechStrainOp_->CalcStrainVec( mechStrainPrev, 1, ent1_ );
//     std::cout <<"EfieldPrev:\n"  << EfieldPrev << std::endl;
//     std::cout <<"Efield:\n"  << Efield << std::endl;

    // compute the nonlinear coupling tensor
    std::string matType = "ElecTensor";
    piezoCoupling_->GetMaterialTensorMicroPiezo( dMat, elecD,
                                                 matType, 
                                                 matDataMech_, 
                                                 matDataElec_,
                                                 matDataCouple_, 
                                                 subTensorType_, 
                                                 Efield, mechStrain, 
                                                 EfieldPrev, mechStrainPrev, 
                                                 ent1_ );

    // std::cout << matType << "\n " << dMat << std::endl;
  }
}

