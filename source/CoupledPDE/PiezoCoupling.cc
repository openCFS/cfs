// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <assert.h>
#include <stddef.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include <utility>

#include "CoupledPDE/BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Composite.hh"
#include "Domain/ansatzFct.hh"
#include "Domain/domain.hh" 
#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "Domain/surfElem.hh"
#include "Driver/assemble.hh"
#include "Driver/formsContext.hh"
#include "Driver/singleDriver.hh"
#include "Driver/transientdriver.hh"
#include "Elements/basefe.hh"
#include "Forms/FlatShellPiezoInt.hh"
#include "Forms/PiezoStressStrain.hh"
#include "Forms/baseForm.hh"
#include "Forms/elecchargeop.hh"
#include "Forms/gradfieldop.hh"
#include "Forms/linPiezoCoupling.hh"
// integrator (bi-)linear forms
#include "Forms/mechStressStrain.hh"
#include "Forms/nLinPiezoHyst.hh"
#include "Forms/nLinPiezoHystRHS.hh"
#include "Forms/nLinPiezoMicro.hh"
#include "Forms/nLinPiezoMicroRHS.hh"
#include "General/Enum.hh"
#include "General/exception.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Materials/baseMaterial.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/mechPDE.hh"
#include "PiezoCoupling.hh"
#include "Utils/StdVector.hh"
#include "Utils/basenodestoresol.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/result.hh"
#include "Utils/tools.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  PiezoCoupling::PiezoCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                PtrParamNode paramNode  )
    : BasePairCoupling( pde1, pde2, paramNode, "piezoDirect")
  {
    materialClass_ = PIEZO;

    // determine subtype from mechanic pde
    pde1_->GetParamNode()->GetValue( "subType", subType_ );


    // read nonlinearity
    nonLin_ = false;
    nonLinHysteresis_ = false;
    nonLinPiezoMicroHF_ = false;
    nonLinPiezoCoupling_=false;

    // Initialize nonlinearities
    InitNonLin();


  }


  // **************
  //   Destructor
  // **************
  PiezoCoupling::~PiezoCoupling() {
  }


  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void PiezoCoupling::InitNonLin() {

    BasePairCoupling::InitNonLin();

    // Check, if nonlinear type is the same for all regions
    if( regionNonLinTypes_.size() > 1 ) {
      std::map<std::string, NonLinType>::iterator it;
      it = nonLinTypes_.begin();
      NonLinType firstType = it->second;
      for( ; it != nonLinTypes_.end(); it++ ) {
        if( it->second != firstType ) {
          EXCEPTION( "Non-linearity should be the same for all regions!" );
        }
      }
    }


    //now do coupling specifics
    std::map<std::string, NonLinType>::iterator it;
    for ( it=nonLinTypes_.begin() ; it != nonLinTypes_.end(); it++ ) {
      if ( (*it).second == HYSTERESIS ) {
        nonLin_ = true;
        nonLinHysteresis_ = true;
      }
      else if ( (*it).second ==  PIEZO_MICRO_HF ) {
        nonLin_ = true;
        nonLinPiezoMicroHF_ = true;
      }
      else if ( (*it).second == GEOMETRIC ) {
        nonLin_ = true;
      }
    }

    if(nonLin_) {
      nonLinPiezoCoupling_=true;
      if ( !nonLinHysteresis_ ) {
        pde1_->SetNonLinearity(nonLin_);
        pde2_->SetNonLinearity(nonLin_);
      }
    }

  }


  // ***************
  //   CalcResults
  // ***************
  void PiezoCoupling::CalcResults( shared_ptr<BaseResult> result ) {

    // neighborregion for electric charge
    RegionIdType neighbor = surfNeighborRegions_[result];

    switch (result->GetResultInfo()->resultType ) {
    case ELEC_CHARGE:
      if ( isComplex_ ) {
      	CalcCharges<Complex>( result, neighbor );
      } else{
        CalcCharges<Double>( result, neighbor );
      }
      break;

    case MECH_STRESS:
      if ( nonLinPiezoMicroHF_ ) 
        CalcStressStrain<Double>( result );
      else {
        if ( isComplex_ ) {
          CalcStressStrain<Complex>( result );
        } else {
          CalcStressStrain<Double>( result );
        }
      }
      break;

    case MECH_STRAIN:
      if ( isComplex_ ) {
        CalcStressStrain<Complex>( result );
      } else {
        CalcStressStrain<Double>( result );
      }
      break;

    case VON_MISES_STRESS:
      if(isComplex_)
        CalcStressStrain<Complex>(result);
      else
        CalcStressStrain<Double>(result);
      break;


    case MECH_STRAIN_IRR:
      if ( isComplex_ ) {
        EXCEPTION("Ireversible Strain makes no sense in Harmonic analysis");
      } else {
        CalcStrainIrr( result );
      }
      break;

    case ELEC_POLARIZATION:
      if ( isComplex_ ) {
        EXCEPTION("Electric Polarization makes no sense in Harmonic analysis");
      } else {
        CalcElecPolarization( result );
        
      }
      break;

    case ELEC_FLUX_DENSITY:
      if ( isComplex_ ) {
         CalcElecFluxDensity<Complex>( result );
      } else {
        CalcElecFluxDensity<Double>( result );
      }
      break;

    default: 
      WARN( "Resulttype not computable by piezoelectric coupling" );
    }
  }


  template <class TYPE>
  void PiezoCoupling::CalcStressStrain(shared_ptr<BaseResult> res)
  {
    SolutionType resultType = res->GetResultInfo()->resultType;
    Result<TYPE> &  actRes = dynamic_cast<Result<TYPE>&>(*res);
    EntityIterator it = actRes.GetEntityList()->GetIterator();

    // get the materials for the subdomain
    std::map<RegionIdType, BaseMaterial*> mechMat = pde1_->getPDEMaterialData();
    BaseMaterial* matPiezo = materials_[it.GetElem()->regionId];
    BaseMaterial* mechMatSD = mechMat[it.GetElem()->regionId];

    bool isMicroModel = matPiezo->GetMicroPiezoModel() != NULL;

    if(!isMicroModel)
    {
      // this is a "shared" mechanic implementation and is more accurate as we do not use the
      // midpoints but integration points. Essential for vonMises stuff and higher order
      // TODO: Manfred, it makes sense to add the micro model there - Fabian
      PiezoStressStrain<TYPE> stress_strain(mechMatSD, pde1_->GetSubTensorType(), this);
      stress_strain.CalcStressStrainResult(dynamic_cast<MechPDE*>(pde1_), res, resultType);
    }
    else
    {

      Vector<Double> intPoint;
      Vector<Double> LCoord;
      Vector<TYPE> TempE;
      Vector<TYPE> TempMechStrain;
      Vector<TYPE> TempDField;

      MechStressStrain<TYPE> * mechStressOp = NULL;

      assert(GetPde1()->GetName() == "mechanic");
      UInt stressDim = dynamic_cast<MechPDE*>(GetPde1())->GetStressStrainDim();
      UInt elecDim   = stressDim == 6 ? 3 : 2;

      // Retrieve solution from nodeStoreSolution Class
      BaseNodeStoreSol * solPDE1 = pde1_->getPDESolution();
      BaseNodeStoreSol * solPDE2 = pde2_->getPDESolution();

      NodeStoreSol<TYPE> * solhelp1 =
          dynamic_cast<NodeStoreSol<TYPE>*>(solPDE1);
      NodeStoreSol<TYPE> * solhelp2 =
          dynamic_cast<NodeStoreSol<TYPE>*>(solPDE2);

      // Determines gradient of electric potential, i.e. E=\grad \phi
      GradientFieldOp<TYPE> * FieldOp2
      = new GradientFieldOp<TYPE>(ptGrid_, pde2_, eqnMap2_,
          *solhelp2, results2_[0]->fctType,
          isaxi_);

      // Determines linear Strain S=Bu, i.e.
      //  partial derivates of mechanical displacement
      Vector<TYPE> elemElecStress, elemStressStrain, sortedStress;
      elemElecStress.Resize(stressDim);
      elemElecStress.Init(0);
      elemStressStrain.Resize(stressDim);
      elemStressStrain.Init(0);
      TempMechStrain.Resize(stressDim);
      TempMechStrain.Init();
      TempDField.Resize(elecDim);
      TempDField.Init();
      sortedStress.Resize(6);

      Matrix<TYPE> stiffnessMat, sTensor;
      Matrix<TYPE> piezoCouplingMat;
      Matrix<TYPE> piezoCouplingMatT;
      Matrix<TYPE> permittivityMat;
      Matrix<TYPE> elemDisp;

      Vector<TYPE> & actVal = actRes.GetVector();
      actVal.Resize( actRes.GetEntityList()->GetSize() * stressDim );

      //transform the type
      SubTensorType type;
      String2Enum(subType_,type);

      // get correct mechanical stress operator and piezoelectric tensor
      mechStressOp = new MechStressStrain<TYPE>(mechMatSD, type);

      // loop over all elements
      for ( it.Begin(); !it.IsEnd(); it++ ) {

        // Calc E - field;
        it.GetElem()->ptElem->GetCoordMidPoint(LCoord);

        FieldOp2->CalcElemGradField( TempE, it, LCoord, 1);

        // Calc linear mechanical stresses
        //get coordinates of element

        // compute strain
        mechMatSD = mechMat[it.GetElem()->regionId];
        mechStressOp->SetMaterial( mechMatSD );
        solhelp1->GetElemSolutionAsMatrix(elemDisp, it);
        mechStressOp->SetActElemSol(elemDisp);
        mechStressOp->SetIntPoint(LCoord);
        mechStressOp->CalcStrainVec(TempMechStrain,1,it);
        mechStressOp->UnsetIntPoint();
        elemStressStrain.Init(0);

        // isMicroModel case!
        //strain has to be reduced by irrebersibel strain
        UInt nrEl  = it.GetElem()->elemNum;
        Vector<Double> actPirr, actSirr;
        actPirr.Resize(elecDim);
        actSirr.Resize(stressDim);

        // get effective irreversible strain and polarization
        matPiezo->GetEffectiveIrreversibleValues( actPirr, actSirr, nrEl,
            false, false );

        Matrix<Double> cTensor, sTensor, dTensor, epsTensor, dTensorTrans, eTensor;
        Vector<Double> actStress(stressDim);
        Vector<Double> actE(elecDim);
        //get current tensors
        matPiezo->GetEffectiveTensors( cTensor, sTensor, epsTensor, dTensor,
            actStress, actE, nrEl,
            false, false );
        //compute actual stress
        dTensor.Transpose( dTensorTrans );
        eTensor = cTensor * dTensorTrans;

        //reduce total strain by irreversible strain
        Vector<TYPE> actStrain;
        actStrain = TempMechStrain - actSirr;

        //compute mechanical part of stress
        elemStressStrain = cTensor * actStrain;

        //electric part of stress
        TempDField = eTensor*TempE;

        //total stress
        elemStressStrain -= TempDField;
      
        for(UInt iDof = 0; iDof < stressDim; iDof++ ) {
          actVal[it.GetPos()*stressDim + iDof] = elemStressStrain[iDof];
        }
      }
      // Delete integrator again (Stressabbau ;-)
      delete mechStressOp;
      delete FieldOp2;
    }
  }


  template <class TYPE>
  void PiezoCoupling::CalcCharges( shared_ptr<BaseResult> res,
                                   RegionIdType neighbor ){

    // do cast of restult and resize solution vector
    Result<TYPE> &  actRes =
      dynamic_cast<Result<TYPE>&>(*res);
    EntityIterator it = actRes.GetEntityList()->GetIterator();

    Vector<TYPE> & actVal = actRes.GetVector();
    actVal.Resize( actRes.GetEntityList()->GetSize() );

    BaseNodeStoreSol * solPDE1 = pde1_->getPDESolution();
    BaseNodeStoreSol * solPDE2 = pde2_->getPDESolution();

    NodeStoreSol<TYPE> * solhelp1 =
      dynamic_cast<NodeStoreSol<TYPE>*>(solPDE1);
    NodeStoreSol<TYPE> * solhelp2 =
      dynamic_cast<NodeStoreSol<TYPE>*>(solPDE2);

    // Determines linear mechanical strains
     MechStressStrain<TYPE> * mechStrainOp = NULL;

     // Determines gradient of electric potential, i.e. E=\grad \phi
     GradientFieldOp<TYPE> * FieldOp2
       = new GradientFieldOp<TYPE>(ptGrid_, pde2_, eqnMap2_,
                                   *solhelp2, results2_[0]->fctType, isaxi_);


     // ------ Calculation of the electric field and linear stresses ------

     TYPE elemNormalD = 0.0;
     TYPE charge = 0.0;
     Elem * ptVolElem;
     BaseFE * ptSurfElemFE, * ptVolElemFE;

     StdVector<Elem*> elemssd;
     StdVector<SurfElem*> surfElems;
     Vector<TYPE> TempE, TempBu;
     Vector<Double> normal,lCoordSurf, lCoordVol, lCoord, Coord;
     Integer regionIndex = 0;
     Double normSign = 0.0;

     //charge operator
     ElecChargeOp<TYPE> * chargeOp;
     chargeOp = new ElecChargeOp<TYPE>(ptGrid_, pde2_, eqnMap2_,
                                       results2_[0], isaxi_);

     Matrix<TYPE> stiffnessMat;
     Matrix<TYPE> piezoCouplingMat;
     Matrix<TYPE> permittivityMat;


     // get material from mechanics
     std::map<RegionIdType, BaseMaterial*> mechMat =
       pde1_->getPDEMaterialData();

     Global::ComplexPart dataType;
    
     // get material from electrostatics
     std::map<RegionIdType, BaseMaterial*>elecMat =
       pde2_->getPDEMaterialData();

     //transform the type
     SubTensorType type;
     String2Enum(subType_,type);

     for ( it.Begin(); !it.IsEnd(); it++ ) {

       // Determine, which volume element is the right neighbour for the
       // calculation
       if ( neighbor ==
            it.GetSurfElem()->ptVolElem1->regionId ) {
         ptVolElem = it.GetSurfElem()->ptVolElem1;
         normSign = -1.0;
       } else {
         ptVolElem = it.GetSurfElem()->ptVolElem2;
         normSign = 1.0;
       }

       normSign *= static_cast<Double>(it.GetSurfElem()->normalSign);

       ptSurfElemFE = it.GetSurfElem()->ptElem;
       ptVolElemFE = ptVolElem->ptElem;
       const StdVector<UInt> & surfConnect = it.GetSurfElem()->connect;
       const StdVector<UInt> & volConnect = ptVolElem->connect;

       // calculate volume integration coordinates from
       // surfe integration coordinates for evaluating the
       // electric flux density on the surface of the volume
       // element
       ptSurfElemFE->GetCoordMidPoint(lCoordSurf);
       ptVolElemFE->GetLocalIntPoints4Surface(surfConnect, volConnect,
                                              lCoordSurf, lCoordVol);

         // Find correct material for volume element
       regionIndex = subdoms_.Find( ptVolElem->regionId );
       if ( regionIndex == -1 ) {
         EXCEPTION( "PiezoPDE::CalcCharges: The region with name '"
                    << ptGrid_->GetRegion().ToString(ptVolElem->regionId)
                    << "' of surface element no. " << it.GetSurfElem()->elemNum
                    << " is not contained in my set of regions!." );
       }

       BaseMaterial* matPiezo  = materials_[ptVolElem->regionId];
       BaseMaterial* mechMatSD = mechMat[ptVolElem->regionId];
       BaseMaterial* elecMatSD = elecMat[ptVolElem->regionId];
       if ( complexMatData_[ptVolElem->regionId] ) {
         dataType = Global::COMPLEX;
       }
       else {
         dataType = Global::REAL;
       }

       // 1.) calculate electric field
       ElemList tempList(ptGrid_);
       tempList.SetElement( ptVolElem );
       EntityIterator tempIt = tempList.GetIterator();
       FieldOp2->CalcElemGradField(TempE, tempIt, lCoordVol,1);


       // 2.) calculate linear mechanical strains
       mechStrainOp = new MechStressStrain<TYPE>(mechMatSD, type);
       matPiezo->GetTensor(piezoCouplingMat,PIEZO_TENSOR,dataType,type);
       elecMatSD->GetTensor(permittivityMat,ELEC_PERMITTIVITY,dataType,type);

       Matrix<Double> Coord;
       ptGrid_->GetElemNodesCoord( Coord, volConnect, false);

       Matrix<TYPE> elemDisp;
       ElemList elemList(ptGrid_);
       elemList.SetElement( ptVolElem );
       EntityIterator itLocal = elemList.GetIterator();

       solhelp1->GetElemSolutionAsMatrix(elemDisp, itLocal);
       mechStrainOp->SetActElemSol(elemDisp);
       mechStrainOp->SetIntPoint(lCoordVol);
       ElemList eList(ptGrid_);
       eList.SetElement( ptVolElem );
       EntityIterator it2 = eList.GetIterator();
       mechStrainOp->CalcStrainVec(TempBu,1,it2);

       Vector<TYPE> DField;
       Vector<TYPE> piezoCouplTimesStrain;

       if (subType_ == "3d"){
         DField.Resize(3);
         piezoCouplTimesStrain.Resize(3);
       }
       else {
         DField.Resize(2);
         piezoCouplTimesStrain.Resize(2);
       }

       piezoCouplTimesStrain = piezoCouplingMat * TempBu;
       DField = permittivityMat * TempE;
       DField+=piezoCouplTimesStrain;

       ptGrid_->CalcSurfNormal(normal, *it.GetSurfElem());
       normal *= normSign;
       elemNormalD = normal * DField;

       // Integrate over DField * normal
       chargeOp->CalcElemCharge(charge, it,
                                lCoordSurf, elemNormalD);

       actVal[it.GetPos()] = charge;
       delete mechStrainOp;

     }


     //  Writes result to StdPDE for later retrieval in SinglePDEs
     // (required by piezoParamIdent)
     std::string analysis =
       param->GetByVal("sequenceStep","index", sequenceStep_)
       ->Get("analysis")->GetChild()->GetName();
     if(analysis == "paramIdent") {

       // calculate sum of vector entries
       TYPE sum = 0.0;
       for( UInt i = 0; i < actVal.GetSize(); i++ ) {
         sum += actVal[i];
       }

       // Hack: SinglePDE expects a vector with size 1
       Vector<Complex> helpChargeSD(1);
       helpChargeSD[0] = sum;
       pde1_->setPDE_complexValuedCharge(helpChargeSD);
       pde2_->setPDE_complexValuedCharge(helpChargeSD);
     }

     delete FieldOp2;

     delete chargeOp;

  } // end CalcCharges


  void PiezoCoupling::CalcStrainIrr( shared_ptr<BaseResult> res ){

    UInt strainDim = dynamic_cast<MechPDE*>(GetPde1())->GetStressStrainDim();

    // Retrieve solution from nodeStoreSolution Class
    BaseNodeStoreSol * solPDE2 = pde2_->getPDESolution();

    NodeStoreSol<Double> * solElec =
      dynamic_cast<NodeStoreSol<Double>*>(solPDE2);

    // Determines gradient of electric potential, i.e. E=\grad \phi
    GradientFieldOp<Double> * ElecFieldOp
      = new GradientFieldOp<Double>( ptGrid_, pde2_, eqnMap2_,
                                     *solElec, results2_[0]->fctType,
                                     isaxi_);

    // Determines linear Strain S=Bu, i.e.
    //  partial derivates of mechanical displacement
    Vector<Double> elemStrain;
    elemStrain.Resize(strainDim);

    // get material objects
    std::map<RegionIdType, BaseMaterial*> elecMat =
      pde2_->getPDEMaterialData();
    std::map<RegionIdType, BaseMaterial*> mechMat =
      pde1_->getPDEMaterialData();

    // get result object
    Result<Double> &  actRes =
      dynamic_cast<Result<Double>&>(*res);
    EntityIterator it = actRes.GetEntityList()->GetIterator();

    Vector<Double> & actVal = actRes.GetVector();
    actVal.Resize( actRes.GetEntityList()->GetSize() * strainDim );

    // Fetch material: As we assume, that all elements belong to
    // one and the same region, we simply take the subdomain of the first
    // element
    it.Begin();

    // get the materials for the subdomain
    BaseMaterial* elecMatSD = elecMat[it.GetElem()->regionId];
    BaseMaterial* mechMatSD = mechMat[it.GetElem()->regionId];
    BaseMaterial* matPiezo = materials_[it.GetElem()->regionId];

    //transform the type
    SubTensorType type;
    String2Enum(subType_,type);

    Vector<Double> actSirr, actPirr, elecField, LCoord;
    UInt nrEl;
    std::string str;
    Double actP, Psat;
    Directions dirP;

    //check for maco-hysteresis-model
    bool isHyst = false;    
    if ( elecMatSD->getHysteresis() != NULL ) 
      isHyst = true;
    
    bool isMicroModel = false;
    if ( matPiezo->GetMicroPiezoModel() != NULL ) 
      isMicroModel = true;

    actSirr.Resize(strainDim);
    actPirr.Resize(dim_);

    // loop over all elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {

      if ( isHyst ) {
        // Calc E - field;
        it.GetElem()->ptElem->GetCoordMidPoint(LCoord);

        ElecFieldOp->CalcElemGradField( elecField, it, LCoord, 1);

        // compute actual and old polarization
        elecMatSD->GetScalar( str, P_DIRECTION );
        String2Enum( str, dirP );
        elecMatSD->GetScalar( Psat, Y_SATURATION, Global::REAL);

        nrEl  = it.GetElem()->elemNum;
        actP  = elecMatSD->ComputeScalarHystVal( nrEl, elecField[dirP] );

        // compute irreversible mechanical strain
        ComputeSirr( actSirr, type, dirP, actP, mechMatSD );

        // Calc linear mechanical stresses
        //get coordinates of element

        for(UInt iDof = 0; iDof < strainDim; iDof++ ) {
          actVal[it.GetPos()*strainDim + iDof] = actSirr[iDof];
        }
      }
      else if ( isMicroModel ) {
        nrEl  = it.GetElem()->elemNum;

        // get effective irreversible strain and polarization
        matPiezo->GetEffectiveIrreversibleValues( actPirr, actSirr, nrEl,
                                                  false, false );
        //        std::cout << "ActSirr:\n " << actSirr << std::endl;
        for(UInt iDof = 0; iDof < strainDim; iDof++ ) {
          actVal[it.GetPos()*strainDim + iDof] = actSirr[iDof];
        }
      }

      
    }
    // Delete integrator again (Stressabbau ;-)
    delete ElecFieldOp;
  }


  void PiezoCoupling::CalcElecPolarization( shared_ptr<BaseResult> res )
  {
    UInt strainDim = dynamic_cast<MechPDE*>(GetPde1())->GetStressStrainDim();

    // get result object 
    Result<Double> &  actRes = 
      dynamic_cast<Result<Double>&>(*res);
    EntityIterator it = actRes.GetEntityList()->GetIterator();
    
    Vector<Double> & actVal = actRes.GetVector();
    actVal.Resize( actRes.GetEntityList()->GetSize() * strainDim );
    
    // Fetch material: As we assume, that all elements belong to
    // one and the same region, we simply take the subdomain of the first 
    // element
    it.Begin();
    BaseMaterial* matPiezo = materials_[it.GetElem()->regionId];

    Vector<Double> actSirr, actPirr;
    UInt nrEl;

    bool isMicroModel = false;
    if ( matPiezo->GetMicroPiezoModel() != NULL ) 
      isMicroModel = true;

    actSirr.Resize(strainDim);
    actPirr.Resize(dim_);

    // loop over all elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      if ( isMicroModel ) {
        nrEl  = it.GetElem()->elemNum;
        
        // get effective irreversible strain and polarization
        matPiezo->GetEffectiveIrreversibleValues( actPirr, actSirr, nrEl,
                                                  false, false );
        for(UInt iDof = 0; iDof<dim_; iDof++ ) {
          actVal[it.GetPos()*dim_ + iDof] = actPirr[iDof];
        }
      }
     
    }
  }


  template <class TYPE>
  void PiezoCoupling::CalcElecFluxDensity( shared_ptr<BaseResult> res ){

    //    SolutionType resultType = res->GetResultInfo()->resultType;


    Vector<Double> intPoint;
    Vector<Double> LCoord;
    Vector<Double> TempE;
    Vector<Double> TempMechStrain;
    Vector<Double> TempDField;

    MechStressStrain<Double> * mechStrainOp = NULL;

    UInt strainDim = dynamic_cast<MechPDE*>(GetPde1())->GetStressStrainDim();
    UInt elecDim = strainDim == 6 ? 3 : 2;

    // Retrieve solution from nodeStoreSolution Class
    BaseNodeStoreSol * solPDE1 = pde1_->getPDESolution();
    BaseNodeStoreSol * solPDE2 = pde2_->getPDESolution();

    NodeStoreSol<Double> * solhelp1 =
      dynamic_cast<NodeStoreSol<Double>*>(solPDE1);
    NodeStoreSol<Double> * solhelp2 =
      dynamic_cast<NodeStoreSol<Double>*>(solPDE2);

    // Determines gradient of electric potential, i.e. E=\grad \phi
    GradientFieldOp<Double> * FieldOp2
      = new GradientFieldOp<Double>(ptGrid_, pde2_, eqnMap2_,
                                  *solhelp2, 
                                  results2_[0]->fctType,
                                  isaxi_);

    // Determines linear Strain S=Bu, i.e.
    //  partial derivates of mechanical displacement
    Vector<Double> elemElecStress, elemStressStrain, sortedStress;
    elemElecStress.Resize(strainDim);
    elemElecStress.Init(0);
    elemStressStrain.Resize(strainDim);
    elemStressStrain.Init(0);
    TempDField.Resize(elecDim);
    TempDField.Init();
    sortedStress.Resize(6);

    Matrix<Double> stiffnessMat, sTensor;
    Matrix<Double> piezoCouplingMat;
    Matrix<Double> piezoCouplingMatT;
    Matrix<Double> permittivityMat;
    Matrix<Double> elemDisp;



    // get material from mechanics
    std::map<RegionIdType, BaseMaterial*> mechMat =
      pde1_->getPDEMaterialData();


    // get
    Result<Double> &  actRes =
      dynamic_cast<Result<Double>&>(*res);
    EntityIterator it = actRes.GetEntityList()->GetIterator();

    Vector<Double> & actVal = actRes.GetVector();
    actVal.Resize( actRes.GetEntityList()->GetSize() * dim_ );

    // Fetch material: As we assume, that all elements belong to
    // one and the same region, we simply take the subdomain of the first
    // element
    it.Begin();

    // get the materials for the subdomain
    BaseMaterial* matPiezo = materials_[it.GetElem()->regionId];
    BaseMaterial* mechMatSD = mechMat[it.GetElem()->regionId];

    //transform the type
    SubTensorType type;
    String2Enum(subType_,type);

    // get correct mechanical strain operator 
    mechStrainOp = new MechStressStrain<Double>(mechMatSD, type);

    // bool isMicroModel = false;  // TODO: Unused variable isMicroModel
    // if ( matPiezo->GetMicroPiezoModel() != NULL ) 
    //   isMicroModel = true;
    //else
    if ( matPiezo->GetMicroPiezoModel() == NULL )
      EXCEPTION("CalcFluxDensity currently just implemented for PiezoMicroBK" );
 

    bool recompute = false;
    bool previousValues = false;
    UInt nrEl;
    Vector<Double> actSirr, actPirr;
    actSirr.Resize(strainDim);
    actPirr.Resize(dim_);

    // loop over all elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {

      // Calc E - field;
      it.GetElem()->ptElem->GetCoordMidPoint(LCoord);

      FieldOp2->CalcElemGradField( TempE, it, LCoord, 1);

      // Calc linear mechanical stresses
      // mechMatSD = mechMat[it.GetElem()->regionId];
      // mechStressOp->SetMaterial( mechMatSD );
      solhelp1->GetElemSolutionAsMatrix(elemDisp, it);
      mechStrainOp->SetActElemSol(elemDisp);
      mechStrainOp->SetIntPoint(LCoord);
      mechStrainOp->CalcStrainVec(elemStressStrain,1,it);
      mechStrainOp->UnsetIntPoint();

      // currrent finite element
      nrEl  = it.GetElem()->elemNum;
      
      // Global::ComplexPart dataType; // TODO: Unused variable dataType
      // if ( complexMatData_[it.GetElem()->regionId] ) {
      //  dataType = Global::COMPLEX;
      // }
      // else {
      //  dataType = Global::REAL;
      // }

      // get effective material tensors (currently just d-Tensor) 
      Matrix<Double> cTensor, sTensor, dTensor, epsTensor, eTensor;
      Matrix<Double> dTensorTrans;

      //get current tensors
      //epsTensor: tensor at constant mechanical stress
      Vector<Double> stress(elemStressStrain.GetSize());
      stress.Init();
      matPiezo->GetEffectiveTensors( cTensor, sTensor, epsTensor, dTensor,   
                                     stress, TempE, nrEl, 
                                     recompute, previousValues );

      matPiezo->GetEffectiveIrreversibleValues( actPirr, actSirr, nrEl,
                                                recompute, previousValues );

      //compute actual stress
      dTensor.Transpose( dTensorTrans );
      eTensor = cTensor * dTensorTrans;
      elemStressStrain -= actSirr;
      stress  = cTensor * elemStressStrain;
      stress -= eTensor * TempE;

      // compute flux denbsity
      TempDField = dTensor * stress + epsTensor * TempE + actPirr;

      for(UInt iDof = 0; iDof < dim_; iDof++ ) {
        actVal[it.GetPos()*dim_ + iDof] = TempDField[iDof];
      }

    }
    // Delete integrator again (Stressabbau ;-)
    delete mechStrainOp;
    delete FieldOp2;
  }


  // *********************
  //   DefineIntegrators
  // *********************
  void PiezoCoupling::DefineIntegrators() {

    Global::ComplexPart matType = Global::REAL;
    RegionIdType actRegion;
    // BaseMaterial * actSDMat = NULL; // TODO: Unused variable actSDMat


    // get material from electrostatics
    std::map<RegionIdType, BaseMaterial*> elecMat =
      pde2_->getPDEMaterialData();


    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      // Set current region and material
      actRegion = it->first;
      // actSDMat = it->second;
      matType = Global::REAL;

      //transform the type
      SubTensorType tensorType;
      String2Enum(subType_,tensorType);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      BiLinFormContext *actContextStiff = NULL;

      if ( nonLinPiezoMicroHF_ ) {
        // Piezoelectric problem using micro-modeling of Belov-Kreher

        // get material from mechanics
        std::map<RegionIdType, BaseMaterial*> mechMat = 
          pde1_->getPDEMaterialData();

        // get material from electrostatics
        std::map<RegionIdType, BaseMaterial*> elecMat = 
          pde2_->getPDEMaterialData();

        // get time step size
        Double dt;
        dt = dynamic_cast<TransientDriver*>(domain->GetSingleDriver())->GetDeltaT();

        //allocate for micro-modeling
        StdVector<Elem*> elemssd;
        ptGrid_->GetElems(elemssd, actRegion);
        UInt numElSD =  elemssd.GetSize();
        materials_[actRegion]->InitPiezoMicro(numElSD, actSDList,
                                              mechMat[actRegion], 
                                              elecMat[actRegion],
                                              tensorType, dt);

        // get solutions of PDEs
        BaseNodeStoreSol * solPDE1 = pde1_->getPDESolution();
        BaseNodeStoreSol * solPDE2 = pde2_->getPDESolution();
        
        NodeStoreSol<Double> * solhelp1 = 
          dynamic_cast<NodeStoreSol<Double>*>(solPDE1);    
        NodeStoreSol<Double> * solhelp2 = 
          dynamic_cast<NodeStoreSol<Double>*>(solPDE2);

        // get previous solution of PDE
        BaseNodeStoreSol * solPrevPDE1 = pde1_->getPDESolutionPrev();
        BaseNodeStoreSol * solPrevPDE2 = pde2_->getPDESolutionPrev();
        
        NodeStoreSol<Double> * solPrevhelp1 = 
          dynamic_cast<NodeStoreSol<Double>*>(solPrevPDE1);    
        NodeStoreSol<Double> * solPrevhelp2 = 
          dynamic_cast<NodeStoreSol<Double>*>(solPrevPDE2);

        //
        //------ add nonlinear coupling stiffness in mechanical equation ---
        //
        nLinMicroPiezoCouple* bilinearCoupleMech =  
          new nLinMicroPiezoCouple( materials_[actRegion], mechMat[actRegion], 
                                    elecMat[actRegion], tensorType, true);
           
        //set solution object for (n+1) and (n)
        bilinearCoupleMech->SetSolutionMech(*solhelp1);
        bilinearCoupleMech->SetSolutionElec(*solhelp2);
        bilinearCoupleMech->SetPrevSolutionMech(*solPrevhelp1);
        bilinearCoupleMech->SetPrevSolutionElec(*solPrevhelp2);

        bilinearCoupleMech->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,  
                                                results2_[0], this);

        bilinearCoupleMech->SetMatDataType( matType );
        
        BiLinFormContext *actContextCouplemech =
          new BiLinFormContext( bilinearCoupleMech, STIFFNESS  );
        
        actContextCouplemech->SetEntryType( matType );
        actContextCouplemech->SetPtPdes( pde1_, pde2_ );
        actContextCouplemech->SetResults( results1_[0], results2_[0],
                                          actSDList, actSDList );
      
        assemble_->AddBiLinearForm( actContextCouplemech );

        //
        //------  add nonlinear coupling stiffness in electric equation ---
        //
        nLinMicroPiezoCouple* bilinearCoupleElec =  
          new nLinMicroPiezoCouple( materials_[actRegion], mechMat[actRegion], 
                                    elecMat[actRegion], tensorType, false);
           
        //set solution object for (n+1) and (n)
        bilinearCoupleElec->SetSolutionMech(*solhelp1);
        bilinearCoupleElec->SetSolutionElec(*solhelp2);
        bilinearCoupleElec->SetPrevSolutionMech(*solPrevhelp1);
        bilinearCoupleElec->SetPrevSolutionElec(*solPrevhelp2);

        bilinearCoupleElec->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,  
                                               results2_[0], this);

        bilinearCoupleElec->SetMatDataType( matType );
        
        BiLinFormContext *actContextCoupleElec =
          new BiLinFormContext( bilinearCoupleElec, STIFFNESS  );
        
        actContextCoupleElec->SetEntryType( matType );
        actContextCoupleElec->SetPtPdes( pde2_, pde1_ );
        actContextCoupleElec->SetResults( results2_[0], results1_[0],
                                   actSDList, actSDList );
      
        assemble_->AddBiLinearForm( actContextCoupleElec );

        //
        //------  add nonlinear mechanical stiffness -------------------
        //
        nLinMicroPiezoMech* bilinearStiffMech =  
          new nLinMicroPiezoMech( materials_[actRegion], mechMat[actRegion], 
                                  elecMat[actRegion], tensorType);
           
        //set soltuion object for (n+1) and (n)
        bilinearStiffMech->SetSolutionMech(*solhelp1);
        bilinearStiffMech->SetSolutionElec(*solhelp2);
        bilinearStiffMech->SetPrevSolutionMech(*solPrevhelp1);
        bilinearStiffMech->SetPrevSolutionElec(*solPrevhelp2);

        bilinearStiffMech->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,  
                                              results2_[0], this);

        bilinearStiffMech->SetMatDataType( matType );
        
        BiLinFormContext *actContextStiffMech =
          new BiLinFormContext( bilinearStiffMech, STIFFNESS  );
        
        actContextStiffMech->SetEntryType( matType );
        actContextStiffMech->SetPtPdes( pde1_, pde1_ );
        actContextStiffMech->SetResults( results1_[0], results1_[0],
                                         actSDList, actSDList );
      
        assemble_->AddBiLinearForm( actContextStiffMech );

        //
        //------  add nonlinear electrostatic stiffness ---
        //
        nLinMicroPiezoElec* bilinearStiffElec =  
          new nLinMicroPiezoElec( materials_[actRegion], mechMat[actRegion], 
                                  elecMat[actRegion], tensorType);
           
        //set soltuion object for (n+1) and (n)
        bilinearStiffElec->SetSolutionMech(*solhelp1);
        bilinearStiffElec->SetSolutionElec(*solhelp2);
        bilinearStiffElec->SetPrevSolutionMech(*solPrevhelp1);
        bilinearStiffElec->SetPrevSolutionElec(*solPrevhelp2);

        bilinearStiffElec->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,  
                                              results2_[0], this);

        bilinearStiffElec->SetMatDataType( matType );
        
        BiLinFormContext *actContextStiffElec =
          new BiLinFormContext( bilinearStiffElec, STIFFNESS  );
        
        actContextStiffElec->SetEntryType( matType );
        actContextStiffElec->SetPtPdes( pde2_, pde2_ );
        actContextStiffElec->SetResults( results2_[0], results2_[0],
                                         actSDList, actSDList );
      
        assemble_->AddBiLinearForm( actContextStiffElec );

        //
        //------ micro-piezo-model: RHS for electric equation ---
        //
        MicroPiezoPolarizationElecRhsInt * elecRHS = 
          new MicroPiezoPolarizationElecRhsInt( materials_[actRegion],
                                                mechMat[actRegion],
                                                elecMat[actRegion],
                                                tensorType );	 
        
        elecRHS->SetSolutionMech(*solhelp1);
        elecRHS->SetSolutionElec(*solhelp2);
        elecRHS->SetPrevSolutionMech(*solPrevhelp1);
        elecRHS->SetPrevSolutionElec(*solPrevhelp2);
        
        elecRHS->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,  
                                    results2_[0], this);
        LinearFormContext * rhsContextElec = 
          new LinearFormContext( elecRHS );
        rhsContextElec->SetPtPde( pde2_ );
        rhsContextElec->SetResult( results2_[0], actSDList );
        assemble_->AddLinearForm( rhsContextElec );
        
        //
        //------ micro-piezo-model: RHS for mechanic equation ---
        //
        MicroPiezoPolarizationMechRhsInt * mechRHS = 
          new MicroPiezoPolarizationMechRhsInt( materials_[actRegion], 
                                                mechMat[actRegion], 
                                                elecMat[actRegion], 
                                                tensorType);
        mechRHS->SetSolutionMech(*solhelp1);
        mechRHS->SetSolutionElec(*solhelp2);
        mechRHS->SetPrevSolutionMech(*solPrevhelp1);
        mechRHS->SetPrevSolutionElec(*solPrevhelp2);
        
        // for computation of electric field
        mechRHS->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,  
                                    results2_[0], this);
        LinearFormContext * rhsContextMech = 
          new LinearFormContext( mechRHS );
        rhsContextMech->SetPtPde( pde1_ );
        rhsContextMech->SetResult( results1_[0], actSDList );
        assemble_->AddLinearForm( rhsContextMech );
      }

      else if (  nonLinHysteresis_ ) {
        //Piezoelectric problem with hysteresis

        //allocate for hystersis modeling
        StdVector<Elem*> elemssd;
        ptGrid_->GetElems(elemssd, actRegion);
        UInt numElSD =  elemssd.GetSize();
        elecMat[actRegion]->InitHyst(numElSD, actSDList);

        // get material from mechanics
        std::map<RegionIdType, BaseMaterial*> mechMat =
          pde1_->getPDEMaterialData();

        // get material from electrostatics
        std::map<RegionIdType, BaseMaterial*> elecMat =
          pde2_->getPDEMaterialData();

        // get solutions of PDEs
        BaseNodeStoreSol * solPDE1 = pde1_->getPDESolution();
        BaseNodeStoreSol * solPDE2 = pde2_->getPDESolution();

        NodeStoreSol<Double> * solhelp1 =
          dynamic_cast<NodeStoreSol<Double>*>(solPDE1);
        NodeStoreSol<Double> * solhelp2 =
          dynamic_cast<NodeStoreSol<Double>*>(solPDE2);

        // get previous solution of PDE
        BaseNodeStoreSol * solPrevPDE1 = pde1_->getPDESolutionPrev();
        BaseNodeStoreSol * solPrevPDE2 = pde2_->getPDESolutionPrev();

        NodeStoreSol<Double> * solPrevhelp1 =
          dynamic_cast<NodeStoreSol<Double>*>(solPrevPDE1);
        NodeStoreSol<Double> * solPrevhelp2 =
          dynamic_cast<NodeStoreSol<Double>*>(solPrevPDE2);

        // add nonlinear coupling stiffness in mechanical equation
        nLinPiezoHystCouple* bilinearStiffNonLin =
          new nLinPiezoHystCouple( materials_[actRegion], mechMat[actRegion],
                                   elecMat[actRegion], tensorType, true);

        //set solution object for (n+1) and (n)
        bilinearStiffNonLin->SetSolutionMech(*solhelp1);
        bilinearStiffNonLin->SetSolutionElec(*solhelp2);
        bilinearStiffNonLin->SetPrevSolutionMech(*solPrevhelp1);
        bilinearStiffNonLin->SetPrevSolutionElec(*solPrevhelp2);

        bilinearStiffNonLin->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,
                                                results2_[0], this);

        bilinearStiffNonLin->SetMatDataType( matType );

        BiLinFormContext *actContextStiffCoupling =
          new BiLinFormContext( bilinearStiffNonLin, STIFFNESS  );

        // explicit definition of counter part!!
        actContextStiffCoupling->SetCounterPart( false );

        actContextStiffCoupling->SetEntryType( matType );
        actContextStiffCoupling->SetPtPdes( pde1_, pde2_ );
        actContextStiffCoupling->SetResults( results1_[0], results2_[0],
                                   actSDList, actSDList );

        assemble_->AddBiLinearForm( actContextStiffCoupling );

        // add nonlinear coupling stiffness in electric equation
        nLinPiezoHystCouple* bilinearCoupleElec =
          new nLinPiezoHystCouple( materials_[actRegion], mechMat[actRegion],
                                   elecMat[actRegion], tensorType, false);

        //set solution object for (n+1) and (n)
        bilinearCoupleElec->SetSolutionMech(*solhelp1);
        bilinearCoupleElec->SetSolutionElec(*solhelp2);
        bilinearCoupleElec->SetPrevSolutionMech(*solPrevhelp1);
        bilinearCoupleElec->SetPrevSolutionElec(*solPrevhelp2);

        bilinearCoupleElec->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,
                                               results2_[0], this);

        bilinearCoupleElec->SetMatDataType( matType );

        BiLinFormContext *actContextCoupleElec =
          new BiLinFormContext( bilinearCoupleElec, STIFFNESS  );

        // explicit definition og counter part!!
        actContextCoupleElec->SetCounterPart( false );

        actContextCoupleElec->SetEntryType( matType );
        actContextCoupleElec->SetPtPdes( pde2_, pde1_ );
        actContextCoupleElec->SetResults( results2_[0], results1_[0],
                                   actSDList, actSDList );

        assemble_->AddBiLinearForm( actContextCoupleElec );


        // add nonlinear electrostattic stiffness
        nLinPiezoHystElec* bilinearStiffElec =
          new nLinPiezoHystElec( materials_[actRegion], mechMat[actRegion],
                                 elecMat[actRegion], tensorType);

        //set soltuion object for (n+1) and (n)
        bilinearStiffElec->SetSolutionMech(*solhelp1);
        bilinearStiffElec->SetSolutionElec(*solhelp2);
        bilinearStiffElec->SetPrevSolutionMech(*solPrevhelp1);
        bilinearStiffElec->SetPrevSolutionElec(*solPrevhelp2);

        bilinearStiffElec->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,
                                              results2_[0], this);

        bilinearStiffElec->SetMatDataType( matType );

        BiLinFormContext *actContextStiffElec =
          new BiLinFormContext( bilinearStiffElec, STIFFNESS  );

        actContextStiffElec->SetEntryType( matType );
        actContextStiffElec->SetPtPdes( pde2_, pde2_ );
        actContextStiffElec->SetResults( results2_[0], results2_[0],
                                         actSDList, actSDList );

        assemble_->AddBiLinearForm( actContextStiffElec );


        //hysteresis formulation: RHS for electric equation
        PiezoPolarizationElecRhsInt * elecRHS =
          new PiezoPolarizationElecRhsInt( materials_[actRegion],
                                           mechMat[actRegion],
                                           elecMat[actRegion],
                                           tensorType );

        elecRHS->SetSolutionMech(*solhelp1);
        elecRHS->SetSolutionElec(*solhelp2);
        elecRHS->SetPrevSolutionMech(*solPrevhelp1);
        elecRHS->SetPrevSolutionElec(*solPrevhelp2);

        elecRHS->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,
                                    results2_[0], this);
        LinearFormContext * rhsContextElec =
          new LinearFormContext( elecRHS );
        rhsContextElec->SetPtPde( pde2_ );
        rhsContextElec->SetResult( results2_[0], actSDList );
        assemble_->AddLinearForm( rhsContextElec );


        // hysteresis formulation: RHS for mechanic equation
        PiezoPolarizationMechRhsInt * mechRHS =
          new PiezoPolarizationMechRhsInt( materials_[actRegion],
                                           mechMat[actRegion],
                                           elecMat[actRegion],
                                           tensorType);
        mechRHS->SetSolutionMech(*solhelp1);
        mechRHS->SetSolutionElec(*solhelp2);
        mechRHS->SetPrevSolutionMech(*solPrevhelp1);
        mechRHS->SetPrevSolutionElec(*solPrevhelp2);

        // for computation of electric field
        mechRHS->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,
                                    results2_[0], this);
        LinearFormContext * rhsContextMech =
          new LinearFormContext( mechRHS );
        rhsContextMech->SetPtPde( pde1_ );
        rhsContextMech->SetResult( results1_[0], actSDList );
        assemble_->AddLinearForm( rhsContextMech );
      }
      else {
        // add linear stiffness
        BaseForm *bilinearStiff =
          new linPiezoCoupling(materials_[actRegion], tensorType);
        
        bilinearStiff->SetMatDataType( matType );
        
        actContextStiff =
          new BiLinFormContext( bilinearStiff, STIFFNESS  );
      

        // We also need to set the transposed of the coupling
        // matrix to the lower diagonal side
        actContextStiff->SetCounterPart( true );
        
        actContextStiff->SetEntryType( matType );
        actContextStiff->SetPtPdes( pde1_, pde2_ );
        actContextStiff->SetResults( results1_[0], results2_[0],
                                     actSDList, actSDList );
        
        assemble_->AddBiLinearForm( actContextStiff );
        
        // check for complex valued material parameter
        if( complexMatData_[actRegion] ) {
          matType = Global::IMAG;
          
          BaseForm * bilinearStiffC =
            new linPiezoCoupling( materials_[actRegion], tensorType);
          
          //GetStiffIntegrator(materialData_[actSD]);
          BiLinFormContext *actComplexContextStiff =
            new BiLinFormContext(bilinearStiffC, STIFFNESS );
          actComplexContextStiff->SetPtPdes(pde1_, pde2_);
          actComplexContextStiff->SetResults( results1_[0], results2_[0],
                                              actSDList, actSDList );
          // We also need to set the transposed of the coupling
          // matrix to the lower diagonal side
          actComplexContextStiff->SetCounterPart( true );
          
          actComplexContextStiff->SetEntryType(matType);
          bilinearStiffC->SetMatDataType(matType);
          assemble_->AddBiLinearForm( actComplexContextStiff );
        }
      }
    }

    // Define integrators for composite materials
    // (only for subType "flatShell")
    std::map<RegionIdType, Composite>::iterator compIt;
    for( compIt=compositeMaterials_.begin(); compIt!=compositeMaterials_.end();
         compIt++ ) {
      // Get current subdomain and composite material
      RegionIdType actRegion = compIt->first;
      Composite * composite = &compIt->second;

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      FlatShellPiezoInt * compPiezoInt = new FlatShellPiezoInt( composite, false);
      BiLinFormContext * stiffContext =
        new BiLinFormContext( compPiezoInt, STIFFNESS);

      // We also need to set the transposed of the coupling
      // matrix to the lower diagonal side
      stiffContext->SetCounterPart( true );

      stiffContext->SetPtPdes( pde1_, pde2_ );
      stiffContext->SetResults( results1_[0], results2_[0],
                                actSDList, actSDList );
      assemble_->AddBiLinearForm( stiffContext );

      // check for complex valued material parameter
      if( complexMatData_[actRegion] ) {
        matType = Global::IMAG;
        FlatShellPiezoInt * compPiezoIntC = new FlatShellPiezoInt( composite, false);
        compPiezoIntC->SetMatDataType(matType);
        BiLinFormContext * stiffContextC =
            new BiLinFormContext( compPiezoIntC, STIFFNESS);

        // We also need to set the transposed of the coupling
        // matrix to the lower diagonal side
        stiffContextC->SetCounterPart( true );

        stiffContextC->SetPtPdes( pde1_, pde2_ );
        stiffContextC->SetResults( results1_[0], results2_[0],
                                  actSDList, actSDList );
        stiffContextC->SetEntryType(matType);
        assemble_->AddBiLinearForm( stiffContextC );

      }
    }

  }


  void PiezoCoupling::DefineAvailResults() {

   // Check for subType
    StdVector<std::string> stressDofNames;
    if( subType_ == "3d" ) {
      stressDofNames = "xx", "yy", "zz", "yz", "xz", "xy";

    } else if( subType_ == "planeStrain" ) {
      stressDofNames = "xx", "yy", "xy";

    } else if( subType_ == "planeStress" ) {
      stressDofNames = "xx", "yy", "xy";

    } else if( subType_ == "axi" ) {
      stressDofNames = "rr", "zz", "rz", "phiphi";

    } else if( subType_ == "flatShell" ) {
      stressDofNames = "";
    }

    // Determine vectorial dofNames
    StdVector<std::string> vecDofNames;
    if( isaxi_) {
      vecDofNames = "r", "z";
    } else if (dim_ == 2) {
      vecDofNames = "x", "y";
    } else {
      vecDofNames = "x", "y", "z";
    }

    // === MECHANIC STRESS ===
    shared_ptr<ResultInfo> stress(new ResultInfo);
    stress->resultType = MECH_STRESS;
    stress->dofNames = stressDofNames;
    stress->unit = "N/m^2";
    stress->entryType = ResultInfo::TENSOR;
    stress->definedOn = ResultInfo::ELEMENT;
    stress->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( stress );

    // === MECHANIC VON MISES STRESS (yield criterion, a scalar value)===
    shared_ptr<ResultInfo> vonMises(new ResultInfo);
    vonMises->resultType = VON_MISES_STRESS;
    vonMises->dofNames = "";
    vonMises->unit =  "";
    vonMises->entryType = ResultInfo::SCALAR;
    vonMises->definedOn = ResultInfo::ELEMENT;
    vonMises->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( vonMises );

    // === MECHANIC STRAIN ===
    shared_ptr<ResultInfo> strain(new ResultInfo);
    strain->resultType = MECH_STRAIN;
    strain->dofNames = stressDofNames;
    strain->unit = "";
    strain->entryType = ResultInfo::TENSOR;
    strain->definedOn = ResultInfo::ELEMENT;
    strain->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( strain );

    // === MECHANIC STRAIN Irreversibel===
    shared_ptr<ResultInfo> strainIrr(new ResultInfo);
    strainIrr->resultType = MECH_STRAIN_IRR;
    strainIrr->dofNames = stressDofNames;
    strainIrr->unit = "";
    strainIrr->entryType = ResultInfo::TENSOR;
    strainIrr->definedOn = ResultInfo::ELEMENT;
    strainIrr->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( strainIrr );

    // === ELECTRIC POLARIZATION ===
    shared_ptr<ResultInfo> pol( new ResultInfo );
    pol->resultType = ELEC_POLARIZATION;
    pol->definedOn = ResultInfo::ELEMENT;
    pol->entryType = ResultInfo::VECTOR;
    pol->dofNames = vecDofNames;
    pol->unit = "C/m^2";
    availResults_.insert( pol );

    // === ELECTRIC Flux Density ===
    shared_ptr<ResultInfo> flux( new ResultInfo );
    flux->resultType = ELEC_FLUX_DENSITY;
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
    flux->dofNames = vecDofNames;
    flux->unit = "C/m^2";
    availResults_.insert( flux );

    // === ELECTRIC CHARGE ===
    shared_ptr<ResultInfo> charge(new ResultInfo);
    charge->resultType = ELEC_CHARGE;
    charge->dofNames = "";
    charge->unit = "C";
    charge->definedOn = ResultInfo::SURF_ELEM;
    charge->entryType = ResultInfo::SCALAR;
    charge->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( charge );
  }



  //==================================   Methods for hysteresis ============================

  void PiezoCoupling::GetNonlinMaterialTensor(Matrix<Double>& matTensor,
                                              Vector<Double>& elecD,
                                              std::string matTensorType,
                                              BaseMaterial* matMech,
                                              BaseMaterial* matElec,
                                              BaseMaterial* matCouple,
                                              SubTensorType subTensorType,
                                              Vector<Double>& elecField,
                                              Vector<Double>& mechStrain,
                                              Vector<Double>& elecFieldPrev,
                                              Vector<Double>& mechStrainPrev,
                                              EntityIterator& ent ) {

    // linear c-Tensor
    Matrix<Double> cTensor;
    matMech->GetTensor( cTensor, MECH_STIFFNESS_TENSOR, Global::REAL, subTensorType );

    // lineaer eps-Tensor at constant strain
    Matrix<Double> epsTensor_cStrain;
    matElec->GetTensor( epsTensor_cStrain, ELEC_PERMITTIVITY, Global::REAL, subTensorType );

    // linear e-Tensor
    Matrix<Double> eTensor;
    matCouple->GetTensor( eTensor, PIEZO_TENSOR, Global::REAL, subTensorType );

    // compute actual and old polarization
    UInt nrEl;
    std::string str;
    Double Psat, actP, prevP, scaleFactor;
    Directions dirP;
    matElec->GetScalar( str, P_DIRECTION );
    String2Enum( str, dirP );
    matElec->GetScalar( Psat, Y_SATURATION, Global::REAL);

    nrEl  = ent.GetElem()->elemNum;
    actP  = matElec->ComputeScalarHystVal( nrEl, elecField[dirP] );
    prevP = matElec->GetScalarHystPrevVal( nrEl );

    Double Pr; // = 0.26;
    matElec->GetScalar( Pr, Y_REMANENCE, Global::REAL);
    Psat += Pr;

    // compute d-Tensor
    Matrix<Double> sTensor, dTensor, epsTensor_cStress;
    Matrix<Double> eTensorTrans;
    cTensor.Invert(sTensor);
    dTensor = eTensor*sTensor;

    // compute eps-tensor at constant stress
    eTensor.Transpose(eTensorTrans);
    epsTensor_cStress = epsTensor_cStrain + dTensor*eTensorTrans;


    // scale the dTensor
    Matrix<Double> dsTensor, dsTensorTrans, dsPrevTensor, dsPrevTensorTrans;
    scaleFactor = (actP+Pr) / (Psat);
    dsTensor = dTensor * scaleFactor;
    dsTensor.Transpose(dsTensorTrans);

    scaleFactor = (prevP+Pr) / (Psat);
    dsPrevTensor = dTensor * scaleFactor;
    dsPrevTensor.Transpose(dsPrevTensorTrans);

    // compute irreversible mechanical strain
    Vector<Double> actSirr, prevSirr;
    ComputeSirr( actSirr, subTensorType, dirP, actP, matMech );
    ComputeSirr( prevSirr, subTensorType, dirP, prevP, matMech );
//     actSirr.Init();
//     prevSirr.Init();

//     std::cout << "\n actSirr:\n" << actSirr << std::endl;
//     std::cout << "\n prevSirr:\n" << prevSirr << std::endl;

    // compute the old mechanical stress
    Vector<Double> help, prevStressVec;
    prevStressVec  = cTensor * mechStrainPrev; 
    prevStressVec -= cTensor * prevSirr;
    help           = dsPrevTensorTrans * elecFieldPrev;
    prevStressVec -= help;

    if ( matTensorType == "CouplingTensorElecEQ" ) {
      Matrix<Double> tmp;
      tmp = dsTensor * cTensor;
      tmp.Transpose( matTensor );
    }
    else if ( matTensorType == "CouplingTensorMechEQ" ||
              matTensorType == "ElecTensor") {
      // compute differential coupling tensor
      Matrix<Double> diffdTensor, diffcTensor;
      diffdTensor = dsTensorTrans;
      ComputeDiffCouplingTensor( diffdTensor, elecField, elecFieldPrev,
                                 actSirr, prevSirr, dirP, subTensorType);

      if ( matTensorType == "CouplingTensorMechEQ" )
        matTensor = cTensor * diffdTensor;
      else {
        Matrix<Double> diffeTensor;
        diffeTensor = cTensor * diffdTensor;

        Double diffEpsVal;
        matTensor  = epsTensor_cStress;
        diffEpsVal = matElec->ComputeScalarDiffVal( nrEl, elecField[dirP] );
        //      std::cout << "diffEpsval: " << diffEpsVal << std::endl;

        if (  diffEpsVal < 0.0 )
          EXCEPTION("Negative effective permittivity");
        matTensor[dirP][dirP] += diffEpsVal;

        Matrix<Double> tmp;
        tmp = dsTensor * diffeTensor;
        matTensor -= tmp;
      }
    }
    else if ( matTensorType == "MechRHS" ) {
      scaleFactor = (actP - prevP) / Psat;
      dsTensor = dTensor * scaleFactor;
      dsTensor.Transpose(dsTensorTrans);
      matTensor = cTensor * dsTensorTrans;
    }
    else if ( matTensorType == "ElecRHS1" ) {
      Matrix<Double> dcTensor;
      dcTensor = dsTensor * cTensor;
      scaleFactor = (actP - prevP) / Psat;
      dsTensor = dTensor * scaleFactor;
      dsTensor.Transpose(dsTensorTrans);
      matTensor = dcTensor * dsTensorTrans;
    }
    else if ( matTensorType == "ElecRHS2" ) {
      scaleFactor = (actP - prevP) / Psat;
      dsTensor = dTensor * scaleFactor;
      elecD = dsTensor * prevStressVec;
    }
  }


  void PiezoCoupling::ComputeSirr( Vector<Double>& VecSirr,
                                   SubTensorType type,
                                   UInt dirP,
                                   Double actP,
                                   BaseMaterial* matMech ) {


    //get basis vector for irreversibel strain
    if ( type == FULL ) {
      VecSirr.Resize(6);
      VecSirr.Init();
      for ( UInt i=0; i<3; i++ )
	VecSirr[i] = -0.5;
      VecSirr[ dirP ] = 1.0;

//       VecSirr[0] = 1.0;
//       VecSirr[1] = -1.0;
//       VecSirr[2] = 0;

    }
    else if ( type ==  AXI ) {
      VecSirr.Resize(4);
      VecSirr.Init();
      for ( UInt i=0; i<4; i++ )
	VecSirr[i] = -0.5;
      VecSirr[2] = 0.0;  // no shear strain
      VecSirr[ dirP ] = 1.0;
    }
    else {
      VecSirr.Resize(3);
      VecSirr.Init();
      for ( UInt i=0; i<1; i++ )
	VecSirr[i] = -0.5;
      VecSirr[ dirP ] = 1.0;
    }

    // compute irreversibel strain
    Vector<Double> coeff;
    matMech->GetVector( coeff, COEFF_STRAIN_IRREVERSIBLE, Global::REAL);

    Double val =  coeff[0] +
      coeff[1] * actP +
      coeff[2] * actP * actP +
      coeff[3] * actP * actP * actP +
      coeff[4] * actP * actP * actP * actP;

    //   std::cout << "actP: " << actP << "  actVal: " << val << "\n Coeff:\n" << coeff << std::endl;
    VecSirr *= val;
  }


  void PiezoCoupling::ComputeDiffCouplingTensor( Matrix<Double>& dMat,
                                                 Vector<Double>& actE,
                                                 Vector<Double>& prevE,
                                                 Vector<Double>& actSirr,
                                                 Vector<Double>& prevSirr,
                                                 Directions dirP,
                                                 SubTensorType subTensorType ) {
    Vector<Double> diffE, diffSirr;
    diffE    = actE - prevE;
    diffSirr = actSirr - prevSirr;

    switch(subTensorType)
    {
    case PLANE_STRAIN:
      if ( abs( diffE[0]) > 1.0 ) {
        dMat[0][0] += diffSirr[0] / diffE[0];
      }
      if ( abs(diffE[1]) > 1.0 ) {
        dMat[1][1] += diffSirr[1] / diffE[1];
      }
      break;

    case AXI:
      if ( abs( diffE[0]) > 1.0 ) {
        dMat[0][0] += diffSirr[0] / diffE[0];
      }
      if ( abs(diffE[1]) > 1.0 ) {
        dMat[1][1] += diffSirr[1] / diffE[1];
      }
      break;

    case FULL:
      if ( abs( diffE[0]) > 1.0 ) {
        dMat[0][0] += diffSirr[0] / diffE[0];
      }
      if ( abs(diffE[1]) > 1.0 ) {
        dMat[1][1] += diffSirr[1] / diffE[1];
      }
      if ( abs(diffE[2]) > 1.0 ) {
        dMat[2][2] += diffSirr[2] / diffE[2];
      }
      break;

    default: EXCEPTION("Problems in ComputeDiffElasticitytensor");
    }
  }



  //==================================   Methods for Micro-Piezo-Model =========================

  void PiezoCoupling::GetMaterialTensorMicroPiezo(Matrix<Double>& matTensor, 
                                                  Vector<Double>& resultVec,
                                                  std::string matTensorType,
                                                  BaseMaterial* matMech,
                                                  BaseMaterial* matElec,
                                                  BaseMaterial* matCouple,
                                                  SubTensorType subTensorType,
                                                  Vector<Double>& elecField,
                                                  Vector<Double>& mechStrain, 
                                                  Vector<Double>& elecFieldPrev,
                                                  Vector<Double>& mechStrainPrev, 
                                                  EntityIterator& ent ) {


    bool recompute = false;
    bool previousValues = false;

    // currrent finite element level
    UInt nrEl  = ent.GetElem()->elemNum;

    Matrix<Double> cTensor, sTensor, dTensor, epsTensor, eTensor;
    Matrix<Double> dTensorTrans;

    //get current tensors
    //epsTensor: tensor at constant mechanical stress
    Vector<Double> actStress(mechStrain.GetSize()), mechStrainLin;
    actStress.Init();
    matCouple->GetEffectiveTensors( cTensor, sTensor, epsTensor, dTensor,   
                                    actStress, elecField, nrEl, 
                                    recompute, previousValues );

    //std::cout << "sTensorGet: \n" << sTensor << std::endl << std::endl;
    // get effective irreversible strain and polarization
    Vector<Double> actPirr( elecField.GetSize() );
    Vector<Double> actSirr( mechStrain.GetSize() );
    matCouple->GetEffectiveIrreversibleValues( actPirr, actSirr, nrEl,
                                               recompute, previousValues );

    //compute actual stress
    dTensor.Transpose( dTensorTrans );
    eTensor       = cTensor * dTensorTrans;
    mechStrainLin = mechStrain - actSirr;
    actStress     = cTensor * mechStrainLin;
    actStress    -= eTensor * elecField;
    //std::cout << "actStress:\n" << actStress << std::endl << std::endl;

    // recompute the tensors
    recompute = true;
    matCouple->GetEffectiveTensors( cTensor, sTensor, epsTensor, dTensor,   
                                    actStress, elecField, nrEl, 
                                    recompute, previousValues );

    // get new effective irreversible strain and polarization
    matCouple->GetEffectiveIrreversibleValues( actPirr, actSirr, nrEl,
                                              recompute, previousValues );

    // get previous Pirr and Sirr
    Vector<Double> prevPirr( elecField.GetSize() );
    Vector<Double> prevSirr( mechStrain.GetSize() );
    previousValues = true;
    matCouple->GetEffectiveIrreversibleValues( prevPirr, prevSirr, nrEl,
                                               false, previousValues );

    // compute previous stress
    Vector<Double> prevStress(mechStrain.GetSize());
    prevStress.Init();
    Matrix<Double> cTensorPrev, sTensorPrev, epsTensorPrev, dTensorPrev;
    matCouple->GetEffectiveTensors( cTensorPrev, sTensorPrev,
                                    epsTensorPrev, dTensorPrev, 
                                    actStress, elecField, nrEl, 
                                    false, previousValues );

    Matrix<Double> eTensorPrev, dTensorPrevTrans;
    dTensorPrev.Transpose( dTensorPrevTrans );
    eTensorPrev   = cTensorPrev * dTensorPrevTrans;
    mechStrainLin = mechStrainPrev - prevSirr;
    prevStress    = cTensorPrev * mechStrainLin;
    prevStress   -= eTensorPrev * elecFieldPrev;
    actStress.Init();
    prevStress.Init();

    // compute new cEffTensor and mech-dEffTensor
    Matrix<Double> dTensorEffMech;
    dTensor.Transpose( dTensorEffMech );
//     ComputeDiffCouplingTensorMicroPiezoMechEQ( dTensorEffMech, elecField, 
//                                                elecFieldPrev, actSirr, 
//                                                prevSirr, subTensorType);

    ComputeEffMechCouplingTensorMicroPiezo( sTensor, dTensorEffMech,
                                            actStress, prevStress,
                                            elecField, elecFieldPrev,
                                            actSirr, prevSirr,
                                            subTensorType );

    sTensor.Invert(cTensor); //cTensor is overwritten by the effective cTensor!!

    // compute new epsEffTensor and elec-dEffTensor
    Matrix<Double> dTensorEffElec, epsEffTensor;
    dTensor.Transpose( dTensorEffElec );
    epsEffTensor = epsTensor;
    ComputeEffElecCouplingTensorMicroPiezo( epsEffTensor, dTensorEffElec,
                                            actStress, prevStress,
                                            elecField, elecFieldPrev,
                                            actPirr, prevPirr,
                                            subTensorType );

    if ( matTensorType == "MechTensor" ) {
      matTensor = cTensor;
      //std::cout << "New cTensor:" << matTensor << std::endl;

    }
    else if ( matTensorType == "CouplingTensorElecEQ" ) {
      Matrix<Double> tmp, dTensorEffTrans;
      dTensorEffElec.Transpose( dTensorEffTrans );
      tmp = dTensorEffTrans * cTensor;
      tmp.Transpose( matTensor );
      //std::cout << "CouplingTensorElecEQ:\n" << matTensor << std::endl;

    }
    else if ( matTensorType == "CouplingTensorMechEQ" ) {
      matTensor = cTensor * dTensorEffMech;
      //        std::cout << "CouplingTensorMechEQ end:\n" << matTensor << std::endl;
    }
    else if ( matTensorType == "ElecTensor") {
      matTensor = epsEffTensor;

      Matrix<Double> eTensorEff;
      eTensorEff = cTensor * dTensorEffMech;
      
      Matrix<Double> tmp, dTensorEffTrans;
      dTensorEffElec.Transpose( dTensorEffTrans );
      tmp = dTensorEffTrans * eTensorEff;
      matTensor -= tmp;
      //   std::cout << "New EPS Tensor:\n" << matTensor << std::endl;
    }
    
    else {
      if ( matTensorType == "MechRHS1" ) {
        Matrix<Double> sTensorDelta; 
        sTensorDelta = sTensor - sTensorPrev;
        matTensor = cTensor * sTensorDelta;

        resultVec  = matTensor * prevStress; 
      }
      else if ( matTensorType == "MechRHS2" ) {
        Matrix<Double> dTensorDelta, dTensorDeltaTrans;

        dTensorDelta = dTensor - dTensorPrev;
        dTensorDelta.Transpose ( dTensorDeltaTrans );

        matTensor = cTensor * dTensorDeltaTrans;
 
      }
      else if ( matTensorType == "ElecRHS1" ) {
          Matrix<Double> dTensorDelta, dTensorDeltaTrans, dTensorEffTrans;
          Matrix<Double> epsTensorDelta;
          Matrix<Double> helpMat;
          
          dTensorDelta   = dTensor - dTensorPrev;
          dTensorDelta.Transpose(dTensorDeltaTrans);    
   
          dTensorEffElec.Transpose( dTensorEffTrans );
          helpMat   = dTensorEffTrans * cTensor;
          matTensor = helpMat * dTensorDeltaTrans;
          
          epsTensorDelta = epsTensor - epsTensorPrev;
          
          matTensor -= epsTensorDelta;

      }
      else if ( matTensorType == "ElecRHS2" ) {
        Matrix<Double> dTensorDelta, sTensorDelta, dTensorEffTrans;
        
        sTensorDelta = sTensor - sTensorPrev;
        dTensorDelta = dTensor - dTensorPrev;
        
        dTensorEffElec.Transpose( dTensorEffTrans );
        matTensor  = dTensorEffTrans * cTensor;
        matTensor *= sTensorDelta;
        
        matTensor -= dTensorDelta;
        
        resultVec  = matTensor * prevStress;
      }
    }
  }



  void PiezoCoupling::ComputeEffMechCouplingTensorMicroPiezo( Matrix<Double>& sEffMat,
                                                              Matrix<Double>& dEffMat,
                                                              Vector<Double>& actStress,
                                                              Vector<Double>& prevStress,
                                                              Vector<Double>& actE,
                                                              Vector<Double>& prevE,
                                                              Vector<Double>& actSirr,
                                                              Vector<Double>& prevSirr,
                                                              SubTensorType subTensorType ) {

    Vector<Double> diffE, diffSirr, diffStress;

    diffE      = actE - prevE;
    diffStress = actStress - prevStress;

    if ( subTensorType == FULL ) {
      UInt dim1 = 6;
      UInt dim2 = 39;
      Matrix<Double> A(dim1,dim2);
      A.Init();

      A[0][0] = diffStress[0];
      A[0][1] = diffStress[1];
      A[0][2] = diffStress[2];
      A[0][3] = diffStress[3];
      A[0][4] = diffStress[4];
      A[0][5] = diffStress[5];
 
      A[1][1]  = diffStress[0];
      A[1][6]  = diffStress[1];
      A[1][7]  = diffStress[2];
      A[1][8]  = diffStress[3];
      A[1][9]  = diffStress[4];
      A[1][10] = diffStress[5];

      A[2][2]  = diffStress[0];
      A[2][7]  = diffStress[1];
      A[2][11] = diffStress[2];
      A[2][12] = diffStress[3];
      A[2][13] = diffStress[4];
      A[2][14] = diffStress[5];

      A[3][3]  = diffStress[0];
      A[3][8]  = diffStress[1];
      A[3][12] = diffStress[2];
      A[3][15] = diffStress[3];
      A[3][16] = diffStress[4];
      A[3][17] = diffStress[5];

      A[4][4]  = diffStress[0];
      A[4][9]  = diffStress[1];
      A[4][13] = diffStress[2];
      A[4][16] = diffStress[3];
      A[4][18] = diffStress[4];
      A[4][19] = diffStress[5];

      A[5][5]  = diffStress[0];
      A[5][10] = diffStress[1];
      A[5][14] = diffStress[2];
      A[5][17] = diffStress[3];
      A[5][19] = diffStress[4];
      A[5][20] = diffStress[5];

      A[0][21] = diffE[0];
      A[0][22] = diffE[1];
      A[0][23] = diffE[2];

      A[1][24] = diffE[0];
      A[1][25] = diffE[1];
      A[1][26] = diffE[2];

      A[2][27] = diffE[0];
      A[2][28] = diffE[1];
      A[2][29] = diffE[2];

      A[3][30] = diffE[0];
      A[3][31] = diffE[1];
      A[3][32] = diffE[2];

      A[4][33] = diffE[0];
      A[4][34] = diffE[1];
      A[4][35] = diffE[2];

      A[5][36] = diffE[0];
      A[5][37] = diffE[1];
      A[5][38] = diffE[2];

      if ( A.NormL2() > 1.0 ) {
        Matrix<Double> Atrans, AA;
        //AA(dim1,dim1);
        //compute A * Atranspose
        //       for (UInt i=0; i<dim1; i++ ) {
        //         for (UInt j=0; i<dim1; i++ ) {
        //           Double val= 0;
        //           for (UInt k=0; k<dim2; k++ )
        //             val += A(i,k)*A(k,j);
        //           AA(i,j) = val;
        //         }
        //       }
        A.Transpose( Atrans );
        AA = A * Atrans;
//         std::cout << "AMech:\n" << A << std::endl << std::endl; 
//         std::cout << "AAMech:\n" << AA << std::endl << std::endl; 

//         Double maxEl = AA[0][0];
//         for (UInt i=0; i<dim1; i++ ) {
//           for (UInt j=0; j<dim1; j++) {
//             if (AA[i][j] > maxEl ) {
//                 maxEl = AA[i][j];
//             }
//           }
//         }
//         for (UInt i=0; i<dim1; i++ ) {
//           AA[i][i] += 5.0*maxEl;
//         }

        //RHS
        diffSirr   = actSirr - prevSirr;
        Matrix<Double> invAA;
        AA.Invert( invAA );
        Vector<Double> z;
        z = invAA * diffSirr;

        Vector<Double> sol;
        sol = Atrans * z;
        
        // now we have to store the result in the tensors
        sEffMat[0][0] += sol[0];
        sEffMat[0][1] += sol[1]; sEffMat[1][0] += sol[1];
        sEffMat[0][2] += sol[2]; sEffMat[2][0] += sol[2];
        sEffMat[0][3] += sol[3]; sEffMat[3][0] += sol[3];
        sEffMat[0][4] += sol[4]; sEffMat[4][0] += sol[4];
        sEffMat[0][5] += sol[5]; sEffMat[5][0] += sol[5];
        sEffMat[1][1] += sol[6];
        sEffMat[1][2] += sol[7]; sEffMat[2][1] += sol[7];
        sEffMat[1][3] += sol[8]; sEffMat[3][1] += sol[8];
        sEffMat[1][4] += sol[9]; sEffMat[4][1] += sol[9];
        sEffMat[1][5] += sol[10]; sEffMat[5][1] += sol[10];
        sEffMat[2][2] += sol[11];
        sEffMat[2][3] += sol[12]; sEffMat[3][2] += sol[12];
        sEffMat[2][4] += sol[13]; sEffMat[4][2] += sol[13];
        sEffMat[2][5] += sol[14]; sEffMat[5][2] += sol[14];
        sEffMat[3][3] += sol[15];
        sEffMat[3][4] += sol[16]; sEffMat[4][3] += sol[16];
        sEffMat[3][5] += sol[17]; sEffMat[5][3] += sol[17];
        sEffMat[4][4] += sol[18]; 
        sEffMat[4][5] += sol[19]; sEffMat[5][4] += sol[19];
        sEffMat[5][5] += sol[20];
//         std::cout << "sEff:\n" << sEffMat << std::endl;
        
//        std::cout << "dMat:\n" << dEffMat << std::endl;
        dEffMat[0][0] += sol[21];
        dEffMat[1][0] += sol[22];
        dEffMat[2][0] += sol[23];
        dEffMat[0][1] += sol[24];
        dEffMat[1][1] += sol[25];
        dEffMat[2][1] += sol[26];
        dEffMat[0][2] += sol[27];
        dEffMat[1][2] += sol[28];
        dEffMat[2][2] += sol[29];
        dEffMat[0][3] += sol[30];
        dEffMat[1][3] += sol[31];
        dEffMat[2][3] += sol[32];
        dEffMat[0][4] += sol[33];
        dEffMat[1][4] += sol[34];
        dEffMat[2][4] += sol[35];
        dEffMat[0][5] += sol[36];
        dEffMat[1][5] += sol[37];
        dEffMat[2][5] += sol[38];
        //        std::cout << "dEffMech:\n" << dEffMat << std::endl << std::endl;
      }
    }
    else 
      EXCEPTION( "Just 3D implementation of ComputeElecCouplingIrrTensorMicroPiezo");



  }

  void PiezoCoupling:: ComputeEffElecCouplingTensorMicroPiezo( Matrix<Double>& epsEffMat,
                                                               Matrix<Double>& dEffMat,
                                                               Vector<Double>& actStress,
                                                               Vector<Double>& prevStress,
                                                               Vector<Double>& actE,
                                                               Vector<Double>& prevE,
                                                               Vector<Double>& actPirr,
                                                               Vector<Double>& prevPirr,
                                                               SubTensorType subTensorType ) {


    Vector<Double> diffE, diffPirr, diffStress;

    diffE      = actE - prevE;
    diffStress = actStress - prevStress;

    if ( subTensorType == FULL ) {
      UInt dim1 = 3;
      UInt dim2 = 24;
      Matrix<Double> A(dim1,dim2);
      A.Init();

      A[0][0] = diffE[0];
      A[0][1] = diffE[1];
      A[0][2] = diffE[2];

      A[1][1] = diffE[0];
      A[1][3] = diffE[1];
      A[1][4] = diffE[2];

      A[2][2] = diffE[0];
      A[2][4] = diffE[1];
      A[2][5] = diffE[2];

      A[0][6]  = diffStress[0];
      A[0][7]  = diffStress[1];
      A[0][8]  = diffStress[2];
      A[0][9]  = diffStress[3];
      A[0][10] = diffStress[4];
      A[0][11] = diffStress[5];
 
      A[1][12] = diffStress[0];
      A[1][13] = diffStress[1];
      A[1][14] = diffStress[2];
      A[1][15] = diffStress[3];
      A[1][16] = diffStress[4];
      A[1][17] = diffStress[5];

      A[2][18] = diffStress[0];
      A[2][19] = diffStress[1];
      A[2][20] = diffStress[2];
      A[2][21] = diffStress[3];
      A[2][22] = diffStress[4];
      A[2][23] = diffStress[5];

      if ( A.NormL2() > 1.0e2 ) {
        Matrix<Double> AA;
        AA.Resize(dim1,dim1);

        //compute A * Atranspose
        for (UInt i=0; i<dim1; i++ ) {
          for (UInt j=0; j<dim1; j++ ) {
            Double val= 0;
            for (UInt k=0; k<dim2; k++ )
              val += A[i][k]*A[j][k];
            AA[i][j] = val;
          }
        }
        
        //some regularization
        Double maxEl = AA[0][0];
        for (UInt i=0; i<dim1; i++ ) {
          for (UInt j=0; j<dim1; j++) {
            if (AA[i][j] > maxEl ) {
                maxEl = AA[i][j];
            }
          }
        }
        for (UInt i=0; i<dim1; i++ ) {
          AA[i][i] += maxEl;
        }

        //RHS
        diffPirr   = actPirr - prevPirr;

        //Solve
//         Vector<Double> z(dim1);
//         AA.DirectSolve( z, diffPirr );
        
        Matrix<Double> invAA;
        AA.Invert( invAA );
        Vector<Double> z; 
        z = invAA * diffPirr;
        
        Vector<Double> sol;
        sol.Resize(dim2);
        for ( UInt i=0; i<dim2; i++) {
          Double val = 0;
          for ( UInt j=0; j<dim1; j++) {
            val += A[j][i] * z[j]; 
          }
          sol[i] = val;
        }

        // now we have to store the result in the tensors
        //        std::cout << "epsMat:\n" << epsEffMat << std::endl << std::endl;
        epsEffMat[0][0] += sol[0];
        epsEffMat[0][1] += sol[1];  epsEffMat[1][0] += sol[1];
        epsEffMat[0][2] += sol[2];  epsEffMat[2][0] += sol[2];
        epsEffMat[1][1] += sol[3];
        epsEffMat[1][2] += sol[4];  epsEffMat[2][1] += sol[4];
        epsEffMat[2][2] += sol[5];
        //        std::cout << "epsEffMat:\n" << epsEffMat << std::endl << std::endl;
                
        //std::cout << "dMat:\n" << dEffMat << std::endl << std::endl;
        Matrix<Double> dOldMat = dEffMat;
        dEffMat[0][0] += sol[6];
        dEffMat[0][1] += sol[7];
        dEffMat[0][2] += sol[8];
        dEffMat[0][3] += sol[9];
        dEffMat[0][4] += sol[10];
        dEffMat[0][5] += sol[11];
        dEffMat[1][0] += sol[12];
        dEffMat[1][1] += sol[13];
        dEffMat[1][2] += sol[14];
        dEffMat[1][3] += sol[15];
        dEffMat[1][4] += sol[16];
        dEffMat[1][5] += sol[17];
        dEffMat[2][0] += sol[18];
        dEffMat[2][1] += sol[19];
        dEffMat[2][2] += sol[20];
        dEffMat[2][3] += sol[21];
        dEffMat[2][4] += sol[22];
        dEffMat[2][5] += sol[23];
        //std::cout << "dEffElec:\n" << dEffMat << std::endl << std::endl;
      }
    }
    else 
      EXCEPTION( "Just 3D implementation of ComputeElecCouplingIrrTensorMicroPiezo");

 
  }


  void PiezoCoupling::ComputeDiffCouplingTensorMicroPiezoMechEQ( Matrix<Double>& dMat, 
                                                                 Vector<Double>& actE,
                                                                 Vector<Double>& prevE,
                                                                 Vector<Double>& actSirr,
                                                                 Vector<Double>& prevSirr,
                                                                 SubTensorType subTensorType ) {

    Vector<Double> diffE, diffSirr;

    diffE    = actE - prevE;
    diffSirr = actSirr - prevSirr;

//     std::cout << "diffE:\n" << diffE << std::endl;
//     std::cout << "diffS:\n" << diffSirr << std::endl << std::endl;

    UInt idx1Max, idx2Max;
    if ( subTensorType == PLANE_STRAIN ) {
      idx1Max = 2;
      idx2Max = 3;
    }
    else if ( subTensorType == FULL ) {
      idx1Max = 3;
      idx2Max = 6;
    }
    else 
      EXCEPTION( "Problems in ComputeDiffechTensorMicroPiezo");

    //do computation
    for (UInt i=0; i<idx1Max; i++) {
      if ( abs( diffE[i]) > 1.0e3 ) {
        for (UInt j=0; j<idx2Max; j++) {
          dMat[i][j] += diffSirr[j] / diffE[i];
        }
      }
    }



    //    std::cout << "dMatEff:\n" << dMat << std::endl;

  }

  void PiezoCoupling::ComputeDiffCouplingTensorMicroPiezoElecEQ( Matrix<Double>& dMat, 
                                                                 Vector<Double>& actStress,
                                                                 Vector<Double>& prevStress,
                                                                 Vector<Double>& actPirr,
                                                                 Vector<Double>& prevPirr,
                                                                 SubTensorType subTensorType ) {

    Vector<Double> diffStress, diffPirr;

    diffStress = actStress - prevStress;
    diffPirr   = actPirr - prevPirr;

    UInt idx1Max, idx2Max;
    if ( subTensorType == PLANE_STRAIN ) {
      idx1Max = 2;
      idx2Max = 3;
    }
    else if ( subTensorType == FULL ) {
      idx1Max = 3;
      idx2Max = 6;
    }
    else 
      EXCEPTION( "Problems in ComputeDiffechTensorMicroPiezo");

    //do computation
    for (UInt i=0; i<idx2Max; i++) {
      if ( abs( diffStress[i]) > 1.0e3 ) {
        for (UInt j=0; j<idx1Max; j++) {
          dMat[j][i] += diffPirr[j] / diffStress[i];
        }
      }
    }

  }


  void PiezoCoupling::ComputeDiffMechTensorMicroPiezo( Matrix<Double>& dMat, 
                                                       Vector<Double>& actStress,
                                                       Vector<Double>& prevStress,
                                                       Vector<Double>& actSirr,
                                                       Vector<Double>& prevSirr,
                                                       SubTensorType subTensorType ) {

    Vector<Double> diffStress, diffSirr;

    diffStress = actStress - prevStress;
    diffSirr   = actSirr - prevSirr;

    std::cout << "StressPrev:\n" << prevStress << std::endl;
    std::cout << "StressAct:\n" << actStress << std::endl << std::endl;

    UInt idxMax;
    if ( subTensorType == PLANE_STRAIN ) 
      idxMax = 3;
    else if ( subTensorType == FULL ) 
      idxMax = 6;
    else 
      EXCEPTION( "Problems in ComputeDiffMechTensorMicroPiezo");

    //do computation
    for (UInt i=0; i<idxMax; i++) {
      if ( abs(diffStress[i]) > 1.0e3 ) {
        for (UInt j=0; j<idxMax; j++) {
          dMat[j][i] += diffSirr[j] / diffStress[i];
        }
      }
    }
  }

  void PiezoCoupling::ComputeDiffElecTensorMicroPiezo( Matrix<Double>& dMat, 
                                                       Vector<Double>& actElec,
                                                       Vector<Double>& prevElec,
                                                       Vector<Double>& actPirr,
                                                       Vector<Double>& prevPirr,
                                                       SubTensorType subTensorType ) {

    Vector<Double> diffElec, diffPirr;

    diffElec   = actElec - prevElec;
    diffPirr   = actPirr - prevPirr;

    UInt idxMax;
    if ( subTensorType == PLANE_STRAIN ) 
      idxMax = 2;
    else if ( subTensorType == FULL ) 
      idxMax = 3;
    else 
      EXCEPTION( "Problems in ComputeDiffechTensorMicroPiezo");

    //do computation
    for (UInt i=0; i<idxMax; i++) {
      if ( abs( diffElec[i]) > 1.0e3 ) {
        for (UInt j=0; j<idxMax; j++) {
          dMat[j][i] += diffPirr[j] / diffElec[i];
        }
      }
    }
  }

  void PiezoCoupling::GetMaterialTensorMicroPiezo2(Matrix<Double>& matTensor, 
                                                  Vector<Double>& resultVec,
                                                  std::string matTensorType,
                                                  BaseMaterial* matMech,
                                                  BaseMaterial* matElec,
                                                  BaseMaterial* matCouple,
                                                  SubTensorType subTensorType,
                                                  Vector<Double>& elecField,
                                                  Vector<Double>& mechStrain, 
                                                  Vector<Double>& elecFieldPrev,
                                                  Vector<Double>& mechStrainPrev, 
                                                  EntityIterator& ent ) {
    
    bool recompute = false;
    bool previousValues = false;
    UInt nrEl;

    // currrent finite element
    nrEl  = ent.GetElem()->elemNum;

    // get effective material tensors 
    Matrix<Double> cTensor, sTensor, dTensor, epsTensor, eTensor;
    Matrix<Double> dTensorTrans;
    Vector<Double> stress;
    matCouple->GetEffectiveTensors( cTensor, sTensor, epsTensor, dTensor,   
                                    stress, elecField, nrEl, 
                                    recompute, previousValues );

    // get effective irreversible strain and polarization
    Vector<Double> actPirr( elecField.GetSize() );
    Vector<Double> actSirr( mechStrain.GetSize() );
    matCouple->GetEffectiveIrreversibleValues( actPirr, actSirr, nrEl,
                                               recompute, previousValues );

    // compute mechanical stress
    dTensor.Transpose( dTensorTrans ); 
    eTensor = cTensor * dTensorTrans;
    stress  = cTensor * mechStrain;
    stress  -= cTensor * actSirr;
    stress -= eTensor * elecField;

    std::cout << "actStress:\n " << stress << std::endl;

    // compute with new stress the volume factions and the corresponding
    // new material tensors
    recompute = true;
    matCouple->GetEffectiveTensors( cTensor, sTensor, epsTensor, dTensor,   
                                    stress, elecField, nrEl, 
                                    recompute, previousValues );

    // now get new irreversibel values
    matCouple->GetEffectiveIrreversibleValues( actPirr, actSirr, nrEl,
                                               recompute, previousValues );

    if ( matTensorType == "MechTensor" ) {
      matTensor = cTensor;
      std::cout << "newCtensor:\n" << cTensor << std::endl;
    }
    else if ( matTensorType == "CouplingTensorElecEQ" ) {
      Matrix<Double> tmp;
      tmp = dTensor * cTensor;
      tmp.Transpose( matTensor );
      //      std::cout << "CouplingTensorElecEQ:\n" << matTensor << std::endl;
    }
    else if ( matTensorType == "CouplingTensorMechEQ" || 
              matTensorType == "ElecTensor") {

      Vector<Double> prevPirr( elecField.GetSize() );
      Vector<Double> prevSirr( mechStrain.GetSize() );
      previousValues = true;
      matCouple->GetEffectiveIrreversibleValues( prevPirr, prevSirr, nrEl,
                                                 false, previousValues );

      // compute differential coupling tensor
      Matrix<Double> diffdTensor, diffcTensor;
      dTensor.Transpose( diffdTensor );
      ComputeDiffCouplingTensorMicroPiezoMechEQ( diffdTensor, elecField, 
                                                 elecFieldPrev, actSirr, 
                                                 prevSirr, subTensorType);
      
      if ( matTensorType == "CouplingTensorMechEQ" ) {
        matTensor = cTensor * diffdTensor;
        //std::cout << "CouplingTensorMechEQ:\n" << matTensor << std::endl;

      }
      else {
        Matrix<Double> diffeTensor;
        diffeTensor = cTensor * diffdTensor;

        // Double diffEpsVal; unused
        matTensor  = epsTensor;
        Double diffE, diffP;
        for ( UInt i=0; i< elecField.GetSize(); i++) {
          diffE = elecField[i] - elecFieldPrev[i];
          diffP = actPirr[i] - prevPirr[i];
          if ( abs(diffE) > 1.0 ) {
            matTensor[i][i] += diffP/diffE;
          }
        }
         
        Matrix<Double> tmp;
        tmp = dTensor * diffeTensor;
        matTensor -= tmp;
      }
    }

    else {
      // compute previous stress
      Vector<Double> stressPrev;

      Vector<Double> prevPirr( elecField.GetSize() );
      Vector<Double> prevSirr( mechStrain.GetSize() );
      previousValues = true;
      matCouple->GetEffectiveIrreversibleValues( prevPirr, prevSirr, nrEl,
                                                 false, previousValues );

      recompute = false;
      Matrix<Double> cTensorPrev, sTensorPrev, epsTensorPrev, dTensorPrev;
      matCouple->GetEffectiveTensors( cTensorPrev, sTensorPrev,
                                      epsTensorPrev, dTensorPrev, stress, 
                                      elecField, nrEl, 
                                      recompute, previousValues );

      // compute mechanical stress
      Matrix<Double> eTensorPrev, dTensorPrevTrans;
      dTensorPrev.Transpose( dTensorPrevTrans );
      eTensorPrev = cTensorPrev * dTensorPrevTrans;
      stressPrev  = cTensorPrev * mechStrainPrev;
      stressPrev -= cTensorPrev * prevSirr;
      stressPrev -= eTensorPrev * elecFieldPrev;

      if ( matTensorType == "MechRHS" ) {
        Matrix<Double> sTensorDelta, sTensorPrev;
        Matrix<Double> dTensorDelta, dTensorDeltaTrans, helpMat1;

        cTensor.Invert( sTensorDelta );
        cTensorPrev.Invert( sTensorPrev );
        sTensorDelta -=  sTensorPrev;
        helpMat1 = cTensor * sTensorDelta;
        helpMat1 *= -1.0;
        resultVec = helpMat1 * stressPrev;

        dTensorDelta = dTensor - dTensorPrev;
        dTensorDelta.Transpose ( dTensorDeltaTrans );

        helpMat1 = cTensor * dTensorDeltaTrans;
        resultVec -= helpMat1 * elecFieldPrev;
      }
      else if ( matTensorType == "ElecRHS" ) {
        Matrix<Double> sTensorDelta, sTensorPrev;
        Matrix<Double> dTensorDelta, epsTensorDelta;
        Matrix<Double> dTensorDeltaTrans, helpMat1, helpMat2;

        cTensor.Invert( sTensorDelta );
        cTensorPrev.Invert( sTensorPrev );
        sTensorDelta -=  sTensorPrev;

        dTensorDelta   = dTensor - dTensorPrev;
        epsTensorDelta = epsTensor - epsTensorPrev; 

        helpMat1 = dTensor * cTensor;
        helpMat2 = helpMat1 * sTensorDelta;
        helpMat1 = -helpMat2;
        helpMat1 += dTensorDelta;
        resultVec = helpMat1 * stressPrev;

        dTensorDelta.Transpose ( dTensorDeltaTrans );
        helpMat1 = dTensor * cTensor;
        helpMat2 = helpMat1 * dTensorDeltaTrans;
        helpMat1 = -helpMat2;
        helpMat1 += epsTensorDelta;
        resultVec += helpMat1 * elecFieldPrev;
      }
    }
  }

}


