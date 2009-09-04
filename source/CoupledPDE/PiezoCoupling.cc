// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "PiezoCoupling.hh"


#include "Driver/assemble.hh"
#include "Materials/baseMaterial.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/WriteInfo.hh"

// integrator (bi-)linear forms
#include "Forms/mechStressStrain.hh"
#include "Forms/gradfieldop.hh"
#include "Forms/elecchargeop.hh"
#include "Forms/nLinPiezoHyst.hh"
#include "Forms/nLinPiezoHyst.hh"
#include "Forms/nLinPiezoHystRHS.hh"
#include "Forms/nLinPiezoCoupling.hh"
#include "Forms/linElecInt.hh"
#include "Forms/FlatShellPiezoInt.hh"
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"

#include "Utils/elemstoresol.hh"

#include "PDE/SinglePDE.hh"
#include "PDE/StdPDE.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  PiezoCoupling::PiezoCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                ParamNode * paramNode  )
    : BasePairCoupling( pde1, pde2, paramNode ) {


    couplingName_ = "piezoDirect";
    materialClass_ = PIEZO;

    // determine subtype from mechanic pde
    pde1_->GetParamNode()->Get( "subType", subType_ );

    // read nonlinearity
    nonLinPiezoCoupling_=false;
    ParamNode * nonLinNode = myParam_->Get("nonLinList", false );
    if( nonLinNode ) {
       StdVector<ParamNode*> nonLinNodes = nonLinNode->GetChildren();
       for( UInt i = 0; i < nonLinNodes.GetSize(); i++ ) {

         std::string nonLinearity = nonLinNodes[i]->GetName();
         if( nonLinearity == "material" || nonLinearity == "hysteresis"
             || nonLinearity == "geo" ) {
           nonLinPiezoCoupling_=true;
           ReadPiezoNonLin();
           break;
         }
       }
    }


    // check for complex material data
    hasComplexMatParams_ = false;
    ParamNode * matNode = myParam_->Get("materialDataType", false);
    if( matNode )
      hasComplexMatParams_ =
        matNode->Get("type")->AsString() == "imagMaterialParameter";

    if( hasComplexMatParams_ ) {
      isComplex_ = true;
    }

  }


  // **************
  //   Destructor
  // **************
  PiezoCoupling::~PiezoCoupling() {
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
      if ( isComplex_ ) {
        CalcStressStrain<Complex>( result );
      } else {
        CalcStressStrain<Double>( result );

      }
      break;

    case MECH_STRAIN:
      if ( isComplex_ ) {
        CalcStressStrain<Complex>( result );
      } else {
        CalcStressStrain<Double>( result );

      }
      break;

    case MECH_STRAIN_IRR:
      if ( isComplex_ ) {
        EXCEPTION("Ireversible Strain makes no sense in Harmonic analysis");
      } else {
        CalcStrainIrr( result );

      }
      break;

    default:
      Warning( "Resulttype not computable by piezoelectric coupling",
               __FILE__, __LINE__ );
    }
  }

  template <class TYPE>
  void PiezoCoupling::CalcStressStrain( shared_ptr<BaseResult> res ){

    SolutionType resultType = res->GetResultInfo()->resultType;

    UInt stressDim=0, elecDim=0;
    Vector<Double> intPoint;
    Vector<Double> LCoord;
    Vector<TYPE> TempE;
    Vector<TYPE> TempMechStress;
    Vector<TYPE> TempDField;

    MechStressStrain<TYPE> * mechStressOp = NULL;
    if (subType_ == "planeStrain") {
      stressDim = 3;
      elecDim   = 2;
    }

    else if (subType_ == "axi") {
      stressDim = 4;
      elecDim   = 2;
    }

    else if (subType_ == "3d") {
      stressDim = 6;
      elecDim   = 3;
    }

    else
      EXCEPTION("StressOp: Unknown subtype in mech PDE!" );

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
                                  *solhelp2, ELEC_POTENTIAL,
                                  results2_[0],
                                  isaxi_);

    // Determines linear Strain S=Bu, i.e.
    //  partial derivates of mechanical displacement
    Vector<TYPE> elemElecStress, elemStressStrain, sortedStress;
    elemElecStress.Resize(stressDim);
    elemElecStress.Init(0);
    elemStressStrain.Resize(stressDim);
    elemStressStrain.Init(0);
    TempMechStress.Resize(stressDim);
    TempMechStress.Init();
    TempDField.Resize(elecDim);
    TempDField.Init();
    sortedStress.Resize(6);

    Matrix<TYPE> stiffnessMat, sTensor;
    Matrix<TYPE> piezoCouplingMat;
    Matrix<TYPE> piezoCouplingMatT;
    Matrix<TYPE> permittivityMat;
    Matrix<TYPE> elemDisp;

    Global::ComplexPart dataType;
    if ( hasComplexMatParams_ ) {
      dataType = Global::COMPLEX;
    }
    else {
      dataType = Global::REAL;
    }

    // get material from mechanics
    std::map<RegionIdType, BaseMaterial*> mechMat =
      pde1_->getPDEMaterialData();


    // get
    Result<TYPE> &  actRes =
      dynamic_cast<Result<TYPE>&>(*res);
    EntityIterator it = actRes.GetEntityList()->GetIterator();

    Vector<TYPE> & actVal = actRes.GetVector();
    actVal.Resize( actRes.GetEntityList()->GetSize() * stressDim );

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

    // get correct mechanical stress operator and piezoelectric tensor
    mechStressOp = new MechStressStrain<TYPE>(mechMatSD, type);


    // loop over all elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {

      // Calc E - field;
      it.GetElem()->ptElem->GetCoordMidPoint(LCoord);

      FieldOp2->CalcElemGradField( TempE, it, LCoord, 1);

      // Calc linear mechanical stresses
      //get coordinates of element

      // c^E S
      mechMatSD = mechMat[it.GetElem()->regionId];
      mechStressOp->SetMaterial( mechMatSD );
      solhelp1->GetElemSolutionAsMatrix(elemDisp, it);
      mechStressOp->SetActElemSol(elemDisp);
      mechStressOp->SetIntPoint(LCoord);
      mechStressOp->CalcStressVec(TempMechStress,1,it);
      mechStressOp->UnsetIntPoint();
      elemStressStrain.Init(0);

      // get correct coupling tensor and transpose it
      matPiezo = materials_[it.GetElem()->regionId];
      matPiezo->GetTensor(piezoCouplingMat,PIEZO_TENSOR,dataType,type);
      piezoCouplingMat.Transpose(piezoCouplingMatT);

      if ( resultType  == MECH_STRAIN ) {
        mechMatSD = mechMat[it.GetElem()->regionId];
        mechMatSD->GetTensor( stiffnessMat,  MECH_STIFFNESS_TENSOR,
                              dataType, type);
        stiffnessMat.Invert( sTensor );
        elemStressStrain  = sTensor * TempMechStress;
      }
      else {
        TempDField = piezoCouplingMatT*TempE;
        elemStressStrain = TempMechStress-TempDField;
      }

      for(UInt iDof = 0; iDof < stressDim; iDof++ ) {
        actVal[it.GetPos()*stressDim + iDof] = elemStressStrain[iDof];
      }

    }
    // Delete integrator again (Stressabbau ;-)
    delete mechStressOp;
    delete FieldOp2;
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
                                   *solhelp2, ELEC_POTENTIAL,
                                   results2_[0], isaxi_);


     // ------ Calculation of the electric field and linear stresses ------

     TYPE elemNormalD = 0.0;
     TYPE charge = 0.0;
     Elem * ptVolElem;
     BaseFE * ptSurfElemFE, * ptVolElemFE;
     SurfElem * ptSurfElem = NULL;

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
     if ( hasComplexMatParams_ ) {
       dataType = Global::COMPLEX;
     }
     else {
       dataType = Global::REAL;
     }

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
         EXCEPTION( "PiezoPDE:CalcCharges The region with Name "
                    << ptGrid_->RegionIdToName(ptVolElem->regionId)
                    << " of surface element Nr. " << ptSurfElem->elemNum
                    << "is not contained in my set of regions!." );
       }

       BaseMaterial* matPiezo  = materials_[ptVolElem->regionId];
       BaseMaterial* mechMatSD = mechMat[ptVolElem->regionId];
       BaseMaterial* elecMatSD = elecMat[ptVolElem->regionId];

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
       param->Get("sequenceStep", "index", GenStr(sequenceStep_) )
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

    UInt strainDim=0;
    if (subType_ == "planeStrain") {
      strainDim = 3;
    }

    else if (subType_ == "axi") {
      strainDim = 4;
    }

    else if (subType_ == "3d") {
      strainDim = 6;
    }

    else
      EXCEPTION("StressOp: Unknown subtype in mech PDE!" );

    // Retrieve solution from nodeStoreSolution Class
    BaseNodeStoreSol * solPDE2 = pde2_->getPDESolution();

    NodeStoreSol<Double> * solElec =
      dynamic_cast<NodeStoreSol<Double>*>(solPDE2);

    // Determines gradient of electric potential, i.e. E=\grad \phi
    GradientFieldOp<Double> * ElecFieldOp
      = new GradientFieldOp<Double>( ptGrid_, pde2_, eqnMap2_,
                                     *solElec, ELEC_POTENTIAL,
                                     results2_[0],
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

    //transform the type
    SubTensorType type;
    String2Enum(subType_,type);

    Vector<Double> actSirr, elecField, LCoord;
    UInt nrEl;
    std::string str;
    Double actP, Psat;
    Directions dirP;

    bool isHyst = false;
    if ( elecMatSD->getHysteresis() != NULL ) {
      isHyst = true;
    }

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

    }
    // Delete integrator again (Stressabbau ;-)
    delete ElecFieldOp;
  }


  void PiezoCoupling::ReadPiezoNonLin(){


    // Check, if "nonLinList" is present
    ParamNode * nonLinListNode = myParam_->Get("nonLinList", false );


    // --------------------------------------
    // @Tom:
    // Shouldn't it be possible to define various types of non-linearities
    // for different regions like in the other PDEs (e.g. acoustic PDE)?
    // Currently the non-linearity type is only determined by the first entry
    // --------------------------------------


    if( nonLinListNode) {
      nonLin_ = false;

      // Get nonlinear types
      StdVector<ParamNode*> nonLinNodes = nonLinListNode->GetChildren();
      for( UInt i = 0; i < nonLinNodes.GetSize(); i++ ) {

        std::string actTypeString = nonLinNodes[i]->GetName();
        std::string actId = nonLinNodes[i]->Get("id")->AsString();

        NonLinType actType;
        String2Enum( actTypeString, actType );
        nonLinIdType_[actId] = actType;

        // check type
        if( actType == HYSTERESIS ) {
          nonLin_ = true;
          nonLinHysteresis_ = true;
        }

        if( actType == MATERIAL ) {
          nonLin_ = true;
          nonLinMaterial_ = true;
        }

        if( actType == GEOMETRIC ) {
          nonLin_ = true;
        }
      }
    }

    // Run over all region and set entry in "regionNonLinId"
    StdVector<ParamNode*> regionNodes =
      myParam_->Get("regionList")->GetChildren();

    RegionIdType actRegionId;
    std::string actRegionName, actNonLinId;

    if( regionNodes.GetSize() > 0 ) {
      Info->PrintF( couplingName_, "Non-linearity in following region(s)\n" );
    }
    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {

      // get data
      regionNodes[i]->Get( "name", actRegionName );
      regionNodes[i]->Get( "nonLinId", actNonLinId );

      if( actNonLinId == "" )
        continue;

      actRegionId = ptGrid_->RegionNameToId( actRegionName );

      // Check nonLinId was already registerd
      if( nonLinIdType_.find( actNonLinId) == nonLinIdType_.end() ) {
        EXCEPTION( "NonLinearity with id '" << actNonLinId
                   << "' was not defined in 'nonLinList'" );
      }

      regionNonLinId_[actRegionId] = actNonLinId;
      regionNonLinType_[actRegionId] = nonLinIdType_[actNonLinId];

      // Log to info file
      std::string nonLinString;
      Enum2String( nonLinIdType_[actNonLinId], nonLinString );
      Info->PrintF( couplingName_, " %s: %s\n", actRegionName.c_str(),
                    nonLinString.c_str() );

    }


    // Check, if nonlinear type is the same for all regions
    if( regionNonLinType_.size() > 1 ) {
      std::map<RegionIdType, NonLinType>::iterator it;
      it = regionNonLinType_.begin();
      NonLinType firstType = it->second;
      for( ; it != regionNonLinType_.end(); it++ ) {
        if( it->second != firstType ) {
          EXCEPTION( "Non-linearity should be the same for all regions!" );
        }
      }
    }


    if(nonLin_) {
      if ( !nonLinHysteresis_ ) {
        pde1_->SetNonLinearity(nonLin_);
        pde2_->SetNonLinearity(nonLin_);
        //    this->SetNonLinearity(nonLin_);

        if ( nonLinMaterial_ ) {
          pde1_->SetMaterialNonLinearity(nonLinMaterial_);
          pde2_->SetMaterialNonLinearity(nonLinMaterial_);
          //      this->SetMaterialNonLinearity(nonLin_);
        }
      }
    }
  }


  // *********************
  //   DefineIntegrators
  // *********************
  void PiezoCoupling::DefineIntegrators() {


    Global::ComplexPart matType = Global::REAL;
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;


    // get material from electrostatics
    std::map<RegionIdType, BaseMaterial*> elecMat =
      pde2_->getPDEMaterialData();


    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;
      matType = Global::REAL;

      //transform the type
      SubTensorType tensorType;
      String2Enum(subType_,tensorType);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      BiLinFormContext *actContextStiff = NULL;

      if (  nonLinHysteresis_ ) {
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
        if (  nonLinMaterial_ ) {
          //material nonlinearity

          ApproxData *nlinFnc;
          std::string nlfnc = materials_[actRegion]->GetNonlinFileName();
          materials_[actRegion]->GetScalar(nlfnc,NONLIN_DATA_NAME);

          nlinFnc = new SmoothSpline(nlfnc);
          // ApproxData *nlinFnc = new SmoothSpline("PZT4_Coupl33NL.dat");
          // oder        ApproxData *nlinFnc = new LinInterpolate(nlfnc);
          nlinFnc->CalcBestParameter();
          nlinFnc->CalcApproximation();

          // get material from mechanics
          std::map<RegionIdType, BaseMaterial*> mechMat =
            pde1_->getPDEMaterialData();

          // get material from electrostatics
          std::map<RegionIdType, BaseMaterial*> elecMat =
            pde2_->getPDEMaterialData();

          // add nonlinear stiffness
          nLinPiezoCoupling * bilinearStiffNonLin;
          bilinearStiffNonLin =  new nLinPiezoCoupling(nlinFnc, materials_[actRegion],
                                                       mechMat[actRegion],
                                                       elecMat[actRegion],
                                                       tensorType);

          BaseNodeStoreSol * solPDE1 = pde1_->getPDESolution();
          BaseNodeStoreSol * solPDE2 = pde2_->getPDESolution();

          NodeStoreSol<Double> * solhelp1 =
            dynamic_cast<NodeStoreSol<Double>*>(solPDE1);
          NodeStoreSol<Double> * solhelp2 =
            dynamic_cast<NodeStoreSol<Double>*>(solPDE2);

          bilinearStiffNonLin->SetSolution1(*solhelp1);
          bilinearStiffNonLin->SetSolution2(*solhelp2);

          bilinearStiffNonLin->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,
                                                  results2_[0]);

          bilinearStiffNonLin->SetMatDataType( matType );

          actContextStiff =
            new BiLinFormContext( bilinearStiffNonLin, STIFFNESS  );

          // We also need to set the transposed of the coupling
          // matrix to the lower diagonal side
          actContextStiff->SetCounterPart( true );
        }
        else {
          // add linear stiffness
          BaseForm *bilinearStiff =
            new linPiezoCoupling(materials_[actRegion], tensorType);

          bilinearStiff->SetMatDataType( matType );

          actContextStiff =
            new BiLinFormContext( bilinearStiff, STIFFNESS  );
        }


        // We also need to set the transposed of the coupling
        // matrix to the lower diagonal side
        actContextStiff->SetCounterPart( true );

        actContextStiff->SetEntryType( matType );
        actContextStiff->SetPtPdes( pde1_, pde2_ );
        actContextStiff->SetResults( results1_[0], results2_[0],
                                     actSDList, actSDList );

        assemble_->AddBiLinearForm( actContextStiff );

        // check for complex valued material parameter
        if( hasComplexMatParams_ ) {
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

    // === MECHANIC STRESS ===
    shared_ptr<ResultInfo> stress(new ResultInfo);
    stress->resultType = MECH_STRESS;
    stress->dofNames = stressDofNames;
    stress->unit = "N/m^2";
    stress->entryType = ResultInfo::TENSOR;
    stress->definedOn = ResultInfo::ELEMENT;
    stress->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( stress );

    // === MECHANIC STRAIN ===
    shared_ptr<ResultInfo> strain(new ResultInfo);
    strain->resultType = MECH_STRAIN;
    strain->dofNames = stressDofNames;
    strain->unit = "N/m^2";
    strain->entryType = ResultInfo::TENSOR;
    strain->definedOn = ResultInfo::ELEMENT;
    strain->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( strain );

    // === MECHANIC STRAIN Irreversibel===
    shared_ptr<ResultInfo> strainIrr(new ResultInfo);
    strainIrr->resultType = MECH_STRAIN_IRR;
    strainIrr->dofNames = stressDofNames;
    strainIrr->unit = "N/m^2";
    strainIrr->entryType = ResultInfo::TENSOR;
    strainIrr->definedOn = ResultInfo::ELEMENT;
    strainIrr->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( strainIrr );

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
    help  = mechStrainPrev;
    help -= dsPrevTensorTrans * elecFieldPrev;
    help -= prevSirr;
    prevStressVec = cTensor * help;

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
      if ( abs( diffE[0]) > 1 ) {
        if ( dirP == X )
          dMat[dirP][dirP] += diffSirr[0] / diffE[0];
      }
      if ( abs(diffE[1]) > 1 ) {
        if ( dirP == Y )
          dMat[dirP][dirP] += diffSirr[1] / diffE[1];
      }
      break;

    case FULL:
      if ( abs( diffE[0]) > 1 ) {
        if ( dirP == X )
          dMat[dirP][dirP] += diffSirr[0] / diffE[0];
      }
      if ( abs(diffE[1]) > 1 ) {
        if ( dirP == Y )
          dMat[dirP][dirP] += diffSirr[1] / diffE[1];
      }
      if ( abs(diffE[2]) > 1 ) {
        if ( dirP == Z )
          dMat[dirP][dirP] += diffSirr[2] / diffE[2];
      }
      break;

    default: EXCEPTION("Problems in ComputeDiffElasticitytensor");
    }
  }

}
