#include "PiezoCoupling.hh"


#include "Driver/assemble.hh"
#include "Materials/baseMaterial.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/writeresults.hh"

// integrator (bi-)linear forms
#include "Forms/linPiezoCoupling.hh"
#include "Forms/forms_header.hh"
#include "Forms/linElecInt.hh"
#include "Forms/FlatShellPiezoInt.hh"

#include "Utils/elemstoresol.hh"

#include "PDE/SinglePDE.hh"
#include "PDE/StdPDE.hh"


namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  PiezoCoupling::PiezoCoupling( SinglePDE *pde1, SinglePDE *pde2 )
    : BasePairCoupling( pde1, pde2 ) {

    ENTER_FCN( "PiezoCoupling::PiezoCoupling" );

    couplingName_ = "piezoDirect";
    materialClass_ = PIEZO;
    params->Get( "subType", subType_, "mechanic");

    if (subType_ == "axi") 
      isaxi_ = true;

  }


  // **************
  //   Destructor
  // **************
  PiezoCoupling::~PiezoCoupling() {
    ENTER_FCN( "PiezoCoupling::~PiezoCoupling" );
  }


  // ***************
  //   PostProcess
  // ***************
  void PiezoCoupling::PostProcess() {
    ENTER_FCN( "PiezoCoupling::PostProcess" );

    //calc charges
    if (calcCharge_.GetSize() !=0 ) {
      
      if (analysisType_ == HARMONIC || analysisType_== MULTIHARMONIC){
        CalcCharges<Complex>();
      }
      else{
        CalcCharges<Double>();
      }
    }

    // calc stresses
    if (calcStress_.GetSize() !=0 ) {
      if (analysisType_ == HARMONIC || analysisType_==MULTIHARMONIC)
        CalcStress<Complex>();
      else
        CalcStress<Double>();
    }
  }
  
  template <class TYPE>
  void PiezoCoupling::CalcStress(){
    ENTER_FCN("PiezoCoupling::CalcStress");

    ShortInt  stressDim=0, elecDim=0;
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
      Info->Error("StressOp: Unknown subtype in mech PDE! "
                  ,__FILE__,__LINE__);  

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
 
 
    Vector<TYPE> elemElecStress, elemStress, sortedStress;
    elemElecStress.Resize(stressDim);
    elemStress.Resize(stressDim);
    TempMechStress.Resize(stressDim);
    TempMechStress.Init();
    TempDField.Resize(elecDim);
    TempDField.Init();
    elemElecStress.Init(0);
    elemStress.Init(0);
    sortedStress.Resize(6);
    UInt pdeElem;

    Matrix<TYPE> stiffnessMat;
    Matrix<TYPE> piezoCouplingMat; 
    Matrix<TYPE> piezoCouplingMatT; 
    Matrix<TYPE> permittivityMat;
    Matrix<TYPE> elemDisp;
 
    // check for complex material data
    bool isComplex = params->HasValue( "type", "imagMaterialParameter", 
                                       "materialDataType" );
    DataType dataType;
    if ( isComplex ) {
      dataType = COMPLEX;
    }
    else {
      dataType = REAL;
    }

    // get material from mechanics
    std::map<RegionIdType, BaseMaterial*> mechMat =    
      pde1_->getPDEMaterialData();


    // loop over all subdomains
    for (UInt isd=0; isd<calcStress_.GetSize(); isd++) {

//       Integer regionIndex = -1;
//       regionIndex = subdoms_.Find( calcStress_[isd] );
//       if ( regionIndex == -1 ) {
//         (*error) << "PiezoPDE:CalcStress:: For the region with Name " 
//                  << ptGrid_->RegionIdToName(calcStress_[isd])
//                  << " no material data was found.!";
//         Error( __FILE__, __LINE__ );
//       }

      // get the materials for the subdomain
      BaseMaterial* matPiezo = materials_[calcStress_[isd]];
      BaseMaterial* mechMatSD = mechMat[calcStress_[isd]];

      //transform the type
      SubTensorType type;
      String2Enum(subType_,type);

      // get correct mechanical stress operator and piezoelectric tensor
      mechStressOp = new MechStressStrain<TYPE>(mechMatSD, type);
      matPiezo->GetTensor(piezoCouplingMat,PIEZO_TENSOR,dataType,type);

      //transpose the piezoelectric tensor
      piezoCouplingMat.Transpose(piezoCouplingMatT);
       
      ElemList actSDList(ptGrid_ );
      actSDList.SetRegion( subdoms_[isd] );
   
      EntityIterator it = actSDList.GetIterator();
      for ( it.Begin(); !it.IsEnd(); it++ ) {
     
        // Calc E - field;
        it.GetElem()->ptElem->GetCoordMidPoint(LCoord);

        FieldOp2->CalcElemGradField( TempE, it, LCoord, 1); 

        pdeElem = eqnMap1_->Mesh2PdeElem(it.GetElem()->elemNum);

        // Calc linear mechanical stresses

        //get coordinates of element

        // c^E S
        solhelp1->GetElemSolutionAsMatrix(elemDisp, it);
        mechStressOp->SetActElemSol(elemDisp);
        mechStressOp->SetIntPoint(LCoord);
        mechStressOp->CalcStressVec(TempMechStress,1,it);

        elemStress.Init(0);
        TempDField = piezoCouplingMatT*TempE;
        elemStress = TempMechStress-TempDField;

        pde1_->sortStresses(elemStress,sortedStress);
        stress_->SetElemResult(pdeElem-1, sortedStress);

      }
      // Delete integrator again (Stressabbau ;-)
      delete mechStressOp;
    }
 
    delete FieldOp2;

    
  
  }

  template <class TYPE>
  void PiezoCoupling::CalcCharges(){
    ENTER_FCN("PiezoCoupling::CalcCharges");

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
    UInt pdeElemNum = 0;
    Elem * ptVolElem;
    BaseFE * ptSurfElemFE, * ptVolElemFE;
    SurfElem * ptSurfElem;

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

    Vector<TYPE> chargeSD(calcCharge_.GetSize());
    chargeSD.Init();

    Matrix<TYPE> stiffnessMat;
    Matrix<TYPE> piezoCouplingMat; 
    Matrix<TYPE> permittivityMat;

    // check for complex material data
    bool isComplex = params->HasValue( "type", "imagMaterialParameter", 
                                       "materialDataType" );


    // get material from mechanics
    std::map<RegionIdType, BaseMaterial*> mechMat = 
      pde1_->getPDEMaterialData();

    DataType dataType;
    if ( isComplex ) {
      dataType = COMPLEX;
    }
    else {
      dataType = REAL;
    }

    // get material from electrostatics
    std::map<RegionIdType, BaseMaterial*>elecMat = 
      pde2_->getPDEMaterialData();

    //transform the type
    SubTensorType type;
    String2Enum(subType_,type);

    // loop over all subdomains
    for (UInt isd=0; isd<calcCharge_.GetSize(); isd++) {

      SurfElemList actSDList( ptGrid_ );
      actSDList.SetRegion( calcCharge_[isd] );
      EntityIterator it = actSDList.GetIterator();
      for ( it.Begin(); !it.IsEnd(); it++ ) {
    
        // Determine, which volume element is the right neighbour for the 
        // calculation
        if ( chargeNeighborRegion_.
             Find(it.GetSurfElem()->ptVolElem1->regionId) != -1 ) {
          ptVolElem = it.GetSurfElem()->ptVolElem1;
          normSign = -1.0;
        } else {
          ptVolElem = it.GetSurfElem()->ptVolElem2;
          normSign = 1.0;
        }
    
        normSign *= (double) it.GetSurfElem()->normalSign;

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
          (*error) << "PiezoPDE:CalcCharges:: The region with Name " 
                   << ptGrid_->RegionIdToName(ptVolElem->regionId)
                   << " of surface element Nr. " << ptSurfElem->elemNum
                   << "is not contained in my set of regions!.";
          Error( __FILE__, __LINE__ );
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
    
        pdeElemNum = eqnMap2_->Mesh2PdeElem(ptVolElem->elemNum);
    
        Vector<TYPE> chargeVec(1);
        chargeVec[0] = charge;
        charges_->SetElemResult(pdeElemNum-1, chargeVec);

        chargeSD[isd] += charge;        
        delete mechStrainOp;
      }
  
    }


    //  Writes result to StdPDE for later retrieval in SinglePDEs 
    // (required by piezoParamIdent)
    std::string analysis;
    params->Get( "type", analysis, "analysis" );
    if(analysis == "paramIdent") {  

      Vector<Complex> helpChargeSD(chargeSD.GetSize());
      for (UInt iC=0;iC<chargeSD.GetSize();iC++)
        helpChargeSD[iC] = (Complex) chargeSD[iC];
      pde1_->setPDE_complexValuedCharge(helpChargeSD);
      pde2_->setPDE_complexValuedCharge(helpChargeSD);
    }    


    Info->PrintF(couplingName_, "Computed surface charge: ");
    Info->PrintVec(chargeSD);

    delete FieldOp2;
    
  } // end CalcCharges




  // *********************
  //   DefineIntegrators
  // *********************
  void PiezoCoupling::DefineIntegrators() {

    ENTER_FCN( "PiezoCoupling::DefineIntegrators" );
    
    DataType matType = REAL;
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;
      matType = REAL;
      
      //transform the type
      SubTensorType tensorType;
      String2Enum(subType_,tensorType);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      // add stiffness
      BaseForm *bilinearStiff = 
        new linPiezoCoupling(materials_[actRegion], tensorType);

      //GetStiffIntegrator( materials_[actSD] );

      BiLinFormContext *actContextStiff =
        new BiLinFormContext( bilinearStiff, STIFFNESS  );

      bilinearStiff->SetMatDataType( matType );

      actContextStiff->SetMatDataType( matType );
      actContextStiff->SetPtPdes( pde1_, pde2_ );
      actContextStiff->SetResults( results1_[0], results2_[0],
                                   actSDList, actSDList );
      assemble_->AddBiLinearForm( actContextStiff );

      // check for complex valued material parameter
      if( params->HasValue( "type", "imagMaterialParameter", 
                            "materialDataType" ) ) {
        matType = IMAG; 

        BaseForm * bilinearStiffC = 
          new linPiezoCoupling( materials_[actRegion], tensorType);

	//GetStiffIntegrator(materialData_[actSD]);
        BiLinFormContext *actComplexContextStiff = 
          new BiLinFormContext(bilinearStiffC, STIFFNESS );
        actComplexContextStiff->SetPtPdes(pde1_, pde2_);
        actComplexContextStiff->SetResults( results1_[0], results2_[0],
                                            actSDList, actSDList );
        actComplexContextStiff->SetMatDataType(matType);
        bilinearStiffC->SetMatDataType(matType);
        assemble_->AddBiLinearForm( actComplexContextStiff );
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
      
      FlatShellPiezoInt * compPiezoInt = new FlatShellPiezoInt( composite );
      BiLinFormContext * stiffContext = 
        new BiLinFormContext( compPiezoInt, STIFFNESS);
      stiffContext->SetPtPdes( pde1_, pde2_ );
      stiffContext->SetResults( results1_[0], results2_[0],
                                actSDList, actSDList );
      assemble_->AddBiLinearForm( stiffContext );
      
    }
    
  }


  void PiezoCoupling::WriteResultsInFile( const UInt kstep,
                                          const Double asteptime,
                                          UInt stepOffset,
                                          Double timeOffset ) {

    ENTER_FCN( "PiezoCoupling::WriteResultsInFile" );


    ComplexFormat complexFormat_ = AMPLITUDE_PHASE; // or REAL_IMAG

    Double actTime = asteptime + timeOffset;
    UInt actStep = kstep + stepOffset;

    //!< pointer to output file
    WriteResults * outFile = pde1_->getPDE_outFile();
        
    if (analysisType_ == STATIC || analysisType_==TRANSIENT) {

      if (calcStress_.GetSize() !=0 ) {
        ElemStoreSol<Double> & stressConverted = 
          dynamic_cast<ElemStoreSol<Double>&>(*stress_);
        outFile->WriteElemSolutionTransient(stressConverted, 
                                            actStep, actTime);
      }
     
      //element results
      if (calcCharge_.GetSize() !=0 ) {
        ElemStoreSol<Double> & chargesConverted = 
          dynamic_cast<ElemStoreSol<Double>&>(*charges_);
        outFile->WriteElemSolutionTransient(chargesConverted, 
                                            actStep, actTime);
      } 
    }
   
    else if (analysisType_ == HARMONIC || analysisType_==MULTIHARMONIC) {
     
      //element results
      if (calcCharge_.GetSize() !=0 ) {
        ElemStoreSol<Complex> & chargesConverted = 
          dynamic_cast<ElemStoreSol<Complex>&>(*charges_);
        outFile->WriteElemSolutionHarmonic(chargesConverted, actStep, 
                                           actTime, complexFormat_);
      }
      if (calcStress_.GetSize() !=0 ) {
        ElemStoreSol<Complex> & stressConverted = 
          dynamic_cast<ElemStoreSol<Complex>&>(*stress_);
        outFile->WriteElemSolutionHarmonic(stressConverted, actStep,
                                           actTime, complexFormat_);
      }
     
    }
  }
  
  // ********************
  //   ReadStoreResults
  // ********************

  //retrieve information on desired output resulta from xml-file 
  void PiezoCoupling::ReadStoreResults() {
    ENTER_FCN( "PiezoCoupling::ReadStoreResults" );

    // Construct vectors for restricted parameter search
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    StdVector<std::string> regionNames;
    std::string quantity;
    
    // *****************************
    // Determine nodal results
    // ***************************** 
    StdVector<std::string> nodeValues;
    keyVec  = couplingName_, "storeResults", "elemResults", "region";
    attrVec = "", "", "type";  
  
    // --- Mechanic Stress ---
    Enum2String(MECH_STRESS, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, regionNames );
    ptGrid_->RegionNameToId ( calcStress_, regionNames );
    
    // If the symbolic name is "all" compute stresses for all regions
    if ( calcStress_.GetSize() == 1 && calcStress_[0] == ALL_REGIONS ) {
      calcStress_ = subdoms_;
    }
    
    //    Log to info file
    if ( calcStress_.GetSize() > 0 ) {
      hasOutput_ = true;
      Info->PrintF(couplingName_,
                   "Computing mechanical stress for regions:\n" );
      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( couplingName_, "%s\n", regionNames[k].c_str() );
      }
      Info->PrintF( "", "\n" );

      if( analysisType_ == HARMONIC || analysisType_==MULTIHARMONIC) {
        stress_ = new ElemStoreSol<Complex>;
      } else {
        stress_ = new ElemStoreSol<Double>;
      }
             
      // Resize solution arrays
      stress_->SetNumSolutions(1);
      stress_->SetSolutionType(MECH_STRESS);
      stress_->SetNumElems(pde1_->getPDE_numElems());
            
      // We always store for six components 
      // (unverg-file-format as capa does)
      stress_->SetNumDofs(6);
            
      // if output is gmv only the three normalstresses are stored
      if (params->HasValue( "format", "gmv", "output"))
        stress_->SetNumDofs(3);
            
      stress_->SetPtrEQNData(pde1_->GetEqnMap().get(), ptGrid_);
      stress_->Init();
            
    }

    // --- Electric Charges ---
    // check for charge computation
    params->GetList( "region", regionNames, couplingName_, "charge" );
    ptGrid_->RegionNameToId( chargeNeighborRegion_, regionNames );

    params->GetList( "element", regionNames, couplingName_, "charge" );
    ptGrid_->RegionNameToId( calcCharge_, regionNames );

    
    if ( calcCharge_.GetSize() > 0 ) {
      
      hasOutput_ = true;
      Info->PrintF( couplingName_, 
                    "Computing electric charges for regions:\n" );
      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( couplingName_, "%s\n", regionNames[k].c_str() );
      }
      Info->PrintF( "", "\n" );
      
      // Resize solution arrays
      if( analysisType_ == HARMONIC || analysisType_==MULTIHARMONIC) {
        charges_ = new ElemStoreSol<Complex>;
      } else {
        charges_ = new ElemStoreSol<Double>;
      }
      charges_->SetNumSolutions(1);
      charges_->SetSolutionType(ELEC_CHARGE);
      charges_->SetNumElems(pde2_->getPDE_numElems());
      charges_->SetNumDofs(1);
      charges_->SetPtrEQNData(pde2_->GetEqnMap().get(), ptGrid_);
      charges_->Init();
    } 

    //  if ( calcDField_.GetSize() > 0 ) {
    //         hasOutput_ = true;
    //         Info->PrintF( couplingName_, "Computing electric field for regions:\n" );

    //         for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
    //           Info->PrintF(couplingName_, "%s\n", regionNames[k].c_str() );
    //         }
    //         Info->PrintF( "", "\n" );

    //         std::cout<<"analysisType_"<<std::endl;
    //         std::cout<<analysisType_<<std::endl;
    //         if ( analysisType_ == HARMONIC || analysisType_==MULTIHARMONIC) {
    //           std::cout<<" Resize solution arrays ... harmonic case"<<std::endl;

    //           DfieldComplex_.SetNumSolutions(1);
    //           DfieldComplex_.SetSolutionType(ELEC_FIELD_INTENSITY);
    //           std::cout<<"pde2_->getPDE_numElems()"<<std::endl;
    //           std::cout<<pde2_->getPDE_numElems()<<std::endl;
    //           DfieldComplex_.SetNumElems(pde2_->getPDE_numElems());
    //           DfieldComplex_.SetNumDofs(pde2_->getPDE_spaceDim());
    //           DfieldComplex_.SetPtrEQNData(eqnData_, ptGrid_);
    //           DfieldComplex_.Init(); 
    //         }
    //         else {
    //           std::cout<<" Resize solution arrays ... static case"<<std::endl;
    //           Dfield_.SetNumSolutions(1);
    //           Dfield_.SetSolutionType(ELEC_FIELD_INTENSITY);
    //           Dfield_.SetNumElems(pde2_->getPDE_numElems());
    //           std::cout<<"dim_:"<<std::endl;
    //           std::cout<<pde2_->getPDE_spaceDim()<<std::endl;
    //           Dfield_.SetNumDofs(pde2_->getPDE_spaceDim());
    //           Dfield_.SetPtrEQNData(eqnData_, ptGrid_);
    //           Dfield_.Init(); 
    //         }
    //       }

    // *****************************
    // Determine element history
    // *****************************
    StdVector<std::string> saveElemHist;
    keyVec  = couplingName_, "storeResults", "elemHistory", "saveElems";
    attrVec = "", "", "";
    valVec = "", "", "";
    params->GetList(keyVec, attrVec, valVec, saveElemHist);

    if (saveElemHist.GetSize() > 0) {
      std::string errMsg = couplingName_;
      errMsg += ": Saving history elements is not implemented yet!\n";
      errMsg += "Meanwhile you can use 'unvtool' to extract element data.";
      Error( errMsg.c_str(), __FILE__, __LINE__);
    }

  }


}
