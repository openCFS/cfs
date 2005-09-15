#include "PiezoCoupling.hh"


#include "Driver/assemble.hh"
#include "DataInOut/MaterialData.hh"
#include "DataInOut/LoadMaterialData.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

// integrator (bi-)linear forms
#include "Forms/linPiezoCoupling3D.hh"
#include "Forms/linPiezoCoupling2DAxi.hh"
#include "Forms/linPiezoCoupling2DPlaneStrain.hh"
#include "Forms/forms_header.hh"
#include "Forms/linElecInt2D.hh"
#include "Forms/linElecInt3D.hh"

#include "Utils/elemstoresol.hh"

#include "PDE/SinglePDE.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  PiezoCoupling::PiezoCoupling( SinglePDE *pde1, SinglePDE *pde2 )
    : BasePairCoupling( pde1, pde2 ) {

    ENTER_FCN( "PiezoCoupling::PiezoCoupling" );

    eqnData_ = pde1-> getPDE_eqnData(); 

    couplingName_ = "piezoDirect";
    materialClass_ = "piezo";
    params->Get( "subType", subType_, "mechanic");
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

    std::string analysis;
    params->Get( "type", analysis, "analysis" );
    String2Enum(analysis,analysistype_);

    if (calcCharge_.GetSize() !=0 ) {
      
      if (analysistype_ == HARMONIC || analysistype_==MULTIHARMONIC){
        CalcComplexValuedCharges();
      }
      else{
        CalcCharges();
      }
    }

    
    if (calcDfield_.GetSize() !=0 ) {
      if (analysistype_==HARMONIC || analysistype_==MULTIHARMONIC){
        std::cout<<"PiezoCoupling::PostProces DField -- needs further implementation! "<<std::endl;
        //        CalcComplexValuedDfie();
      }
      else{
        std::cout<<"PiezoCoupling::PostProces CalcDField -- needs further implementation! "<<std::endl;
        //        CalcDfield();
      }
    }
  }


  void PiezoCoupling::CalcCharges(){
    ENTER_FCN("PiezoCoupling::CalcCharges");

    BaseNodeStoreSol * solPDE1 = pde1_->getPDESolution();
    BaseNodeStoreSol * solPDE2 = pde2_->getPDESolution();
    eqnData_ = pde2_-> getPDE_eqnData(); 

    NodeStoreSol<Double> * solhelp1 = dynamic_cast<NodeStoreSol<Double>*>(solPDE1);    
    NodeStoreSol<Double> * solhelp2 = dynamic_cast<NodeStoreSol<Double>*>(solPDE2);
    
    // Determines gradient of electric potential, i.e. E=\grad \phi
    GradientFieldOp<Double> * FieldOp2
      = new GradientFieldOp<Double>(ptGrid_, pde2_,pde2_->getPDE_eqnData(),
                                    *solhelp2, ELEC_POTENTIAL, 
                                    isaxi_);

    // Determines linear Strain S=Bu, i.e.
    //  partial derivates of mechanical displacement

    LinStrainOp<Double> * FieldOp1 
      = new LinStrainOp<Double>(ptGrid_, pde1_,pde1_->getPDE_eqnData(),
                                *solhelp1,MECH_DISPLACEMENT, 
                                isaxi_);

    // ------ Calculation of the electric field and linear Strain ------

    Double elemNormalD = 0.0;
    Double charge = 0.0;
    UInt pdeElemNum = 0;
    Elem * ptVolElem;
    BaseFE * ptSurfElemFE, * ptVolElemFE;
    SurfElem * ptSurfElem;
   
    StdVector<Elem*> elemssd;
    StdVector<SurfElem*> surfElems;
    Vector<Double> TempE, TempBu, normal,lCoordSurf, lCoordVol, lCoord, Coord;
    Integer regionIndex = 0;
    UInt pdeElem;
    Double normSign = 0.0;

    //charge operator  
    ElecChargeOp * chargeOp;
    chargeOp = new ElecChargeOp(ptGrid_, pde2_, pde2_->getPDE_eqnData(), isaxi_);

    Vector<Double> chargeSD(calcCharge_.GetSize());
    chargeSD.Init(0);

    Matrix<Double> stiffnessMat;
    Matrix<Double> piezoCouplingMat; 
    Matrix<Double> permittivityMat;

    if (subType_ == "axi") 
      isaxi_=TRUE;

    Matrix<Double> *matDat=pde1_->getPDEMaterialData()->GetMatrix();
    calcMaterialMatrices(stiffnessMat, piezoCouplingMat, permittivityMat, matDat);

    // loop over all subdomains
    for (UInt isd=0; isd<calcCharge_.GetSize(); isd++)
      {
        // get surface and according volume elements
        ptGrid_->GetSurfElems( surfElems, calcCharge_[isd] );

        // loop over all surface elements
        for (UInt iel=0; iel<surfElems.GetSize(); iel++){

          // Determine, which volume element is the right neighbour for the 
          // calculation
          if ( chargeNeighborRegion_.
               Find(surfElems[iel]->ptVolElem1->regionId) != -1 ) {
            ptVolElem = surfElems[iel]->ptVolElem1;
            normSign = -1.0;
          } else {
            ptVolElem = surfElems[iel]->ptVolElem2;
            normSign = 1.0;
          }
          
          normSign *= (double) surfElems[iel]->normalSign;

          ptSurfElemFE = surfElems[iel]->ptElem;
          ptVolElemFE = ptVolElem->ptElem;
          const StdVector<UInt> & surfConnect = surfElems[iel]->connect;
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
        
          // calculate electric field
          FieldOp2->CalcElemGradField(TempE, ptVolElem, lCoordVol,1);

          StdVector<UInt> connecth = ptVolElem->connect;
          //get coordinates of element
          Matrix<Double> Coord;
          pde1_->GetElemCoords(connecth, Coord);

          FieldOp1->CalcElemLinearStrain(TempBu, ptVolElem, Coord);

          Vector<Double> DField;
          Vector<Double> piezoCouplTimesStrain;

          if (subType_ == "3d"){
            DField.Resize(3);
            piezoCouplTimesStrain.Resize(3);
          }
          else
            {
              DField.Resize(2);
              piezoCouplTimesStrain.Resize(2);
            }

          piezoCouplingMat.Mult(TempBu,piezoCouplTimesStrain);
          permittivityMat.Mult(TempE,DField);
          DField+=piezoCouplTimesStrain;

          //pdeElem = pde2_-> getPDE_eqnData()->Mesh2PDEElem(elemssd[iel]->elemNum);
            
          if(FALSE) // if calc DField:
            Dfield_.SetElemResult(pdeElem-1,DField);

          ptGrid_->CalcSurfNormal(normal, *surfElems[iel]);
          normal *= normSign;
          elemNormalD = normal * DField;

          // Integrate over DField * normal
          chargeOp->CalcElemCharge(charge, surfElems[iel], 
                                   lCoordSurf, elemNormalD);

          pdeElemNum = eqnData_->Mesh2PDEElem(ptVolElem->elemNum);
            
          Vector<Double> chargeVec(1);
          chargeVec[0] = charge;
          charges_.SetElemResult(pdeElemNum-1, chargeVec);
          chargeSD[isd] += charge;        

        }
       
      }

    Info->PrintF(couplingName_, "Computed surface charge: ");
    Info->PrintVec(chargeSD);
    
    delete FieldOp1;
    delete FieldOp2;
    
  } // end CalcCharges


  void PiezoCoupling::calcMaterialMatrices(Matrix<Double> & stiffnessMat,
                                           Matrix<Double> & piezoCouplingMat,
                                           Matrix<Double> & permittivityMat,
                                           Matrix<Double> * matDat){
    
    ENTER_FCN("PiezoCoupling::caclMaterialMatrices");
      
    UInt  stressDim, elecDim;    
    orientation2D actOrientation;

    // create material matrices

    if (subType_ == "planeStrain") 
      {
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
      Info->Error("StressOp: Unknown subtype in mech PDE! ",__FILE__,__LINE__);  

    stiffnessMat.Resize(stressDim,stressDim);
    piezoCouplingMat.Resize(elecDim,stressDim);
    permittivityMat.Resize(elecDim,elecDim);

    if (subType_ == "3d") {
      matDat->GetSubMatrix(stiffnessMat,0,0);
      matDat->GetSubMatrix(piezoCouplingMat,stressDim,0);
      matDat->GetSubMatrix(permittivityMat,stressDim,stressDim);
    } 

    else if (subType_ == "axi") {
      //stiffnessMat


      UInt rowPtrXY[]={1,2,6,3,7,8};
      UInt rowPtrYZ[]={2,3,4,1,8,9};
      UInt rowPtrXZ[]={1,3,5,2,7,9};
      UInt * rowPtr;   //alte Version

      switch(actOrientation) {

      case xy:
        {
          rowPtr = rowPtrXY;
          break;
        }
      case yz:
        {
          rowPtr = rowPtrYZ;
          break;
        }
      case xz:
        {
          rowPtr = rowPtrXZ;
          break;
        }
      default:    //if no orientation was specified
        {
          rowPtr = rowPtrYZ;
          break;
        }
      }
      

      Matrix<Double> tempMat;
      tempMat.Resize(matDat->GetSizeRow(),matDat->GetSizeCol());
      
      // Copy entries from material matrix object into temporary matrix
      for ( UInt i = 0; i < 6 ; i++ ) {
        for( UInt j = 0; j < 6; j++ ) {
          tempMat[i][j] = (*matDat)[rowPtr[i]-1][rowPtr[j]-1];
        }
      }
      tempMat.GetSubMatrix(stiffnessMat,0,0);
      tempMat.GetSubMatrix(piezoCouplingMat,stressDim,0);
      tempMat.GetSubMatrix(permittivityMat,stressDim,stressDim);
    }

    else if (subType_ == "planeStrain"){
      std::cout<<"PiezoCoupling::calcMaterialMatrices - plane strain ... still in progress"<<std::endl;
    
    }
  
}// end calcMaterialMatrices


  void PiezoCoupling::CalcComplexValuedCharges(){
    ENTER_FCN("PiezoCoupling::CalcComplexValuedCharges");


    BaseNodeStoreSol * solPDE1 = pde1_->getPDESolution();
    BaseNodeStoreSol * solPDE2 = pde2_->getPDESolution();
    eqnData_ = pde2_-> getPDE_eqnData(); 

    NodeStoreSol<Complex> * solhelp1 = dynamic_cast<NodeStoreSol<Complex>*>(solPDE1);    
    NodeStoreSol<Complex> * solhelp2 = dynamic_cast<NodeStoreSol<Complex>*>(solPDE2);
    
    // Determines gradient of electric potential, i.e. E=\grad \phi
    GradientFieldOp<Complex> * FieldOp2
      = new GradientFieldOp<Complex>(ptGrid_, pde2_,pde2_->getPDE_eqnData(),
                                     *solhelp2, ELEC_POTENTIAL, 
                                     isaxi_);

    // Determines linear Strain S=Bu, i.e.
    //  partial derivates of mechanical displacement

    LinStrainOp<Complex> * FieldOp1 
      = new LinStrainOp<Complex>(ptGrid_, pde1_,pde1_->getPDE_eqnData(),
                                 *solhelp1,MECH_DISPLACEMENT, 
                                 isaxi_);

    // ------ Calculation of the electric field and linear Strain ------

    Complex elemNormalD = 0.0;
    Complex  charge = 0.0;
    UInt pdeElemNum = 0;
    Elem * ptVolElem;
    BaseFE * ptSurfElemFE, * ptVolElemFE;
    SurfElem * ptSurfElem;
   
    StdVector<Elem*> elemssd;
    StdVector<SurfElem*> surfElems;
    Vector<Double>  normal,lCoordSurf, lCoordVol, lCoord, Coord;
    Vector<Complex> TempE, TempBu;
    Integer regionIndex = 0;
    UInt pdeElem;
    Double normSign = 0.0;

    //charge operator  
    ElecChargeOp * chargeOp;
    chargeOp = new ElecChargeOp(ptGrid_, pde2_, pde2_->getPDE_eqnData(), isaxi_);

    Vector<Complex> chargeSD(calcCharge_.GetSize());
    chargeSD.Init(0);

    if (subType_ == "axi") 
      isaxi_=TRUE;

    //! Matrices containing real parts of material
    Matrix<Double> stiffnessMatC;
    Matrix<Double> piezoCouplingMatC; 
    Matrix<Double> permittivityMatC;

    //! Matrices containing imaginary parts of material
    Matrix<Double> stiffnessMatD;
    Matrix<Double> piezoCouplingMatD; 
    Matrix<Double> permittivityMatD;

    //! complex material matrices
    Matrix<Complex> stiffnessMat;
    Matrix<Complex> piezoCouplingMat; 
    Matrix<Complex> permittivityMat;
    
    Matrix<Double> *matDat;

    if( params->HasValue( "type", "imagMaterialParameter", "materialDataType" ) ) {
      // calc real - parts
      matDat=pde1_->getPDEMaterialData()->GetMatrix();
      calcMaterialMatrices(stiffnessMatD, piezoCouplingMatD, permittivityMatD, matDat);
      // calc imag - parts
      matDat=pde1_->getPDEMaterialData()->GetMatrixC();
      calcMaterialMatrices(stiffnessMatC, piezoCouplingMatC, permittivityMatC, matDat);
      // resize Matrices
      stiffnessMat.Resize(stiffnessMatD.GetSizeRow(),stiffnessMatD.GetSizeRow());
      piezoCouplingMat.Resize(piezoCouplingMatD.GetSizeRow(), piezoCouplingMatD.GetSizeCol());
      permittivityMat.Resize(permittivityMatD.GetSizeRow(),permittivityMatD.GetSizeCol());

      // put them together
      for (UInt it=0;it<stiffnessMatC.GetSizeRow();it++)
        for (UInt jt=0;jt<stiffnessMatC.GetSizeRow();jt++){
          stiffnessMat[it][jt]=Complex(stiffnessMatD[it][jt], stiffnessMatC[it][jt]);
          
          if (it<piezoCouplingMatC.GetSizeRow()){
            piezoCouplingMat[it][jt] = Complex(piezoCouplingMatD[it][jt], piezoCouplingMatC[it][jt]);
            if (jt<permittivityMatC.GetSizeRow())
              permittivityMat[it][jt]=Complex(permittivityMatD[it][jt],permittivityMatC[it][jt]);
          }
        } // end for jt ...
    } // end if params->HasValue
    
    else { 
      matDat=pde1_->getPDEMaterialData()->GetMatrix();
      calcMaterialMatrices(stiffnessMatD, piezoCouplingMatD, permittivityMatD, matDat);
      
      stiffnessMat=stiffnessMatD;
      piezoCouplingMat = piezoCouplingMatD;
      permittivityMat = permittivityMatD;
    }
    
       // loop over all subdomains
    for (UInt isd=0; isd<calcCharge_.GetSize(); isd++)
      {       
        // get surface and according volume elements
        ptGrid_->GetSurfElems( surfElems, calcCharge_[isd] );

        complexValuedCharge_.Resize(surfElems.GetSize());

        // loop over all surface elements
        for (UInt iel=0; iel<surfElems.GetSize(); iel++){

          // Determine, which volume element is the right neighbour for the 
          // calculation
          if ( chargeNeighborRegion_.
               Find(surfElems[iel]->ptVolElem1->regionId) != -1 ) {
            ptVolElem = surfElems[iel]->ptVolElem1;
            normSign = -1.0;
          } else {
            ptVolElem = surfElems[iel]->ptVolElem2;
            normSign = 1.0;
          }
          
          normSign *= (double) surfElems[iel]->normalSign;

          ptSurfElemFE = surfElems[iel]->ptElem;
          ptVolElemFE = ptVolElem->ptElem;
          const StdVector<UInt> & surfConnect = surfElems[iel]->connect;
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
        
          // calculate electric field
          FieldOp2->CalcElemGradField(TempE, ptVolElem, lCoordVol,1);

          StdVector<UInt> connecth = ptVolElem->connect;
          //get coordinates of element
          Matrix<Double> Coord;
          pde1_->GetElemCoords(connecth, Coord);
        
          FieldOp1->CalcElemLinearStrain(TempBu, ptVolElem, Coord);

          Vector<Complex> DField;
          Vector<Complex> piezoCouplTimesStrain;

          if (subType_ == "3d"){
            DField.Resize(3);
            piezoCouplTimesStrain.Resize(3);
          }
          else
            {
              DField.Resize(2);
              piezoCouplTimesStrain.Resize(2);
            }

          piezoCouplingMat.Mult(TempBu,piezoCouplTimesStrain);
          permittivityMat.Mult(TempE,DField);
          DField+=piezoCouplTimesStrain;

//           std::cout<<"direct - piezo: DField" <<std::endl;
//           std::cout<<DField<<std::endl;

          //pdeElem = pde2_-> getPDE_eqnData()->Mesh2PDEElem(elemssd[iel]->elemNum);
            
          if(FALSE) // if calc DField:
            Dfield_.SetElemResult(pdeElem-1,DField);
        
          ptGrid_->CalcSurfNormal(normal, *surfElems[iel]);
          normal *= normSign;

          for(UInt i=0;i<DField.GetSize();i++)
            elemNormalD+=normal[i]*DField[i];                
        
          // Integrate over DField * normal
          chargeOp->CalcElemCharge(charge, surfElems[iel], 
                                   lCoordSurf, elemNormalD);

          pdeElemNum = eqnData_->Mesh2PDEElem(ptVolElem->elemNum);
            
          Vector<Complex> chargeVec(1);
          chargeVec[0] = charge;
          complexValuedCharge_[iel]=charge;
          chargesComplex_.SetElemResult(pdeElemNum-1, chargeVec);        
          chargeSD[isd] += charge;        
        
        }

      }
    // ! Writes result to StdPDE for later retrieval in SinglePDEs
    pde1_->setPDE_complexValuedCharge(complexValuedCharge_);

     Info->PrintF(couplingName_, " Computed surface charge:");
    Info-> PrintVec(chargeSD);

    delete FieldOp1;
    delete FieldOp2;
  }

  

  // *********************
  //   DefineIntegrators
  // *********************
  void PiezoCoupling::DefineIntegrators() {

    ENTER_FCN( "PiezoCoupling::DefineIntegrators" );
    
    piezoMaterialType matType = REALMATERIALPARAMETER;

    // iterate over all subdomains
    for ( UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++ ) {
      
      // add stiffness
      BaseForm *bilinearStiff = GetStiffIntegrator( &materialData_[actSD] );

      IntegratorDescriptor *actIntDescrStiff =
        new IntegratorDescriptor( bilinearStiff, STIFFNESS );

      bilinearStiff->SetPiezoMaterialType( matType );

      actIntDescrStiff->SetPiezoMaterialType( matType );
      actIntDescrStiff->SetPDEIds( pde1_, pde2_ );
      assemble_->AddIntegrator( actIntDescrStiff, subdoms_[actSD] );


      // Check for damping:

      //! list of damping types for all regions
      std::map<RegionIdType,DampingType> dampingList = pde1_->getPDE_dampingList();

      if ( dampingList[subdoms_[actSD]] == RAYLEIGH ) {

        Boolean isdamping = TRUE;
        Boolean reducedIntegration = FALSE; //is currently not supported
        BaseForm * dampStiff = 
           GetStiffIntegrator(&materialData_[actSD], reducedIntegration,isdamping );
        dampStiff->SetRaylDamping();
        
        IntegratorDescriptor *actIntDescrDamp =
          new IntegratorDescriptor(dampStiff, DAMPING);
        actIntDescrDamp->SetPDEIds(pde1_, pde1_);
        
        dampStiff->SetPiezoMaterialType(matType);
        actIntDescrDamp->SetPiezoMaterialType(matType);
        assemble_->AddIntegrator(actIntDescrDamp, subdoms_[actSD]);
      }

        // check for complex valued material parameter
      if( params->HasValue( "type", "imagMaterialParameter", "materialDataType" ) ) {
        matType = IMAGMATERIALPARAMETER; 

        BaseForm * bilinearStiffC = GetStiffIntegrator(&materialData_[actSD]);
        IntegratorDescriptor *actComplexIntDescrStiff = 
          new IntegratorDescriptor(bilinearStiffC, STIFFNESS);
        actComplexIntDescrStiff->SetPDEIds(pde1_, pde2_);
        actComplexIntDescrStiff->SetPiezoMaterialType(matType);
        bilinearStiffC->SetPiezoMaterialType(matType);
        assemble_->AddIntegrator(actComplexIntDescrStiff, subdoms_[actSD]);
      }
      
    }
  }


  // **********************
  //   GetStiffIntegrator
  // **********************
  BaseForm* PiezoCoupling::GetStiffIntegrator( MaterialData *actSDMat,
                                               Boolean reducedInt,
                                               Boolean isdamping ) {

    ENTER_FCN( "PiezoCoupling::GetStiffIntegrator" );

    // Get problem geometry and mechanic subtype
    std::string probGeo, subType;
    params->Get( "subType", subType, "mechanic" );
    params->Get( "type", probGeo, "geometry" );

    BaseForm *bilinearStiff;

    if (subType == "planeStrain") {
      bilinearStiff = new linPiezoCoupling2DPlaneStrain();
    }
    else if (subType == "axi") { 
      bilinearStiff = new linPiezoCoupling2DAxi();
    }
    else if (subType == "3d") {
      bilinearStiff = new linPiezoCoupling3D(); 
    }
    else {
      (*error) << "PiezoCoupling::GetStiffIntegrator:\n "
               << "Don't kn  ow how to handle subType = "
               << subType;
      Error( __FILE__, __LINE__ );
    }

    // Set pointer to material type
    bilinearStiff->SetMaterial( actSDMat );

    return bilinearStiff;
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

    std::string analysis;
    params->Get( "type", analysis, "analysis" );
    String2Enum(analysis,analysistype_);


    // *****************************
    // Determine nodal results
    // ***************************** 
    StdVector<std::string> nodeValues;
    keyVec  = couplingName_, "storeResults", "nodeResults", "region";
    attrVec = "", "", "type";  
  
    
    //  if ( calcDField_.GetSize() > 0 ) {
    //         hasOutput_ = TRUE;
    //         Info->PrintF( couplingName_, "Computing electric field for regions:\n" );

    //         for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
    //           Info->PrintF(couplingName_, "%s\n", regionNames[k].c_str() );
    //         }
    //         Info->PrintF( "", "\n" );

    //         std::cout<<"analysistype_"<<std::endl;
    //         std::cout<<analysistype_<<std::endl;
    //         if ( analysistype_ == HARMONIC || analysistype_==MULTIHARMONIC) {
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

    // --- Mechanic Stress ---
    //   Enum2String(MECH_STRESS, quantity);
    //       valVec  = "", "", quantity;
    //       params->GetList( keyVec, attrVec, valVec, regionNames );
    //       ptGrid_->RegionNameToId ( calcStress_, regionNames );
  
    //       // If the symbolic name is "all" compute stresses for all regions
    //       if ( calcStress_.GetSize() == 1 && calcStress_[0] == ALL_REGIONS ) {
    //         calcStress_ = subdoms_;
    //       }

    // Log to info file
    //   if ( calcStress_.GetSize() > 0 ) {
    //         hasOutput_ = TRUE;
    //         Info->PrintF(couplingName_,
    //                      "Computing mechanical stress for regions:\n" );
    //         for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
    //           Info->PrintF( couplingName_, "%s\n", regionNames[k].c_str() );
    //         }
    //         Info->PrintF( "", "\n" );

    //         if( analysistype_ == HARMONIC || analysistype_==MULTIHARMONIC) {

    //           // Resize solution arrays
    //           stressComplex_.SetNumSolutions(1);
    //           stressComplex_.SetSolutionType(MECH_STRESS);
    //           stressComplex_.SetNumElems(pde2_->getPDE_numElems());

    //           // We always store for six components (unverg-file-format as capa does
    //           stressComplex_.SetNumDofs(6);
    //           stressComplex_.SetPtrEQNData(eqnData_, ptGrid_);
    //           stressComplex_.Init(0);
    //         }

    //         else {

    //           // Resize solution arrays
    //           stress_.SetNumSolutions(1);
    //           stress_.SetSolutionType(MECH_STRESS);
    //           stress_.SetNumElems(pde2_->getPDE_numElems());

    //           // We always store for six components (unverg-file-format as capa does
    //           stress_.SetNumDofs(6);
    //           stress_.SetPtrEQNData(eqnData_, ptGrid_);
    //           stress_.Init(0);

    //         }
    //       }

    // --- Electric Charges ---
    // check for charge computation
    params->GetList( "region", regionNames, couplingName_, "charge" );
    ptGrid_->RegionNameToId( chargeNeighborRegion_, regionNames );

    params->GetList( "element", regionNames, couplingName_, "charge" );
    ptGrid_->RegionNameToId( calcCharge_, regionNames );

    
    if ( calcCharge_.GetSize() > 0 ) {
      
      hasOutput_ = TRUE;
      Info->PrintF( couplingName_, "Computing electric charges for regions:\n" );
      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( couplingName_, "%s\n", regionNames[k].c_str() );
      }
      Info->PrintF( "", "\n" );
      
      // Resize solution arrays
      if( analysistype_ == HARMONIC || analysistype_==MULTIHARMONIC) {
        chargesComplex_.SetNumSolutions(1);
        chargesComplex_.SetSolutionType(ELEC_CHARGE);
        chargesComplex_.SetNumElems(pde1_->getPDE_numElems());
        chargesComplex_.SetNumDofs(1);
        chargesComplex_.SetPtrEQNData(eqnData_, ptGrid_);
        chargesComplex_.Init();
      }
      else {
        charges_.SetNumSolutions(1);
        charges_.SetSolutionType(ELEC_CHARGE);
        charges_.SetNumElems(pde2_->getPDE_numElems());
        charges_.SetNumDofs(1);
        charges_.SetPtrEQNData(eqnData_, ptGrid_);
        charges_.Init();
      }
    } 

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
