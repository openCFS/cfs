#include <stdio.h>

#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "Forms/forms_header.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "Driver/assemble.hh"
#include "Driver/piezoParamIdent.hh"
#include "newmark.hh"
#include "Elements/basefe.hh"
#include "blocknodeEQN.hh"
#include "superblockEQN.hh"
#include "scalarblockEQN.hh"
#include "Forms/forms_header.hh"

#include "piezoPDE.hh" 

namespace CoupledField {

  PiezoPDE::PiezoPDE( Grid *aptgrid, BCs *aptbcs,  TimeFunc *aptTimeFunc,
		      FileType *aptFileType, WriteResults *aptOut ) :
    BasePDE( aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc ) {

    ENTER_FCN( "PiezoPDE::PiezoPDE" );

    // =====================================================================
    // set solution information
    // =====================================================================

    // Set name of the PDE and its material class

    pdename_ = "piezo";
    pdematerialclass_ = "piezo";

#ifndef XMLPARAMS
    // Determine dimension of problem
    conf->getstr( "subtype", subType_, pdename_ ); 
    if (subType_ == "3d") {
      dofspernode_ = 4;
      Info->PrintF("", "=== 3D PROBLEM\n");
    }
    else if (subType_ == "axi") {
      isaxi_ = TRUE;
      dofspernode_ = 3;
      Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
    }
    else {
      // default is planeStrain 
      dofspernode_ = 3;
      Info->PrintF("", "=== PLANE STRAIN PROBLEM\n");
    }

//     conf->getsubdompde(subdoms_,pdename_);
#else

    // Get problem geometry and PDE subtype
    params->Get( "subtype", subType_, pdename_ );
    std::string probGeo;
    params->Get( "type", probGeo, "geometry" );

    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry
    if ( subType_ == "3d" && probGeo == "3d" ) {
      dofspernode_ = 4;
      Info->PrintF("", "=== 3D PROBLEM\n");
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = TRUE;
      dofspernode_ = 3;
      Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
    }
    else if ( subType_ == "planeStrain" && probGeo == "plane" ) {
	dofspernode_ = 3;
	Info->PrintF("", "=== PLANE STRAIN PROBLEM\n");
    }
    else {
      std::string errmsg = "Subtype '" + subType_;
      errmsg += "' of PDE '" + pdename_ + "' does not fit to problem ";
      errmsg += "geometry '" + probGeo + "'\n";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

#endif

    // =====================================================================
    // set solution information
    // =====================================================================
    solTypes_ = MECH_DISPLACEMENT, ELEC_POTENTIAL;
    solDofs_ = dofspernode_-1, 1;
    
#ifndef XMLPARAMS
   effectiveMass_ = FALSE;
   if (conf->get_option("effMass",  pdename_ ))
     effectiveMass_ = TRUE;

   //check for damping model
   std::string dampstr;
   conf->ifget("damping",dampstr,pdename_);
   if (dampstr == "rayleigh")
     {
       //       Error("Currenrly Rayleigh damping not woirking for PiezoPDE",__FILE__,__LINE__);       
       dampingType_ = RAYLEIGH;
     }
   else
     dampingType_ = NONE;
   
   if (dampingType_)
     {
       needsDampingMatrix_ = TRUE;
     }
#else

   // Use effective mass approach?
   effectiveMass_ = params->IsSet( "effMass" );

   // Do we use damping?
   if( params->HasValue( "type", "rayleigh", pdename_, "damping" ) ) {
     dampingType_ = RAYLEIGH;
     Info->PrintF( pdename_, " Using RAYLEIGH damping\n" );
   }
   else {
     dampingType_ = NONE;
     Info->PrintF( pdename_, " Using no damping\n" );
   }
   if( dampingType_ != NONE ) {
     needsDampingMatrix_ = TRUE;
   }

#endif

#ifndef XMLPARAMS
#else

    //check for pressure loads
    params->GetList( "name"    , pressSurf_ , pdename_, "pressure" );
    params->GetList( "value"   , pressVals_ , pdename_, "pressure" );
    params->GetList( "dynamics", pressFnc_  , pdename_, "pressure" );

    // Check consistency
    if ( pressSurf_.GetSize() != pressVals_.GetSize() )
      {
	std::string errmsg = "PressureLoads: ";
	errmsg += "#name = " + Info->GenStr(pressSurf_.GetSize());
	errmsg += ", #value = " + Info->GenStr(pressVals_.GetSize());
	errmsg += ", #dynamics = " + pressFnc_.GetSize() + '\n';
	Info->Error( errmsg, __FILE__, __LINE__ );
      }

    // We need not have as many function/filenames as pressureloads!
    for ( Integer k = pressFnc_.GetSize(); k < pressSurf_.GetSize(); k++ )
      {
	pressFnc_.Push_back( "none" );
      }
#endif

  }


  // *********************
  //   DefineIntegrators
  // *********************
  void PiezoPDE::DefineIntegrators( const Integer level ) {
    ENTER_FCN( "PiezoPDE::DefineIntegerators" );

    Integer posOfElectricPot;
    if ( subType_ == "3d" )
      posOfElectricPot = 4;
    else
      posOfElectricPot = 3;

    Boolean nonLin = FALSE;
    for ( int actSD = 0; actSD < subdoms_.GetSize(); actSD++ ) {

      // ==============  add stiffness ========================================

      MaterialData actSDMat(materialData_[actSD]);

      // ==============  add "standard" stiffness =============================

      BaseForm * bilinearStiff = GetStiffIntegrator(actSDMat);
      IntegratorDescriptor *actIntDescrStiff =
	new IntegratorDescriptor(bilinearStiff, STIFFNESS);
      assemble_->AddIntegrator(actIntDescrStiff, subdoms_[actSD]);
      
      //check for damping
      if ( dampingType_ == RAYLEIGH ) {
	Boolean isdamping = TRUE;
	Boolean reducedIntegration = FALSE; //is currently not supported
	BaseForm * dampStiff = GetStiffIntegrator( actSDMat,reducedIntegration,
						   isdamping );
	IntegratorDescriptor *actIntDescrDamp =
	  new IntegratorDescriptor(dampStiff, DAMPING);
	assemble_->AddIntegrator(actIntDescrDamp, subdoms_[actSD]);
      }

      // ==============  add mass =============================================
      Double density = actSDMat.GetDensity();    
      BaseForm * bilinearMass  = new MassInt(density, dofspernode_,
					     posOfElectricPot, isaxi_);

      IntegratorDescriptor * actIntDescrMass =
	new IntegratorDescriptor(bilinearMass, MASS);

      //check for damping (mass part)
      if (dampingType_ == RAYLEIGH)    
	actIntDescrMass->SetSecondaryMat(DAMPING, actSDMat.GetDampingAlfa(),
					 analysistype_);
      assemble_->AddIntegrator(actIntDescrMass, subdoms_[actSD]);
    }

    //surface integrators
    //RHS-part
    Integer nonlin = 0;
    for (Integer actSF = 0; actSF < pressSurf_.GetSize(); actSF++) {
      BaseForm * rhsSrcSurf = new PressureLinForm(pressVals_[actSF], isaxi_);
      rhsSrcSurf->SetDofZero(posOfElectricPot);
      assemble_->AddRhsSrcSurfIntegrator(rhsSrcSurf, pressSurf_[actSF], 
					 pressFnc_[actSF],nonlin);
    }

  }


  // **********************
  //   GetStiffIntegrator
  // **********************
  BaseForm *PiezoPDE::GetStiffIntegrator( MaterialData& actSDMat,
					  Boolean reducedInt,
					  Boolean isdamping ) {

    ENTER_FCN( "PiezoPDE::GetStiffIntegrator" );
  
    BaseForm * bilinearStiff;
    if (subType_ == "planeStrain")
      bilinearStiff = new piezoPlainStrainInt(actSDMat,isdamping);

    else if (subType_ == "axi")      
      bilinearStiff = new piezoAxiInt(actSDMat,isdamping);
    
    else if (subType_ == "3d")
      bilinearStiff = new linPiezo3DInt(actSDMat,isdamping);

    else
      Error("Unknown subtype in mech PDE! ",__FILE__,__LINE__);

    return bilinearStiff;
  }


// ======================================================
// TRANSIENT SOLVING SECTION
// ======================================================


  void PiezoPDE::InitTimeStepping() {
    ENTER_FCN( "PiezoPDE::InitTimeStepping" );

    if (effectiveMass_)
      TS_alg_ = new NewmarkEffMass(pdename_, algsys_, 
				   eqnData_, needsDampingMatrix_);
    else
      TS_alg_ = new Newmark(pdename_, algsys_, eqnData_, needsDampingMatrix_);
  }


  void PiezoPDE::WriteResultsInFile(Integer stepOffset,
				    Double timeOffset)
  {
    ENTER_FCN( "PiezoPDE::WriteResultsInFile" );

    NodeStoreSol<Double> * solTransient;
    NodeStoreSol<Complex> * solHarmonic;

    Double actTime = lasttimecalc_ + timeOffset;
    Integer actStep = laststepcalc_ + stepOffset;
 
    if (analysistype_ == STATIC ||
	analysistype_ == TRANSIENT) {
	solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);      

	if (saveSol_)
	  outFile_->WriteNodeSolutionTransient(*solTransient, actStep, actTime);

	if (saveSolHist_)
	  outFile_->WriteNodeSolutionTransient(*solTransient, actStep, actTime);
    
	
	// Write derivatives
	if (analysistype_== TRANSIENT) {
	  if (saveDeriv1_ == TRUE)
	    {
	      solDeriv1_.SetAlgSysVector(getS1());
	      outFile_->WriteNodeSolutionTransient(solDeriv1_, actStep, actTime);
	      
	      if (saveDeriv1Hist_ == TRUE)
		outFile_->WriteNodeHistoryTransient(solDeriv1_, actStep, actTime);
	    }
	  
	  if (saveDeriv2_ == TRUE)
	    {
	      solDeriv2_.SetAlgSysVector(getS2());
	      outFile_->WriteNodeSolutionTransient(solDeriv2_, actStep, actTime);
	    }
	  if (saveDeriv2Hist_ == TRUE)
	    outFile_->WriteNodeHistoryTransient(solDeriv2_, actStep, actTime);
	}
	
	//element results
	if (calcEfield_.GetSize() !=0 ) 
	  outFile_->WriteElemSolutionTransient(Efield_, actStep, actTime);
	
	
	if (calcStress_.GetSize() !=0 ) 
	  outFile_->WriteElemSolutionTransient(stress_, actStep, actTime);
	
	
	//element results
	if (calcCharge_.GetSize() !=0 ) {
	  outFile_->WriteElemSolutionTransient(charges_, actStep, actTime);
	} 
	
    } else if (analysistype_ == HARMONIC) {
      
      if (saveSol_)
	{
	  solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);
	  outFile_->WriteNodeSolutionHarmonic(*solHarmonic,  actFreqStep_, 
					      actFrequency_, complexFormat_);
	}
      //element results
         if (calcCharge_.GetSize() !=0 ) {
      	outFile_->WriteElemSolutionHarmonic(chargesComplex_, actFreqStep_,  actFrequency_, complexFormat_);
       }
    }
  }



// ***********************************************************************
//   Obtain information on desired output quantities from parameter file
// ***********************************************************************
#ifdef XMLPARAMS
void PiezoPDE::ReadStoreResults() {

  ENTER_FCN( "PiezoPDE::ReadStoreResults" );


  // Construct vectors for restricted parameter search
  StdVector<std::string> keyVec;
  StdVector<std::string> attrVec;
  StdVector<std::string> valVec;
  std::string quantity;
  
  // *****************************
  // Determine nodal results
  // ***************************** 
  StdVector<std::string> nodeValues;
  keyVec  = pdename_, "storeResults", "nodeResults", "region";
  attrVec = "", "", "type";  

  // --- Mechanic Velocity ---
  Enum2String(MECH_VELOCITY, quantity);
  valVec = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, nodeValues);
  if (nodeValues.GetSize() > 0) {
    saveDeriv1_ = TRUE;
    
    // intialize corresponding storesolution object
    solDeriv1_.SetNumSolutions(1);
    solDeriv1_.SetNumNodes(numPDENodes_);
    solDeriv1_.SetSolutionType(MECH_VELOCITY);
    solDeriv1_.SetNumDofs(dofspernode_);
    solDeriv1_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
    solDeriv1_.Init();
  }
  
  // --- Mechanic Acceleration ---
  Enum2String(MECH_ACCELERATION, quantity);
  valVec = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, nodeValues);
  if (nodeValues.GetSize() > 0) {
    saveDeriv2_ = TRUE;
    
    // intialize corresponding storesolution object
    solDeriv2_.SetNumSolutions(1);
    solDeriv2_.SetNumNodes(numPDENodes_);
    solDeriv2_.SetSolutionType(MECH_ACCELERATION);
    solDeriv2_.SetNumDofs(dofspernode_);
    solDeriv2_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
    solDeriv2_.Init();
  }


  // Construct vectors for restricted search parameter

  //  keyVec  = "piezo", "storeResults", "elemResults", "region";
  attrVec = "", "", "type";

  // *****************************
  // Determine element results
  // *****************************
  StdVector<std::string> elemResults;
  keyVec  = pdename_, "storeResults", "elemResults", "region";


  // --- Mechanic Stress ---
  Enum2String(MECH_STRESS, quantity);
  valVec  = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, calcStress_ );


  // If the symbolic name is "all" compute electric field for all regions
  if ( calcStress_.GetSize() == 1 && calcStress_[0] == "all" ) {
    calcStress_ = subdoms_;
  }

  // Log to info file
  if ( calcStress_.GetSize() > 0 ) {
    Info->PrintF( pdename_,
		  " Computing mechanical stress for regions:");
    for ( Integer k = 0; k < calcStress_.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s", calcStress_[k].c_str() );
    }
  }


  // --- Electric Field Intensity ---
  Enum2String(ELEC_FIELD_INTENSITY, quantity);
  valVec  = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, calcEfield_ );

  // If the the symbolic name is "all" compute electric field for all regions
  if ( calcEfield_.GetSize() == 1 && calcEfield_[0] == "all" ) {
    calcEfield_ = subdoms_;
  }

  if ( calcEfield_.GetSize() > 0 ) {
    Info->PrintF( pdename_, " Computing electric field for regions:" );
    for ( Integer k = 0; k < calcEfield_.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s", calcEfield_[k].c_str() );
    }
    if (analysistype_ == HARMONIC){
    EfieldComplex_.SetNumSolutions(1);
    EfieldComplex_.SetSolutionType(ELEC_FIELD_INTENSITY);
    EfieldComplex_.SetNumElems(numElems_);
    EfieldComplex_.SetNumDofs(dim_);
    EfieldComplex_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
    EfieldComplex_.Init(); 
    }
    else{
    Efield_.SetNumSolutions(1);
    Efield_.SetSolutionType(ELEC_FIELD_INTENSITY);
    Efield_.SetNumElems(numElems_);
    Efield_.SetNumDofs(dim_);
    Efield_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
    Efield_.Init(); 
    }
  }

  // --- Mechanic Stress ---
  Enum2String(MECH_STRESS, quantity);
  valVec  = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, calcStress_ );
  
  // If the symbolic name is "all" compute stresses for all regions
  if ( calcStress_.GetSize() == 1 && calcStress_[0] == "all" ) {
    calcStress_ = subdoms_;
  }

  // Log to info file
  if ( calcStress_.GetSize() > 0 ) {
    Info->PrintF( pdename_,
		  " Computing mechanical stress for regions:");
    for ( Integer k = 0; k < calcStress_.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s", calcStress_[k].c_str() );
    }
    if(analysistype_ == HARMONIC){
      // Resize solution arrays
      stressComplex_.SetNumSolutions(1);
      stressComplex_.SetSolutionType(MECH_STRESS);
      stressComplex_.SetNumElems(numElems_);
      //we always store for six components (unverg-file-format as capa does
      stressComplex_.SetNumDofs(6);
      stressComplex_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
      stressComplex_.Init(0);
    }
    else{

      // Resize solution arrays
      stress_.SetNumSolutions(1);
      stress_.SetSolutionType(MECH_STRESS);
      stress_.SetNumElems(numElems_);
      //we always store for six components (unverg-file-format as capa does
      stress_.SetNumDofs(6);
      stress_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
      stress_.Init(0);

    }
  }
  
  // --- Electric Charges ---
  //check for charge computation
  params->GetList( "region", chargeNeighborRegion_, pdename_, "charge" );
  params->GetList( "element", calcCharge_, pdename_, "charge" );

  if (calcCharge_.GetSize() > 0)
    {
      Info->PrintF( pdename_,
		    " Computing electric charges for regions:");
      for ( Integer k = 0; k < calcCharge_.GetSize(); k++ ) {
	Info->PrintF( pdename_, " %s", calcCharge_[k].c_str() );
      } 
      // Resize solution arrays
      if(analysistype_==HARMONIC)
	{
	  chargesComplex_.SetNumSolutions(1);
	  chargesComplex_.SetSolutionType(ELEC_CHARGE);
	  chargesComplex_.SetNumElems(numElems_);
	  chargesComplex_.SetNumDofs(1);
	  chargesComplex_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
	  chargesComplex_.Init();
	}
      else
       {
	charges_.SetNumSolutions(1);
	charges_.SetSolutionType(ELEC_CHARGE);
	charges_.SetNumElems(numElems_);
	charges_.SetNumDofs(1);
	charges_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
	charges_.Init();
      }
    } 
  
  // *****************************
  // Determine nodal history
  // *****************************
  StdVector<std::string> saveNodeHist; 
  keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
  attrVec = "", "", "type";
  
  // --- mechDisplacement / elecPotential ---
  Enum2String(MECH_DISPLACEMENT, quantity);
  valVec  = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
  if (saveNodeHist.GetSize() > 0) {
    if (saveSol_ == FALSE) {
      std::string errMsg = pdename_;
      errMsg += ": History of ";
      errMsg += quantity + " can only be written, if it is activated ";
      errMsg += "in section 'nodalResults', too.";
      Error(errMsg.c_str(), __FILE__, __LINE__);
      }
    saveSolHist_ = TRUE;
    Info->PrintF( pdename_, " Saving mechDisplacement for Nodes:" );
    for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s", saveNodeHist[k].c_str() );
    }
  }
  
  // --- mechVelocity ---
  Enum2String(MECH_VELOCITY, quantity);
  valVec  = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
  if (saveNodeHist.GetSize() > 0) {
    if (saveDeriv1_ == FALSE) {
      std::string errMsg = pdename_;
      errMsg += ": History of ";
      errMsg += quantity + " can only be written, if it is activated ";
      errMsg += "in section 'nodalResults', too.";
      Error(errMsg.c_str(), __FILE__, __LINE__);
      }
    saveDeriv1Hist_ = TRUE;
    Info->PrintF( pdename_, " Saving mechVelocity for Nodes:" );
    for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s", saveNodeHist[k].c_str() );
    }
  }

  // --- mechAcceleration ---
  Enum2String(MECH_ACCELERATION, quantity);
  valVec  = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
  if (saveNodeHist.GetSize() > 0) {
    if (saveDeriv2_ == FALSE) {
      std::string errMsg = pdename_;
      errMsg += ": History of ";
      errMsg += quantity + " can only be written, if it is activated ";
      errMsg += "in section 'nodalResults', too.";
      Error(errMsg.c_str(), __FILE__, __LINE__);
    }
    saveDeriv1Hist_ = TRUE;
    Info->PrintF( pdename_, " Saving mechAcceleration for Nodes:" );
    for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s", saveNodeHist[k].c_str() );
    }
  } 

  // *****************************
  // Determine element history
  // *****************************
  StdVector<std::string> saveElemHist;
  keyVec  = pdename_, "storeResults", "elemHistory", "saveElems";
  attrVec = "", "", "";
  valVec = "", "", "";
  params->GetList(keyVec, attrVec, valVec, saveElemHist);

  if (saveElemHist.GetSize() < 0) {
    std::string errMsg = pdename_;
    errMsg += ": Saving history elements is not implemented yet!\n";
    errMsg += "Meanwhile you can use 'unvtool' to extract element data.";
    Error( errMsg.c_str(), __FILE__, __LINE__);
  }

}

#endif

// ************************************************************
//   PostProcess
// ************************************************************

void PiezoPDE::PostProcess(const Integer level) {

  ENTER_FCN( "PiezoPDE::PostProcess" );

  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);


  // calc electric field
  if (calcEfield_.GetSize() !=0 ) {
    if (analysistype_==HARMONIC)
      CalcComplexValuedEfield();
    else
      CalcEfield();
  }
  
  // calc stresses
  if (calcStress_.GetSize() !=0 ) {
    if (analysistype_ == HARMONIC)
      CalcComplexValuedStress();
    else
      CalcStress();
  }
    

  //calc charges
  if (calcCharge_.GetSize() !=0 ) {
   
    if (analysistype_ == HARMONIC)
      CalcComplexValuedCharges();
    else
      CalcCharges();
  }      
}


void PiezoPDE::CalcEfield(){
  ENTER_FCN( "PiezoPDE::CalcEfield" );
  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

   GradientFieldOp<Double> * FieldOp = new GradientFieldOp<Double>(ptgrid_, this, eqnData_,
						  *solhelp, ELEC_POTENTIAL, 
						  actlevel_, isaxi_);


  // ------ Calculation of the electric field ------

  Vector<Double> LCoord;
      
  StdVector<Elem*> elemssd;
  Integer counterElems=0;
  Vector<Double> TempE;
  Integer pdeElem;
  
  // loop over all subdomains
  for (Integer isd=0; isd<calcEfield_.GetSize(); isd++)
    {
      // get vector of Elem of subdomain with color: subdoms[isd]
      ptgrid_->GetElemSD(elemssd,calcEfield_[isd],actlevel_);
      
      // loop over elements of subdomain
      for (Integer iel=0; iel< elemssd.GetSize(); iel++,counterElems++)
	{
	  elemssd[iel]->ptElem->GetCoordMidPoint(LCoord);
	  FieldOp->CalcElemGradField( TempE, elemssd[iel], LCoord, 1); 
	  pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
	  Efield_.SetElemResult(pdeElem-1,TempE);
	}
    }
  delete FieldOp;
}


void PiezoPDE::CalcComplexValuedEfield(){
  ENTER_FCN( "PiezoPDE::CalcComplexValuedEfield" );
  NodeStoreSol<Complex> * solhelp = dynamic_cast<NodeStoreSol<Complex>*>(sol_);

  GradientFieldOp<Complex> * FieldOp = new GradientFieldOp<Complex>(ptgrid_, this, eqnData_,
						  *solhelp, ELEC_POTENTIAL, 
						  actlevel_, isaxi_);


  // ------ Calculation of the electric field ------

  Vector<Double> LCoord;
 
      
  StdVector<Elem*> elemssd;
  Integer counterElems=0;
  Vector<Complex> TempE;
  Integer pdeElem;
  
  // loop over all subdomains
  for (Integer isd=0; isd<calcEfield_.GetSize(); isd++)
    {
      // get vector of Elem of subdomain with color: subdoms[isd]
      ptgrid_->GetElemSD(elemssd,calcEfield_[isd],actlevel_);
      
      // loop over elements of subdomain
      for (Integer iel=0; iel< elemssd.GetSize(); iel++,counterElems++)
	{
	  elemssd[iel]->ptElem->GetCoordMidPoint(LCoord);
	  FieldOp->CalcElemGradField(TempE, elemssd[iel], LCoord, 1); 
	  pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
	  EfieldComplex_.SetElemResult(pdeElem-1,TempE);
	}
    }
  delete FieldOp;
}


void PiezoPDE::CalcStress(){
  ENTER_FCN( "PiezoPDE::CalcSress" );
   NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
  
  
  //get the correct bilinearform
  ShortInt stressElecDim, stressDim, elecDim;
  Vector<Double> intPoint;
  
#ifndef XMLPARAMS 
  if (subType_ == "plainStrain") 
#else
    if (subType_ == "planeStrain") 
#endif
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
  
  
  Vector<Double> elemElecStress, elemStress, sortedStress;
  elemElecStress.Resize(stressDim+elecDim);
  elemStress.Resize(stressDim);
  elemElecStress.Init(0);
  elemStress.Init(0);
  sortedStress.Resize(6);

  // loop over all subdomains
  for (Integer isd=0; isd<subdoms_.GetSize(); isd++) {
    
    MaterialData actSDMat(materialData_[isd]);
    linPiezoInt * stress;

    if (subType_ == "planeStrain")
      stress = new piezoPlainStrainInt(actSDMat);

    else if (subType_ == "axi")      
      stress = new piezoAxiInt(actSDMat);
    
    else if (subType_ == "3d")
      stress = new linPiezo3DInt(actSDMat);
    
    // get vector of Elements of subdomains
    StdVector<Elem*> elemssd;     
    ptgrid_->GetElemSD(elemssd,subdoms_[isd],actlevel_);
    
    // loop over elements of subdomain
    for (Integer iel=0; iel< elemssd.GetSize(); iel++) {
      Integer pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
      elemssd[iel]->ptElem->GetCoordMidPoint(intPoint);

      //set element pointer
      BaseFE * ptEl = elemssd[iel]->ptElem;
      stress->SetElemPtr(ptEl);
      
      //set element solution	
      Matrix<Double> elSol;
      StdVector<Integer> connecth = elemssd[iel]->connect;
      sol_->GetElemSolutionAsMatrix(elSol, connecth);
      stress->SetActElemSol(elSol);
      
      //get coordinates of element
      Matrix<Double> ptCoord;
      GetElemCoords(connecth, ptCoord, actlevel_);
      
      Vector<Double> actStress;	
      
      //set the integration point
      stress->SetIntPoint(intPoint);
      
      //calculates the stress
      stress->CalcStressVec(elemElecStress,1,ptCoord);
      elemStress = elemElecStress.Part(0,stressDim-1);
      sortStresses(elemStress,sortedStress);
      stress_.SetElemResult(pdeElem-1, sortedStress);
    }
  }
}
  
void PiezoPDE::CalcComplexValuedStress(){
  ENTER_FCN( "PiezoPDE::CalcComplexValuedStress" );
  NodeStoreSol<Complex> * solhelp = dynamic_cast<NodeStoreSol<Complex>*>(sol_);
  
  //get the correct bilinearform
  ShortInt stressElecDim, stressDim, elecDim;
  Vector<Double> intPoint;
  
#ifndef XMLPARAMS 
  if (subType_ == "plainStrain") 
#else
    if (subType_ == "planeStrain") 
#endif
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
  
  
  Vector<Complex> elemElecStress, elemStress, sortedStress;
  elemElecStress.Resize(stressDim+elecDim);
  elemStress.Resize(stressDim);
  elemElecStress.Init(0);
  elemStress.Init(0);
  sortedStress.Resize(6);

  // loop over all subdomains
  for (Integer isd=0; isd<subdoms_.GetSize(); isd++) {
    
    MaterialData actSDMat(materialData_[isd]);
    linPiezoInt * stress;

    if (subType_ == "planeStrain")
      stress = new piezoPlainStrainInt(actSDMat);

    else if (subType_ == "axi")      
      stress = new piezoAxiInt(actSDMat);
    
    else if (subType_ == "3d")
      stress = new linPiezo3DInt(actSDMat);
    
    // get vector of Elements of subdomains
    StdVector<Elem*> elemssd;     
    ptgrid_->GetElemSD(elemssd,subdoms_[isd],actlevel_);
    
    // loop over elements of subdomain
    for (Integer iel=0; iel< elemssd.GetSize(); iel++) {
      Integer pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
      elemssd[iel]->ptElem->GetCoordMidPoint(intPoint);
	    
      //set element pointer
      BaseFE * ptEl = elemssd[iel]->ptElem;
      stress->SetElemPtr(ptEl);
      
      //set element solution	
      Matrix<Complex> elSol;
      StdVector<Integer> connecth = elemssd[iel]->connect;
      sol_->GetElemSolutionAsMatrix(elSol, connecth);
      stress->SetActElemSol(elSol);
      
      //get coordinates of element
      Matrix<Double> ptCoord;
      GetElemCoords(connecth, ptCoord, actlevel_);
      
      Vector<Double> actStress;	
      
      //set the integration point
      stress->SetIntPoint(intPoint);
      
      //calculates the stress
  
      stress->CalcStressVec(elemElecStress,1,ptCoord);

      elemStress = elemElecStress.Part(0,stressDim-1);
      sortStresses(elemStress,sortedStress);
      stressComplex_.SetElemResult(pdeElem-1, sortedStress);

    }
  }
}
  
  
void PiezoPDE::CalcCharges(){
  ENTER_FCN( "PiezoPDE::CalcCharges" );
    
  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
  StdVector<Elem*> surfElems, volElems;
  Vector<Double> lCoordSurf, lCoordVol, elemDField;
  GradientFieldOp<Double> * dFieldOp;
  ElecChargeOp * chargeOp;
  BaseFE * ptSurfElem, * ptVolElem;
  Double permittivity = 0.0;
  Double elemNormalD = 0.0;
  Double charge = 0.0;
  Integer pdeElemNum = 0;

  ShortInt stressElecDim, stressDim, elecDim;
  Vector<Double> intPoint;
  
#ifndef XMLPARAMS 
  if (subType_ == "plainStrain") 
#else
    if (subType_ == "planeStrain") 
#endif
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
  
  
  Vector<Double> elemElecStress;
  elemElecStress.Resize(stressDim+elecDim);
  elemElecStress.Init(0);

    
  // Create vector with interpolation coordinate.
  // For simplicity we only evaluate the integral
  // in coordinate origin
  lCoordSurf.Resize(dim_-1);
  lCoordSurf.Init(0);

  //charge operator  
  chargeOp = new ElecChargeOp(ptgrid_, this, eqnData_, actlevel_, isaxi_);
			      
  // loop over all subdomains
  for (Integer iSD=0; iSD<calcCharge_.GetSize(); iSD++){
    
    // get surface and acoording volume elements
    if (dim_ == 3)
      surfElems = ptBCs_->getFacesBC(calcCharge_[iSD], actlevel_);
    else if (dim_ == 2)
      surfElems = ptBCs_->getEdgesBC(calcCharge_[iSD], actlevel_);
    
    // get neighbouring volume elements of
    // surface elements
    ptgrid_->GetVolNeighboursForSurf(surfElems,chargeNeighborRegion_,
				     volElems, actlevel_);
    
    //get correct stress-
    MaterialData actSDMat(materialData_[iSD]);
    linPiezoInt * stress;

    if (subType_ == "planeStrain")
      stress = new piezoPlainStrainInt(actSDMat);

    else if (subType_ == "axi")      
      stress = new piezoAxiInt(actSDMat);
    
    else if (subType_ == "3d")
      stress = new linPiezo3DInt(actSDMat);

    // loop over all surface elements
    for (Integer iElem=0; iElem<surfElems.GetSize(); iElem++)
      {
	
	ptSurfElem = surfElems[iElem]->ptElem;
	ptVolElem = volElems[iElem]->ptElem;
	const StdVector<Integer> & surfConnect = surfElems[iElem]->connect;
	const StdVector<Integer> & volConnect = volElems[iElem]->connect;

	// calculate volume integration coordinates from
	// surfe integration coordinates for evaluating the 
	// electric flux density on the surface of the volume
	// element
	ptVolElem->GetLocalIntPoints4Surface(surfConnect, volConnect,
					     lCoordSurf, lCoordVol);

	// Get the right material parameter for actual volume element
	//	for (Integer i=0; i<subdoms_.GetSize(); i++)
	//	  {
	//	    if (subdoms_[i] == volElems[iElem]->namesd)
	//	      permittivity  = materialData_[iSD].GetPermittivity(2,2);
	//	  }
	
	// Calc electric flux density
	//	dFieldOp->CalcElemGradField(elemDField, volElems[iElem], 
	//				    lCoordVol, permittivity);
	
	//elemNormalD = lCoordVol * elemDField;


	//set volume element
	stress->SetElemPtr(ptVolElem);
      
	//set element solution	
	Matrix<Double> elSol;
	sol_->GetElemSolutionAsMatrix(elSol, volConnect);
	stress->SetActElemSol(elSol);
    
	//set the integration point
	stress->SetIntPoint(lCoordVol);

	//get coordinates of element
	Matrix<Double> ptCoord;
	GetElemCoords(volConnect, ptCoord, actlevel_);

	//calculates the stress
	Vector<Double> actElecD;
	
	stress->CalcStressVec(elemElecStress,1,ptCoord);

	actElecD = elemElecStress.Part(stressDim,elemElecStress.GetSize()-1);
      
	//scalar product with normal vector
	elemNormalD = lCoordVol * actElecD;


	chargeOp->CalcElemCharge(charge, surfElems[iElem], 
				 lCoordSurf, elemNormalD);

	pdeElemNum = eqnData_->Mesh2PDEElem(volElems[iElem]->elemNum);
	
	// Create temporary vector, since SetElemResult only
	// can handle these
	Vector<Double> chargeVec(1);
	chargeVec[0] = charge;
	charges_.SetElemResult(pdeElemNum-1, chargeVec);
	
      }
  }
  //  Warning ("Charges are written to unv/gmv file, although capapost \
// can not draw them", __FILE__, __LINE__);
 }



void PiezoPDE::CalcComplexValuedCharges(){
  ENTER_FCN( "PiezoPDE::CalcComplexValuedCharges" );
  NodeStoreSol<Complex> * solhelp = dynamic_cast<NodeStoreSol<Complex>*>(sol_);
  StdVector<Elem*> surfElems, volElems;
  Vector<Double> lCoordSurf, lCoordVol;
  Vector<Complex> elemDField;
  GradientFieldOp<Complex> * dFieldOp;
  ElecChargeOp * chargeOp;
  BaseFE * ptSurfElem, * ptVolElem;
  Double permittivity = 0.0;
  Complex   elemNormalD = Complex(0.0,0.0);
  Complex charge = Complex(0.0,0.0);
  Integer pdeElemNum = 0;

  ShortInt stressElecDim, stressDim, elecDim;
  Vector<Double> intPoint;
  
#ifndef XMLPARAMS 
  if (subType_ == "plainStrain") 
#else
    if (subType_ == "planeStrain") 
#endif
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
  
  Vector<Complex> elemElecStress;
  elemElecStress.Resize(stressDim+elecDim);
  elemElecStress.Init(0);
    
  // Create vector with interpolation coordinate.
  // For simplicity we only evaluate the integral
  // in coordinate origin

  lCoordSurf.Resize(dim_-1);
  lCoordSurf.Init(0);

  //charge operator  
  chargeOp = new ElecChargeOp(ptgrid_, this, eqnData_, actlevel_, isaxi_);
			      
  // loop over all subdomains
  for (Integer iSD=0; iSD<calcCharge_.GetSize(); iSD++){
    
    // get surface and acoording volume elements
    if (dim_ == 3)
      surfElems = ptBCs_->getFacesBC(calcCharge_[iSD], actlevel_);
    else if (dim_ == 2)
      surfElems = ptBCs_->getEdgesBC(calcCharge_[iSD], actlevel_);
    
    // get neighbouring volume elements of
    // surface elements
    ptgrid_->GetVolNeighboursForSurf(surfElems,chargeNeighborRegion_,
				     volElems, actlevel_);
    
    //get correct stress-
    MaterialData actSDMat(materialData_[iSD]);
    linPiezoInt * stress;

    if (subType_ == "planeStrain")
      stress = new piezoPlainStrainInt(actSDMat);

    else if (subType_ == "axi")      
      stress = new piezoAxiInt(actSDMat);
    
    else if (subType_ == "3d")
      stress = new linPiezo3DInt(actSDMat);

    // loop over all surface elements

  complexValuedCharge_.Resize(surfElems.GetSize());

    for (Integer iElem=0; iElem<surfElems.GetSize(); iElem++)
      {
	
	ptSurfElem = surfElems[iElem]->ptElem;
	ptVolElem = volElems[iElem]->ptElem;
	const StdVector<Integer> & surfConnect = surfElems[iElem]->connect;
	const StdVector<Integer> & volConnect = volElems[iElem]->connect;

	// calculate volume integration coordinates from
	// surfe integration coordinates for evaluating the 
	// electric flux density on the surface of the volume
	// element
	ptVolElem->GetLocalIntPoints4Surface(surfConnect, volConnect,
					     lCoordSurf, lCoordVol);

	// Get the right material parameter for actual volume element
	//	for (Integer i=0; i<subdoms_.GetSize(); i++)
	//	  {
	//	    if (subdoms_[i] == volElems[iElem]->namesd)
	//	      permittivity  = materialData_[iSD].GetPermittivity(2,2);
	//	  }
	
	// Calc electric flux density
	//	dFieldOp->CalcElemGradField(elemDField, volElems[iElem], 
	//				    lCoordVol, permittivity);
	
	//elemNormalD = lCoordVol * elemDField;


	//set volume element
	stress->SetElemPtr(ptVolElem);
      
	//set element solution	
	Matrix<Complex> elSol;
	sol_->GetElemSolutionAsMatrix(elSol, volConnect);

	stress->SetActElemSol(elSol);

	//set the integration point
	stress->SetIntPoint(lCoordVol);

	//get coordinates of element
	Matrix<Double> ptCoord;
	GetElemCoords(volConnect, ptCoord, actlevel_);

	//calculates the stress
	Vector<Complex> actElecD;	
	stress->CalcStressVec(elemElecStress,1,ptCoord);
	actElecD = elemElecStress.Part(stressDim,elemElecStress.GetSize()-1);

	elemNormalD = Complex();
	
	// scalarProduct lCoordCal*actElecD
	for (Integer i=0; i<lCoordVol.GetSize(); i++)
	     elemNormalD += lCoordVol[i] * actElecD[i];


	//lCoordVol.Mult(actElecD, elemNormalD);
	chargeOp->CalcElemCharge(charge, surfElems[iElem], 
				 lCoordSurf, elemNormalD);

	pdeElemNum = eqnData_->Mesh2PDEElem(volElems[iElem]->elemNum);

	
	// Create temporar vector, since SetElemResult only
	// can handle these
	Vector<Complex> chargeVec(1);
	chargeVec[0] = charge;
	complexValuedCharge_[iElem]=charge;
       	chargesComplex_.SetElemResult(pdeElemNum-1, chargeVec);	

      }
  }
  //  Warning ("Charges are written to unv/gmv file, although capapost \
//can not draw them", __FILE__, __LINE__);

}


  
  

} // end of namespace
