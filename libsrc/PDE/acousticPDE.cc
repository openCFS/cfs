#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acousticPDE.hh" 
#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "Forms/forms_header.hh"
#include "Estimator/spaceerror.hh"
#include "newmark.hh"
#include "newmarkFracDamp.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/scalarnodeEQN.hh"
#include "Utils/mathfunctions.hh"
#include "Utils/nodestoresol.hh"
#include "Driver/solveStepAcoustic.hh"

namespace CoupledField {

  // =====================================================================
  // set solution information
  // =====================================================================
  AcousticPDE::AcousticPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc,
			   FileType *aptFileType, WriteResults *aptOut)
    :BasePDE(aptgrid,aptbcs,aptFileType,aptOut,aptTimeFunc) {

    ENTER_FCN( "AcousticPDE::AcousticPDE" );

    dofspernode_ = 1;
    solDofs_ = 1;
    pdename_          = "acoustic";
    pdematerialclass_ = "fluid";
    laststepcalc_ = 0;

    coarsealpha_ = 0.01; // solver parameter, see basePDE.cc

    // PDE formulation either in acoustic potential or pressure
    std::string s_;
    params->Get( "formulation", s_, "pdeList", pdename_ );
    String2Enum( s_, formulation_ );
    s_ = "Using * " + s_ + " as state variable in formulation of PDE\n";
    Info->PrintF( pdename_, s_.c_str() );

    // class NodeStoreSol will be initialized with acoustic potential
    solTypes_ = ACOU_POTENTIAL;

    // timestepping formulation
    params->Get( "timeSteppingFormulation", s_, "pdeList", pdename_ );
    if ( s_ == "effMassMatrix" ) {
      effectiveMass_ = TRUE;
      Info->PrintF( pdename_, "      * effective mass matrix timestepping\n");
    } 
    else {
      effectiveMass_ = FALSE;
      Info->PrintF( pdename_, "      * effective stiffness matrix timestepping\n");
    }

    // axisymmetric setup    
    isaxi_ = params->HasValue( "type", "axi", "geometry" );

    // *********************************************
    //   Check what type of damping should be used
    // *********************************************

    dampingType_ = NONE;
    params->GetList( "type", dampingList_, pdename_, "damping");

    if ( dampingList_.GetSize() == 0 )
      Info->PrintF( pdename_, "No Damping information detected!\n" );

    Boolean sorted_ = TRUE;
    for ( Integer k = 0; k < dampingList_.GetSize(); k++) {
      if (dampingList_[k] == "rayleigh") {
	dampingType_ = RAYLEIGH;
	Info->PrintF( pdename_, "      * RAYLEIGH damping for region: %d\n", k );
      }
      else if (dampingList_[k] == "thermoViscous") {
	dampingType_ = THERMOVISCOUS;
	Info->PrintF( pdename_, "      * THERMOVISCOUS damping for region: %d\n", k );
      }
      else if (dampingList_[k] == "fractional") {
	dampingType_ = FRACTIONAL;
	Info->PrintF( pdename_, "      * FRACTIONAL damping for region: %d\n", k );
      }
      else if (dampingList_[k] == "no") {
	dampingType_ = NONE;
	Info->PrintF( pdename_, "      * NO damping at all for region: %d\n", k );
      }
      if ( k > 1 )
	if ( dampingList_[k] != dampingList_[k-1] )
	  sorted_ = FALSE;
    }
    if ( sorted_ == FALSE )
      Error("Specify same type of damping for all regions!",__FILE__,__LINE__);

    // get additional information for fractional damping model
    if ( dampingType_ == FRACTIONAL ) {

      StdVector<std::string> fracAlgList_;
      params->GetList( "fracAlg", fracAlgList_, pdename_, "damping" );
      StdVector<Integer> fracMemoryList_;
      params->GetList( "fracMemory", fracMemoryList_, pdename_, "damping" );
      StdVector<std::string> interpolationList_;
      params->GetList( "interpolation", interpolationList_, pdename_, "damping" );

      if( fracAlgList_.GetSize() == 0 || fracMemoryList_.GetSize() == 0
	  || interpolationList_.GetSize() == 0 )
	Error("Specify attributes fracAlg, fracMemory and interpolation!",__FILE__,__LINE__);
      // up to now take values from first subdomain
      else {
	if ( fracAlgList_[0] == "gl" )
	  Info->PrintF( pdename_, "         with Gruenwald-Letnikov algorithm,\n");
	else if (fracAlgList_[0] == "blank")
	  Info->PrintF( pdename_, "         with Blanks algorithm,\n");

	fracMemory_ = fracMemoryList_[0];
	Info->PrintF( pdename_, "         memory size is: %d,\n", fracMemory_);

	if ( interpolationList_[0] == "lin1pt")
	  inType_ = LIN1PT;
	else
	  inType_ = NOTUSED;
	Info->PrintF( pdename_, "         %s interpolation of past values\n\n"
		      , interpolationList_[0].c_str() );
      }
      // modify dampingList, so that fracAlg is included
      for ( Integer k = 0; k < dampingList_.GetSize(); k++) {
	if ( fracAlgList_[k] == "gl" )
	  dampingList_[k] = "fractional_gl";
	else if (fracAlgList_[k] == "blank")
	  dampingList_[k] = "fractional_blank";
      }
    }

    // *************************************************************
    //   Check what type of nonlinear PDE formulation should be used
    // *************************************************************

    nonLin_ = FALSE; //declaration in basePDE.hh
    params->GetList( "nonLinear", nonLinPDEType_, pdename_, "region" );

    for ( Integer k = 0; k < nonLinPDEType_.GetSize(); k++ ) {
      if ( nonLinPDEType_[k] != "no" )
	nonLin_ = TRUE;

      if ( nonLinPDEType_[k] != "kuznetsov" && solTypes_ == ACOU_PRESSURE )
	Error("Acoustic pressure formulation not supported for Kuznetsov equation!"
	      ,__FILE__,__LINE__);

      if ( nonLinPDEType_[k] != "westervelt" && solTypes_ == ACOU_POTENTIAL )
	Error("Acoustic potential formulation not supported for Westervelt equation!"
	      ,__FILE__,__LINE__);
    }
    if( nonLin_ ) {
      // solution method
      params->Get("method", nonLinMethod_, pdename_, "nonLinear" );
      // perform logging?
      nonLinLogging_ = params->IsSet( "logging", pdename_, "nonLinear" );
      // type of line search
      params->Get("type", lineSearch_, pdename_, "lineSearch" );
      // incremental stopping criterion
      params->Get("incStopCrit", incStopCrit_, pdename_, "nonLinear" );
      // residual stopping criterion
      params->Get("resStopCrit", residualStopCrit_, pdename_, "nonLinear" );
      // maximal number of NL-iterations
      params->Get("maxNumIters", nonLinMaxIter_, pdename_, "nonLinear");
    }

    // ***************************************************************
    //   If no other damping type is specified and we have absorbing
    //   boundary conditions, then use ABCDAMP
    // ***************************************************************
    absorbingBCs_ = FALSE;
    params->GetList( "name", absBCs_, pdename_, "absorbingBCs" );
    if ( absBCs_.GetSize() ) {
      absorbingBCs_ = TRUE;
      Info->PrintF( pdename_, " Apply Absorbing Boundary Conditions\n" );
      surfdoms_ = absBCs_;
    }

    needsDampingMatrix_ = FALSE;
    // two matrices formulation
    if ( absorbingBCs_ == TRUE || (dampingType_ != NONE  && dampingType_ != FRACTIONAL) )
      needsDampingMatrix_ = TRUE;
    // conservative formulation
    //if ( absorbingBCs_ == TRUE || dampingType_ != NONE )
    //  needsDampingMatrix_ = TRUE;  
  }


  void AcousticPDE::DefineIntegrators(const Integer level) {

    ENTER_FCN( "AcousticPDE::DefineIntegerators" );

    Boolean nonLin = FALSE;
    Double density, compressibility, c0, alpha, beta, BoverA;
    Double coeffmass, coeffdamp;

    for (Integer actSD = 0; actSD < subdoms_.GetSize(); actSD++) {

      density         = materialData_[actSD].GetDensity();
      compressibility = materialData_[actSD].GetCompressibility();
      c0 = sqrt(compressibility/density);

      alpha  = materialData_[actSD].GetDampingAlfa();
      beta   = materialData_[actSD].GetDampingBeta();
      BoverA = materialData_[actSD].GetBoverA();

      // **********************************************************************
      //   Linear wave equation
      // **********************************************************************

      // stiffness integrator
      BaseForm * bilinearStiff = new LaplaceInt(density,isaxi_);        
      IntegratorDescriptor * stiffIntDescr = new IntegratorDescriptor(bilinearStiff, STIFFNESS);
      assemble_->AddIntegrator(stiffIntDescr, subdoms_[actSD]);

      // mass integrator
      coeffmass = density / (c0*c0);
      BaseForm * bilinearMass  = new MassInt(coeffmass, dofspernode_, isaxi_);
      IntegratorDescriptor * massIntDescr = new IntegratorDescriptor(bilinearMass, MASS);
      assemble_->AddIntegrator(massIntDescr, subdoms_[actSD]);

      // **********************************************************************
      //   Additional terms for damping
      // **********************************************************************

      if ( !dampingList_.IsEmpty() ) {

	if (dampingList_[actSD] == "rayleigh") {
	  // This works even after assemble_->AddIntegrator() is executed
	  //   because of the pointers...
	  // stiffness part
	  std::cout << "beta=" << beta << std::endl;

	  stiffIntDescr->SetSecondaryMat(DAMPING, beta, analysistype_);
                
	  // mass part
	  std::cout << "alpha=" << alpha << std::endl;
	  massIntDescr->SetSecondaryMat(DAMPING, alpha, analysistype_);
	}
          
	else if ( dampingList_[actSD] == "thermoViscous" ) {
	  coeffdamp  =  density * 2.0 * alpha * c0;
	  BaseForm * bilinearStiff  = new LaplaceInt(coeffdamp, isaxi_);  
	  IntegratorDescriptor * dampIntDescr = new IntegratorDescriptor(bilinearStiff, DAMPING);
                
	  assemble_->AddIntegrator(dampIntDescr, subdoms_[actSD]);
	}

	else if ( dampingList_[actSD] == "fractional_gl" ) {
	  coeffdamp = - density * 2.0 * alpha / c0 / sin((beta-1.0)*PI/2.0);
	  BaseForm * bilinearDamp  = new MassInt(coeffdamp, dofspernode_, isaxi_);
	  bilinearDamp->SetFracDamping();
	  // conservative formulation
	  // adapt NewmarkFracDamp::Init and BasePDE::GetFracDampMatrixCoeff
	  // IntegratorDescriptor * dampIntDescr = new IntegratorDescriptor(bilinearDamp, DAMPING);

	  // two matrices formulation
	  IntegratorDescriptor * dampIntDescr = new IntegratorDescriptor(bilinearDamp, STIFFNESS);

	  assemble_->AddIntegrator(dampIntDescr, subdoms_[actSD]);
	}

	else if  ( dampingList_[actSD] == "fractional_blank" ) {
	  coeffdamp =  - density * 2.0 * alpha / c0 / sin((beta-1.0)*PI/2.0);
	  coeffdamp *= exp(-gammaln(1.0- (beta- 1.0)) ); // prefactor of blank alg
	  coeffdamp *= 1.0/(1.0- (beta- 1.0));           // weight factor of index 0
	  BaseForm * bilinearDamp  = new MassInt(coeffdamp, dofspernode_, isaxi_);
	  bilinearDamp->SetFracDamping();
	  // conservative formulation
	  // adapt NewmarkFracDamp::Init and BasePDE::GetFracDampMatrixCoeff              
	  // IntegratorDescriptor * dampIntDescr = new IntegratorDescriptor(bilinearDamp, DAMPING);

	  // two matrices formulation
	  IntegratorDescriptor * dampIntDescr = new IntegratorDescriptor(bilinearDamp, STIFFNESS);

	  assemble_->AddIntegrator(dampIntDescr, subdoms_[actSD]);
	}
      }

      //      if ( nonLinPDEType_[actSD] == "kuznetsov" ) {
      //        Double coeffN1 = density * BoverA / pow(c0,4);
      //        BaseForm * N1 = new nLinAcoustic1(coeffN1, isaxi_);
      //        assemble_->AddRhsIntegrator(N1, subdoms_[actSD], nonLin_);
      //        Double coeffN2 = density * 2 / (c0*c0);
      //        BaseForm * N2 = new nLinAcoustic2(coeffN2, isaxi_);
      //        assemble_->AddRhsIntegrator(N2, subdoms_[actSD], nonLin_);
      //      }
      //      else if ( nonLinPDEType_[actSD] == "westervelt" ) {
      //        Double coeffN1 = (1.0+0.5*BoverA) / pow(c0,4);
      //        BaseForm * N1 = new nLinAcoustic1(coeffN1, isaxi_);
      //        assemble_->AddRhsIntegrator(N1, subdoms_[actSD], nonLin_);  
      //      }

    }

    // **********************************************************************
    //   surface-integration: Absorbing boundaries
    // **********************************************************************

    if ( absorbingBCs_ == TRUE && analysistype_ != HARMONIC ) {
      for (Integer actSD = 0; actSD < absBCs_.GetSize(); actSD++) {
	//currently hard-coded!!
	Double density = materialData_[0].GetDensity();
	Double compressibility = materialData_[0].GetCompressibility();
	Double coeffdamp = density/((sqrt(compressibility/density)));
        
	BaseForm * bilinear_damp = new MassInt(coeffdamp,dofspernode_, isaxi_);
	assemble_->AddSurfIntegrator(bilinear_damp,  absBCs_[actSD],
				     DAMPING, nonLin);
      }
    }
 
  }

  void AcousticPDE::DefineSolveStep()
  {
    ENTER_FCN( "AcousticPDE::DefineSolveStep" );
  
    solveStep_ = new SolveStepAcoustic(*this);
  }



  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================

  void AcousticPDE::InitTimeStepping() {
    ENTER_FCN( "AcousticPDE::InitTimeStepping" );

    if ( dampingType_!=FRACTIONAL ) {
      // this includes rayleigh and thermoviscous damping
      if ( effectiveMass_ == FALSE )
	TS_alg_ = new Newmark(pdename_, algsys_, eqnData_, needsDampingMatrix_);
      else
	TS_alg_ = new NewmarkEffMass(pdename_, algsys_, eqnData_, needsDampingMatrix_);
    }
    
    if ( dampingType_ == FRACTIONAL ) {
      if ( effectiveMass_ == FALSE )
	TS_alg_ = new NewmarkFracDamp(pdename_, algsys_, eqnData_, ptgrid_, this,
				      subdoms_, dampingList_, 
				      fracMemory_, inType_, isaxi_);
      else
	; // needs to be implemented
    }

    // Needed for fractional damping model, see Assemble::AssembleMatrices
    assemble_->SetPDEPointer(this);

    // Needed for nonlinear wave equation, get memory for linear part of RHS
    if ( nonLin_ )
      RhsLinVal_.Resize(numPDENodes_);

  }



  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void AcousticPDE::InitCoupling(PDECoupling * Coupling) {
    
    ENTER_FCN( "AcousticPDE::InitCoupling" );
    
    pdeIsCoupled_ = TRUE;
    ptCoupling_   = Coupling;
    
    // Intialize the memory of the coupling values
    for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
      if (ptCoupling_->GetOutputQuantity(i) == ACOU_FORCE)    {
	ptCoupling_->CreateCouplingVector(i,isComplex_);
	
	//now since we need a incremental formulation, initialize some necessary vectors
	isIncrFormulation_ = TRUE;
	solIncr_.Resize(eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN());
	actSol_.Resize(eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN());
	solIncr_.Init(0);
	actSol_.Init(0);
      }
    }
  }
  

  void AcousticPDE::CalcOutputCoupling() {

    ENTER_FCN( "AcousticPDE::CalcOutputCoupling" );

    Integer dof;
    SolutionType quantity;
    StdVector<Elem*> * couplingElems = NULL;
    StdVector<Elem*> * interfaceVolElems = NULL;
    StdVector<Integer> * couplingNodes = NULL;
    StdVector<MaterialData*> * couplingMaterials = NULL;
    CFSVector * temp_values = NULL;
  
    // loop over all output coupling quantities
    for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
      quantity = ptCoupling_->GetOutputQuantity(i);
      ptCoupling_->GetOutputValues(i, temp_values);
      
      Vector<Double> * values = dynamic_cast<Vector<Double>*>(temp_values);

      switch(ptCoupling_->GetOutputType(i)) {

      case NODE:
	if (quantity == ACOU_FORCE) {
	  ptCoupling_->GetOutputElements(i, couplingElems);
	  ptCoupling_->GetOutputNodes(i, couplingNodes);
	  ptCoupling_->GetOwnMaterials(i, couplingMaterials);
	  ptCoupling_->GetInputNeighbourElems(i, interfaceVolElems);
	  dof = ptCoupling_->GetOutputDof(i);
            
	  CalcMechCouplingRHS(couplingElems, *couplingNodes,
			      couplingMaterials, *values, dof,
			      interfaceVolElems);
	}
	break;

      case ELEM:
	Error("No Element coupling output", __FILE__,__LINE__);
      }
    }
  }


  void AcousticPDE::CalcMechCouplingRHS(StdVector<Elem*> * couplingElems, 
					StdVector<Integer> & couplingNodes,
					StdVector<MaterialData*> * couplingMaterials,
					Vector<Double>& elemCouplingSols,
					Integer couplingdof,
					StdVector<Elem*> * interfaceVolElems)
  {
    ENTER_FCN( "AcousticPDE::CalcMechCouplingRHS" );

    Double density=0;
   
    elemCouplingSols.Init(0.0);

    for (Integer actElem=0; actElem<couplingElems->GetSize(); actElem++) {
      Elem * actCoupleElem = (*couplingElems)[actElem];
      BaseFE * ptElem          = actCoupleElem->ptElem;

      StdVector<Integer> connecth = actCoupleElem->connect;

      Matrix<Double> ptCoord; 
      GetElemCoords(connecth, ptCoord, actlevel_);

      Boolean found = FALSE;
      
      // get correct density belonging to the neighbouring element
      // of the interface
      for (Integer actSD = 0; actSD < subdoms_.GetSize(); actSD++)    {  
	if ((*interfaceVolElems)[actElem]->namesd ==  subdoms_[actSD]) {
	  density = materialData_[actSD].GetDensity();
	  found = TRUE;
	}
      }
      
      if (found ==FALSE) {
	std::string msg = "Cannot find correct density to compute acoustic ";
	msg += "pressure forces! Using density of acoustic subdomain 1";
	Info->Warning( msg );
	density = materialData_[0].GetDensity();
      }
          
      BaseForm * bilinear_mass = new MassInt(ptElem, density, isaxi_);
      Matrix<Double> elemmat;
      bilinear_mass->CalcElementMatrix(ptCoord, elemmat);
      delete bilinear_mass;     
      Vector<Double> sol;
      GetDerivSolVecOfElement(sol, connecth);
      
      Vector<Double> forceOnElem = elemmat * sol;
      // force has to be added on RHS with negative sign
      forceOnElem *= -1;

      // the normal vector points outwards of the MECHANICAL domain
      // (see. Kaltenbacher, "Num. Sim. of Mechatr. Act. & Sens." chapter 8.2)
      Vector<Double> n;

      // points outward own domain
      ptgrid_->CalcSurfNormalOutOfVol(n, *actCoupleElem, *(*interfaceVolElems)[actElem]); 
      n *=-1;
      //std::cerr << "acoustic Normal =\n" << n << std::endl;

      for (Integer actNode=0; actNode<ptCoord.GetSizeRow(); actNode++) {
	Integer nodePos = 0;
          
	while(connecth[actNode] != couplingNodes[nodePos] &&
	      nodePos < couplingNodes.GetSize()) {
	  nodePos++;
	}
          
	for (Integer actDof=0; actDof < couplingdof ; actDof++) {
	  elemCouplingSols[nodePos*couplingdof +actDof] += 
	    forceOnElem[actNode] * n[actDof];
	  //std::cerr << "acouforceonElem +=" << forceOnElem[actNode] * n[actDof] << std::endl;
	}
      }
    }
  }


  Boolean AcousticPDE::HasOutput(SolutionType output) {
    ENTER_FCN( "AcousticPDE::HasOutput" );
    if (output == ACOU_FORCE) {
      return TRUE;
    }
    return FALSE;
  }


  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================

  void AcousticPDE::WriteResultsInFile(Integer stepOffset,
				       Double timeOffset) {
    ENTER_FCN( "AcousticPDE::WriteResultsInFile" );

#ifdef PARALLEL //only one thread should write the output
    int commrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
    if (!commrank) {
#endif
      NodeStoreSol<Double> solIm_mesh;
      NodeStoreSol<Double> * solTransient;
      NodeStoreSol<Complex> * solHarmonic;
       
      Double actTime = lasttimecalc_ + timeOffset;
      Integer actStep = laststepcalc_ + stepOffset;
      
      if (analysistype_==HARMONIC) {
	solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);

	if (saveSol_)
	  outFile_->WriteNodeSolutionHarmonic(*solHarmonic,  actFreqStep_, 
					      actFrequency_, complexFormat_);
	if (saveSolHist_)
	  outFile_->WriteNodeHistoryHarmonic(*solHarmonic,  actFreqStep_, 
					     actFrequency_, complexFormat_);
      }
      else {  
	solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);

	if (saveSol_){
	  outFile_->WriteNodeSolutionTransient(*solTransient, actStep, actTime);
          
	  if (saveSolHist_)
	    outFile_->WriteNodeHistoryTransient(*solTransient, actStep, actTime);
	}
        
	if (saveDeriv1_) {
	  solDeriv1_.SetAlgSysVector(getS1()); 
	  outFile_->WriteNodeSolutionTransient(solDeriv1_, actStep, actTime);
        
	  if (saveDeriv1Hist_)
	    outFile_->WriteNodeHistoryTransient(solDeriv1_, actStep, actTime);
	}

	if (saveDeriv2_) {
	  solDeriv2_.SetAlgSysVector(getS2());
	  outFile_->WriteNodeSolutionTransient(solDeriv2_, actStep, actTime);
          
	  if (saveDeriv2Hist_)
	    outFile_->WriteNodeHistoryTransient(solDeriv2_, actStep, actTime);
	}
      }
#ifdef PARALLEL
    }//!commrank
#endif
  }


  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void AcousticPDE::ReadStoreResults() {  
    ENTER_FCN( "AcousticPDE::ReadStoreResults" );
    
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

    if ( formulation_ == ACOU_POTENTIAL ) {
      // --- pressure ---
      Enum2String(ACOU_PRESSURE, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
	Error("It makes no sense to have a PDE in acoustic potential and retrieve \n \
node results in pressure. Try element results!", __FILE__,__LINE__);
      }

      // --- acoustic potential ---
      Enum2String(ACOU_POTENTIAL, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
	saveSol_ = TRUE;
	hasOutput_ = TRUE;
      }

      // --- acoustic potential, 1. Deriv ---
      Enum2String(ACOU_POTENTIAL_DERIV_1, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
	saveDeriv1_ = TRUE;
	hasOutput_ = TRUE;
	// intialize corresponding storesolution object
	solDeriv1_.SetNumSolutions(1);
	solDeriv1_.SetNumNodes(numPDENodes_);
	solDeriv1_.SetSolutionType(ACOU_POTENTIAL_DERIV_1);
	solDeriv1_.SetNumDofs(1);
	solDeriv1_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_); 
	solDeriv1_.Init(0);
      }

      // --- acoustic potential, 2. Deriv ---
      Enum2String(ACOU_POTENTIAL_DERIV_2, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
	saveDeriv2_ = TRUE;
	hasOutput_ = TRUE;
	// intialize corresponding storesolution object
	solDeriv2_.SetNumSolutions(1);
	solDeriv2_.SetNumNodes(numPDENodes_);
	solDeriv2_.SetSolutionType(ACOU_POTENTIAL_DERIV_2);
	solDeriv2_.SetNumDofs(1);
	solDeriv2_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_); 
	solDeriv2_.Init(0);
      }
    }
    else if ( formulation_ == ACOU_PRESSURE ) {
      // --- acoustic potential ---
      Enum2String(ACOU_POTENTIAL, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
	Error("It makes no sense to have a PDE in pressure and retrieve \n \
node results in acoustic potential.", __FILE__,__LINE__);
      }

      // --- pressure ---
      Enum2String(ACOU_PRESSURE, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
	sol_->SetSolutionType(ACOU_PRESSURE);
	sol_->SetNumDofs(1);
	sol_->Init();
	saveSol_ = TRUE;
	hasOutput_ = TRUE;
      }
    }
    
    
    // *****************************
    // Determine element results
    // *****************************
    StdVector<std::string> elementValues;
    keyVec  = pdename_, "storeResults", "elementResults", "region";
    
    // --- much to do here ---

    // *****************************
    // Determine nodal history
    // *****************************
    StdVector<std::string> saveNodeHist; 
    keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
    attrVec = "", "", "type";

    if ( formulation_ == ACOU_POTENTIAL ) {
      // pressure
      Enum2String(ACOU_PRESSURE, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist);
      if (saveNodeHist.GetSize() > 0) {
	Error("It makes no sense to have a PDE in acoustic potential and retrieve \n \
history node results in pressure.", __FILE__,__LINE__);
      }

      // --- acoustic potential  ---
      Enum2String(ACOU_POTENTIAL, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
	saveSolHist_ = TRUE;
	hasOutput_ = TRUE;
	Info->PrintF( pdename_, "Saving acouPotential for Nodes:\n" );
	for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
	  Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
	}
      }
    
      // --- acoustic potential, 1. Deriv ---
      Enum2String(ACOU_POTENTIAL_DERIV_1, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
	saveDeriv1Hist_ = TRUE;
	hasOutput_ = TRUE;
	Info->PrintF( pdename_, "Saving acouPotentialD1 for Nodes:\n" );
	for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
	  Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
	}
      }

      // --- acoustic potential, 1. Deriv ---
      Enum2String(ACOU_POTENTIAL_DERIV_2, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
	saveDeriv1Hist_ = TRUE;
	hasOutput_ = TRUE;
	Info->PrintF( pdename_, "Saving acouPotetentialD2 for Nodes:\n" );
	for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
	  Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
	}
      }
    }

    else if ( formulation_ == ACOU_PRESSURE ) {
      // --- acoustic potential ---
      Enum2String(ACOU_POTENTIAL, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
	Error("It makes no sense to have a PDE in pressure and retrieve \n \
results in acoustic potential.", __FILE__,__LINE__);
      }

      // --- pressure ---
      Enum2String(ACOU_PRESSURE, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
	sol_->SetSolutionType(ACOU_PRESSURE);
	sol_->SetNumDofs(1);
	sol_->Init();
	saveSolHist_ = TRUE;
	hasOutput_ = TRUE;
	Info->PrintF( pdename_, "Saving acouPotential for Nodes:\n" );
	for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
	  Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
	}
      }
    }
  }

} // end of namespace
