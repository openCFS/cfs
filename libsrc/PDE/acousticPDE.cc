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
#include "newmarkdamp.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/scalarnodeEQN.hh"

namespace CoupledField {

  AcousticPDE::AcousticPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc,
			   FileType *aptFileType, WriteResults *aptOut)
    :BasePDE(aptgrid,aptbcs,aptFileType,aptOut,aptTimeFunc) {

    ENTER_FCN( "AcousticPDE::AcousticPDE" );

    // =====================================================================
    // set solution information
    // =====================================================================
    dofspernode_ = 1;
    solTypes_ = ACOU_POTENTIAL;
    solDofs_ = 1;
    pdename_          = "acoustic";
    pdematerialclass_ = "fluid";

    absorbingBCs_ = FALSE;

#ifndef XMLPARAMS
    isaxi_ = FALSE;
    std::string subtype;
    conf->ifget("subtype",subtype,pdename_);
    if (subtype == "axi")
      isaxi_ = TRUE;
#else
    isaxi_ = params->HasValue( "type", "axi", "geometry" );
#endif

    laststepcalc_ = 0;
    dampingType_  = NONE;

#ifndef XMLPARAMS
    std::string dampstr;
    conf->ifget("damping",dampstr,pdename_);
    if (dampstr == "fractional") {
      conf->get( "frac_memory", fracMemory_, pdename_);
      Info->PrintF(pdename_,"\n Attenuation according to power law, number of memory is %g\n\n", fracMemory_);
      dampingType_ = FRACTIONAL;
    }
    else if (dampstr == "rayleigh") {
      dampingType_ = RAYLEIGH;
      Info->PrintF(pdename_, " Using RAYLEIGH damping\n" );
    }

#else

    // *********************************************
    //   Check what type of damping should be used
    // *********************************************

    // Rayleigh damping
    if ( params->HasValue( "type", "rayleigh", pdename_, "damping" )) {
      dampingType_ = RAYLEIGH;
      Info->PrintF( pdename_, "Using RAYLEIGH damping" );
    }

    // Fractional damping
    else if ( params->HasValue( "type", "fractional", pdename_, "damping" )) {
      dampingType_ = FRACTIONAL;
      Info->PrintF( pdename_, "Using FRACTIONAL damping" );
      params->Get( "fracMemory", fracMemory_, pdename_, "damping" );
      std::string msg = "\n Attenuation according to power law, number of ";
      msg += "memory is %g\n\n";
      Info->PrintF( pdename_, msg.c_str(), fracMemory_);
    }

    // Thermoviscous damping
    else if ( params->HasValue( "type", "thermoViscous", pdename_, "damping")){
      Info->PrintF( pdename_, "Using THERMOVISCOUS damping" );
      dampingType_ = ABCDAMP;
    }

    // No damping
    else {
      dampingType_ = NONE;
    }

#endif

    // ***************************************************************
    //   If no other damping type is specified and we have absorbing
    //   boundary conditions, then use ABCDAMP
    // ***************************************************************

#ifndef XMLPARAMS
    conf->ifgetliststr("bnd_for_absBCs",absBCs_,pdename_); 
    if (absBCs_.GetSize() > 0) {
       absorbingBCs_ = TRUE;
       Info->PrintF( pdename_, " Apply Absorbing Boundary Conditions\n" );
    }
    
#else
    params->GetList( "name", absBCs_, pdename_, "absorbingBCs" );
    if ( absBCs_.GetSize() ) {
      absorbingBCs_ = TRUE;
      Info->PrintF( pdename_, " Apply Absorbing Boundary Conditions\n" );
      surfdoms_ = absBCs_;
    }
#endif

    needsDampingMatrix_ = FALSE;
    if ( absorbingBCs_ == TRUE || dampingType_ != NONE ) {
      needsDampingMatrix_ = TRUE;
    }
  
  }


  void AcousticPDE::DefineIntegrators(const Integer level) {

    ENTER_FCN( "AcousticPDE::DefineIntegerators" );

    Boolean nonLin = FALSE;

    for (Integer actSD = 0; actSD < subdoms_.GetSize(); actSD++) {

      Double density = materialData_[actSD].GetDensity();
      Double compressibility = materialData_[actSD].GetCompressibility();

      //stiffness integrator
      BaseForm * bilinearStiff = new LaplaceInt(density,isaxi_);	  
      IntegratorDescriptor * stiffIntDescr = new IntegratorDescriptor(bilinearStiff, STIFFNESS);

      //check for damping (mass part)
      if (dampingType_ == RAYLEIGH)    
      	stiffIntDescr->SetSecondaryMat(DAMPING, materialData_[actSD].GetDampingBeta(), analysistype_);

      assemble_->AddIntegrator(stiffIntDescr, subdoms_[actSD]);


      //mass integrator
      Double coeffmass  = density*density/compressibility;
      BaseForm * bilinearMass  = new MassInt(coeffmass, dofspernode_, isaxi_);

      IntegratorDescriptor * massIntDescr = new IntegratorDescriptor(bilinearMass, MASS);

      //check for damping (mass part)
      if (dampingType_ == RAYLEIGH)    
      	massIntDescr->SetSecondaryMat(DAMPING, materialData_[actSD].GetDampingAlfa(), analysistype_);

      assemble_->AddIntegrator(massIntDescr, subdoms_[actSD]);
    }

    //surface-integration
    // BEGIN DAMPING MATRIX PART: Absorbing boundaries
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

  // ======================================================
  // SOLVING SECTION
  // ======================================================

  void AcousticPDE::InitTimeStepping(const Double dt) {
    ENTER_FCN( "AcousticPDE::InitTimeStepping" );
    
    if ( dampingType_ == FRACTIONAL ) {

      // currently the parameter y is taken from the first subdomain
      // => currently just one subdomain makes sense
      // Double y = materialData_[0].GetDampingBeta();
      // TS_alg_ = new NewmarkDamp(pdename_, algsys_, dofspernode_,
      // numPDENodes_, dampingType_, fracMemory_,y);
    }

    else {
      TS_alg_ = new Newmark(pdename_, algsys_, eqnData_, needsDampingMatrix_);
    }

    TS_alg_->Init(matrix_factor_, dt);

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
      if (ptCoupling_->GetOutputQuantity(i) == ACOU_FORCE)	{
	ptCoupling_->CreateCouplingVector(i,isComplex_);
      }
    }

    iterCoupledCounter_ = 0;
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
      for (Integer actSD = 0; actSD < subdoms_.GetSize(); actSD++)	{  
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
      CalcLineNormalVec(n, *actCoupleElem, *(*interfaceVolElems)[actElem]);

      for (Integer actNode=0; actNode<ptCoord.GetSizeRow(); actNode++) {
	Integer nodePos = 0;
	  
	while(connecth[actNode] != couplingNodes[nodePos] &&
	      nodePos < couplingNodes.GetSize()) {
	  nodePos++;
	}
	  
	for (Integer actDof=0; actDof < couplingdof ; actDof++) {
	  elemCouplingSols[nodePos*couplingdof +actDof] += 
	    forceOnElem[actNode] * n[actDof];
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
#ifdef XMLPARAMS
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
    
    // --- Acou Potential, 1. Deriv ---
    Enum2String(ACOU_POTENTIAL_DERIV_1, quantity);
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveDeriv1_ = TRUE;
      
      // intialize corresponding storesolution object
      solDeriv1_.SetNumSolutions(1);
      solDeriv1_.SetNumNodes(numPDENodes_);
      solDeriv1_.SetSolutionType(ACOU_POTENTIAL_DERIV_1);
      solDeriv1_.SetNumDofs(1);
      solDeriv1_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_); 
      solDeriv1_.Init(0);
    }

    // --- Acou Potential, 2. Deriv ---
    Enum2String(ACOU_POTENTIAL_DERIV_2, quantity);
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveDeriv2_ = TRUE;
      
      // intialize corresponding storesolution object
      solDeriv2_.SetNumSolutions(1);
      solDeriv2_.SetNumNodes(numPDENodes_);
      solDeriv2_.SetSolutionType(ACOU_POTENTIAL_DERIV_2);
      solDeriv2_.SetNumDofs(1);
      solDeriv2_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_); 
      solDeriv2_.Init(0);
    }
    
    // *****************************
    // Determine element results
    // *****************************
    
    // --- nothding to do here ---

    // *****************************
    // Determine nodal history
    // *****************************
    StdVector<std::string> saveNodeHist; 
    keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
    attrVec = "", "", "type";
    
    // --- Acoustic Potential  ---
    Enum2String(ACOU_POTENTIAL, quantity);
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
    Info->PrintF( pdename_, " Saving acouPotential for Nodes:" );
    for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s", saveNodeHist[k].c_str() );
    }
    }
    
    // --- Acoustic Potential, 1. Deriv ---
    Enum2String(ACOU_POTENTIAL_DERIV_1, quantity);
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
      Info->PrintF( pdename_, " Saving acouPotentialD1 for Nodes:" );
      for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
	Info->PrintF( pdename_, " %s", saveNodeHist[k].c_str() );
      }
    }
    // --- Acoustic Potential, 1. Deriv ---
    Enum2String(ACOU_POTENTIAL_DERIV_2, quantity);
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
      Info->PrintF( pdename_, " Saving acouPotetentialD2 for Nodes:" );
      for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
	Info->PrintF( pdename_, " %s", saveNodeHist[k].c_str() );
      }
    }
    
  }
#endif

}
