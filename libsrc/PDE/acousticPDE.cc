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

    dofspernode_      = 1;
    pdename_          = "acoustic";
    pdematerialclass_ = "fluid";

    absorbingBCs_ = FALSE;

#ifndef XMLPARAMS
    isaxi_ = FALSE;
    std::string subtype;
    conf->ifget("subtype",subtype,pdename_);
    if (subtype == "axi")
      isaxi_ = TRUE;
    conf->getsubdompde(subdoms_,pdename_);
#else
    isaxi_ = params->HasValue( "type", "axi", "geometry" );
    params->GetList( "name", subdoms_, pdename_, "region" );
    Info->PrintF( pdename_, " PDE lives on regions:" );
    for ( Integer k = 0; k < subdoms_.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s", subdoms_[k].c_str() );
    }
#endif

    laststepcalc_=0;

    AssignPDENodeNumbers(mesh2PDENode_, pde2MeshNode_, subdoms_);  
    AssignPDEElemNumbers(mesh2PDEElem_, pde2MeshElem_, subdoms_);
    numPDENodes_ = pde2MeshNode_.GetSize();
    numElems_ = pde2MeshElem_.GetSize();

    size_ = numPDENodes_;

    


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

    dampingType_ = NONE;
    
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

    // Absorbing boundary conditions
    if ( params->HasValue( "type", "absorbingBC", pdename_, "damping" )) {
      Info->PrintF( pdename_, "Apply Absorbing Boundary Conditions\n" );
      absorbingBCs_ = TRUE;
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
       Info->PrintF( pdename_, "Apply Absorbing Boundary Conditions\n" );
    }
    
#else
    params->GetList( "name", absBCs_, pdename_, "absorbingBCs" );
    if ( absBCs_.GetSize() > 0 && dampingType_ == NONE ) {
      absorbingBCs_ = TRUE;
      Info->PrintF( pdename_, "Re-setting damping type to ABCDAMP" );
    }
#endif
    
    ReadBCs(pdename_);

    // set analysis parameters
    assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, subdoms_,
				absBCs_);
    assemble_->SetGraphType(NODEGRAPH);
    assemble_->SetMesh2PDENode(&mesh2PDENode_);

#ifdef USE_OLAS
    if (analysistype_==HARMONIC) {
      assemble_->SetMatrixEntryType(OLAS::COMPLEX);
      assemble_->SetMatrixStorageType(OLAS::SPARSE_NONSYM);
    }
    else {
      assemble_->SetMatrixEntryType(OLAS::DOUBLE);
      assemble_->SetMatrixStorageType(OLAS::SPARSE_NONSYM);
    }


#else
    if (analysistype_==HARMONIC) 
      assemble_->SetMatrixType(CSCALAR);
    else
      assemble_->SetMatrixType(RSCALAR);
#endif


    assemble_->SetNumDirichlet(GetNumRestraints(actlevel_));
    assemble_->SetPtrBCs(ptBCs_);
    assemble_->SetPtr2Sol(sol_);
    assemble_->SetPtr2TimeFnc(ptTimeFunc_);

    needsDampingMatrix_ = FALSE;
    if ( absorbingBCs_ == TRUE || dampingType_ == FRACTIONAL ) {
      assemble_->NeedDampingMatrix();
      needsDampingMatrix_ = TRUE;
    }

    ReadMaterialData();
   
    DefineIntegrators(actlevel_);

#ifndef XMLPARAMS
    ReadSavings();
#else
    ReadStoreResults();
#endif
    
    // initialize eqation data object
    eqnData_  = new ScalarNodeEQN(ptgrid_, ptBCs_, subdoms_, actlevel_, dofspernode_);
    eqnData_->SetHomoDirichletBCs(bcs_hd_, homDirichDof_);
    eqnData_->CalcMapping();
    //eqnData_->Print(std::cerr);
    assemble_->SetPtr2EQNData(eqnData_); 

    // Initalize solution class
    sol_->SetNumSolutions(1);
    sol_->SetSolutionType(ACOU_POTENTIAL);
    sol_->SetNumNodes(numPDENodes_);
    sol_->SetNumDofs(dofspernode_);
    sol_->SetPtrEQNData(eqnData_);
    sol_->Init();
    
    if (savederiv1_) {
      sol_der1Array_.SetNumSolutions(1);
      sol_der1Array_.SetNumNodes(numPDENodes_);
      sol_der1Array_.SetSolutionType(ACOU_POTENTIAL_DERIV1);
      sol_der1Array_.SetNumDofs(1);    
      sol_der1Array_.SetPtrEQNData(eqnData_);
      sol_der1Array_.Init(0);
    }

    if (savederiv2_) {
      sol_der2Array_.SetNumSolutions(1);
      sol_der2Array_.SetNumNodes(numPDENodes_);
      sol_der2Array_.SetSolutionType(ACOU_POTENTIAL_DERIV2);
      sol_der2Array_.SetNumDofs(1);    
      sol_der2Array_.SetPtrEQNData(eqnData_);
      sol_der2Array_.Init(0);
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
      if (ptCoupling_->GetOutputQuantity(i) == "acousticforce")	{
	ptCoupling_->CreateCouplingVector(i,isComplex_);
      }
    }

    iterCoupledCounter_ = 0;
  }


  void AcousticPDE::CalcOutputCoupling() {

    ENTER_FCN( "AcousticPDE::CalcOutputCoupling" );

    Integer dof;
    std::string quantity;
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
	  if (quantity == "acousticforce") {
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


  Boolean AcousticPDE::HasOutput(std::string output) {
    ENTER_FCN( "AcousticPDE::HasOutput" );
    if (output == "acousticforce") {
      return TRUE;
    }
    return FALSE;
  }


  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================

  void AcousticPDE::WriteResultsInFile() {
    ENTER_FCN( "AcousticPDE::WriteResultsInFile" );

#ifdef PARALLEL //only one thread should write the output
    int commrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
    if (!commrank) {
#endif
      NodeStoreSol<Double> solIm_mesh;
      
      NodeStoreSol<Double> const & solConverted = 
	dynamic_cast<NodeStoreSol<Double>&>(*sol_);
 
      if (analysistype_==HARMONIC) {
	Error(" Has to be adapted due to new StoreSol class", __FILE__, __LINE__);
	NodeStoreSol<Complex> solmesh;
	Vector<Complex> solution;
	solmesh.GetCompleteVector(solution);
	for (Integer k=0; k<solution.GetSize(); k++)
	  std::cout << "node: " << k+1 << ":  " << solution[k] << std::endl;

	//	outFile_->WriteNodeSolution(sol_mesh,laststepcalc_,lasttimecalc_,"fluid potential, cw realpart,");
	//	outFile_->WriteNodeSolution(solIm_mesh,laststepcalc_,lasttimecalc_,"fluid potential, cw imagpart, ");
      }
      else {  

	//if (savesol_) {
	//  sol_->TransformNodeSolution(sol_mesh, ptgrid_, actlevel_);
	//}
	if (savederiv1_) {
	  sol_der1Array_.SetSolVector(ACOU_VELOCITY,getS1());
	}

	if (savederiv2_) {
	  sol_der2Array_.SetSolVector(ACOU_VELOCITY,getS2());
	  //sol_der2Array_.TransformNodeSolution(solder2_mesh, ptgrid_,actlevel_);
	}
      
	if (outFile_->IsGMV()) {
	  if (savesol_)
	    outFile_->WriteNodeSolution(solConverted,laststepcalc_,lasttimecalc_,"vp");
	  if (savederiv1_)
	    outFile_->WriteNodeSolution(sol_der1Array_,laststepcalc_,
					lasttimecalc_,"vp_1der");
	  if (savederiv2_)
	    outFile_->WriteNodeSolution(sol_der2Array_,laststepcalc_,
					lasttimecalc_,"vp_2der");
	}
	else {
	  if (savesol_)
	    outFile_->WriteNodeSolution(solConverted,laststepcalc_,lasttimecalc_,
					"fluid potential");
	  
	  if (savederiv1_)
	    outFile_->WriteNodeSolution(sol_der1Array_,laststepcalc_,
					lasttimecalc_,
					"fluid potential, 1st deriv.");
	  if (savederiv2_)
	    outFile_->WriteNodeSolution(sol_der2Array_,laststepcalc_,
					lasttimecalc_,
					"fluid potential, 2nd deriv.");
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
    
    // NOTE: This must be changed soon!!!
    // By default we only save the solution at nodal values
    savesol_    = TRUE;
    savederiv1_ = FALSE;
    savederiv2_ = FALSE;

  }
#endif

}
