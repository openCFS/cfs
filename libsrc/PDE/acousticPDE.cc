#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acousticPDE.hh" 
#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "Forms/forms_header.hh"
#include "Forms/abcInt.hh"
#include "Estimator/spaceerror.hh"
#include "newmark.hh"
#include "newmarkFracDamp.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/scalarnodeEQN.hh"
#include "Utils/mathfunctions.hh"
#include "Utils/nodestoresol.hh"
#include "Driver/solveStepAcoustic.hh"
#include "Driver/solveStepAcousticBubble.hh"
#include "CoupledPDE/pdecoupling.hh"

#ifdef PARALLEL
#include <mpi.h>
#endif

namespace CoupledField {

  // =========================================================================
  // set solution information
  // =========================================================================
  AcousticPDE::AcousticPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, 
                           WriteResults *aptOut)
    :SinglePDE(aptgrid,aptOut,aptTimeFunc) {

    ENTER_FCN( "AcousticPDE::AcousticPDE" );

    dofspernode_ = 1;
    solDofs_ = 1;
    pdename_          = "acoustic";
    pdematerialclass_ = "fluid";

    isMechCoupled_ = FALSE;
    saveRHSval_ = FALSE;
    saveRHSvalHist_ = FALSE;

    coarsealpha_ = 0.01; // solver parameter, see basePDE.cc

    // PDE formulation either in acoustic potential or pressure
    std::string str;
    params->Get( "formulation", str, "pdeList", pdename_ );
    String2Enum( str, formulation_ );
    str = "Using * " + str + " as state variable in formulation of PDE\n";
    Info->PrintF( pdename_, str.c_str() );

    // class NodeStoreSol will be initialized with acoustic potential
    solTypes_ = ACOU_POTENTIAL;

    // timestepping formulation
    params->Get( "timeSteppingFormulation", str, "pdeList", pdename_ );
    if ( str == "effMassMatrix" ) {
      effectiveMass_ = TRUE;
      Info->PrintF( pdename_, 
                    "      * effective mass matrix timestepping\n");
    } 
    else {
      effectiveMass_ = FALSE;
      Info->PrintF( pdename_, 
                    "      * effective stiffness matrix timestepping\n");
    }

    // axisymmetric setup    
    isaxi_ = params->HasValue( "type", "axi", "geometry" );


    // *************************************************************
    //   Check what type of nonlinear PDE formulation should be used
    // *************************************************************
    nonLin_ = FALSE; //declaration in basePDE.hh
    StdVector<std::string> strVec;
    params->GetList( "nonLinear", strVec, pdename_, "region" );
    nonLinPDEName_.Resize(strVec.GetSize());

    for ( UInt k = 0; k < strVec.GetSize(); k++ ) {

      if ( strVec[k] == "no" )
        ;
      else if ( strVec[k] == "westervelt" ) {
        nonLin_ = TRUE;
        nonLinPDEName_[k] = WESTERVELT;
        Info->PrintF(pdename_, "      * Westervelt equation for region: %d\n",
					 k );
        if ( solTypes_ == ACOU_POTENTIAL )
          Error("Acoustic potential formulation not supported for \
Westervelt equation!" ,__FILE__,__LINE__);
      }
      else if ( strVec[k] == "kuznetsov" ) {
        nonLin_ = TRUE;
        nonLinPDEName_[k] = KUZNETSOV;
        Info->PrintF(pdename_, "      * Kuznetsov equation for region: %d\n",
					 k );
        if ( solTypes_ == ACOU_PRESSURE )
          Error("Acoustic pressure formulation not supported for \
Kuznetsov equation!" ,__FILE__,__LINE__);
      }
    }

    // **************************************************
    //   Check what type of bubble model should be used
    // **************************************************
    StdVector<std::string> auxVec;
    params->GetList( "bubbleType", auxVec, "acoustic", "bubbles" );
    if ( auxVec.GetSize() == 1 ) {
      String2Enum( auxVec[0], bubbleDynType_ );

      //Set bubbledensity
      params->Get("bubbleNumDensity", bubbleDensity_, "acoustic", "bubbles");
    }
    else if ( auxVec.GetSize() > 1 ) {
      Error("Specification of bubble type not unique in xml-file", __FILE__,
			__LINE__ );
    }
    else {
      bubbleDynType_ = NOBUBBLETYPE;
    }
  

    if( nonLin_ || bubbleDynType_ != NOBUBBLETYPE ) {
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
    params->GetList( "name", auxVec, pdename_, "absorbingBCs" );
    ptgrid_->RegionNameToId( absBCs_, auxVec );

    if ( absBCs_.GetSize() ) {
      matrixTypes_.insert(DAMPING);
	  needsDampingMatrix_ = TRUE;
      absorbingBCs_ = TRUE;
      Info->PrintF( pdename_, " Apply Absorbing Boundary Conditions\n" );
      surfdoms_ = absBCs_;
    }
  }


  // *********************************************
  //   Check what type of damping should be used
  // *********************************************
  void AcousticPDE::ReadDampingInformation( Grid *aptgrid ) {

    ENTER_FCN( "AcousticPDE::ReadDampingInformation" );

    dampingType_ = NONE;
	dampingList_.Resize(subdoms_.GetSize());
	fracMemory_ = 0;
    Boolean identical = TRUE; // i.e. same type of damping for all regions
	Integer firstFrac=-1;

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

	for (UInt k = 0; k < subdoms_.GetSize(); k++) {

	  std::string actRegion;
	  actRegion = aptgrid->RegionIdToName( k );
	  keyVec = "acoustic" , "region" , "damping" , "type";
	  attrVec= ""         , "name"   , "";
	  valVec = ""         , actRegion, "";
	  StdVector<std::string> dampInfo;
	  params->GetList( keyVec, attrVec, valVec, dampInfo);

	  if ( dampInfo.IsEmpty() ) {
		dampingList_[k] = NONE;
		Info->PrintF( pdename_, 
					  "No information specifying damping detected!\n" );
	  }
	  else if (dampInfo[0] == "no") {
		dampingList_[k] = NONE;
		Info->PrintF( pdename_, 
					  "      * NO damping at all for region: %d\n", k );
	  }
	  else if (dampInfo[0] == "rayleigh") {
		dampingList_[k] = RAYLEIGH;
		needsDampingMatrix_ = TRUE;
		Info->PrintF( pdename_, 
					  "      * RAYLEIGH damping for region: %d\n", k );
	  }
	  else if (dampInfo[0] == "thermoViscous") {
		dampingList_[k] = THERMOVISCOUS;
		needsDampingMatrix_ = TRUE;
		Info->PrintF( pdename_, 
					  "      * THERMOVISCOUS damping for region: %d\n", k );
	  }
	  else if (dampInfo[0] == "fractional") {
		dampingType_ = FRACTIONAL;
		Info->PrintF( pdename_, 
					  "      * FRACTIONAL damping for region: %d\n", k );

		// Find first region containing fractional damping
		if ( firstFrac < 0 )
		  firstFrac = k;

		keyVec = "acoustic" , "region" , "damping" , "fracAlg";
		StdVector<std::string> fracAlg;
		params->GetList( keyVec, attrVec, valVec, fracAlg );

		keyVec = "acoustic" , "region" , "damping" , "fracMemory";
		StdVector<UInt> fracMem;
		params->GetList( keyVec, attrVec, valVec, fracMem );

		keyVec = "acoustic" , "region" , "damping" , "interpolation";
		StdVector<std::string> interpol;
		params->GetList(  keyVec, attrVec, valVec, interpol );


		if( fracAlg.IsEmpty() || fracMem.IsEmpty() || interpol.IsEmpty() ) {
		  (*error) << "Specify attributes fracAlg, fracMemory " 
				   << "and interpolation!";
		  Error( __FILE__, __LINE__ ); 
		}
		else if  ( fracAlg[0] == "gl" ) {

		  // Include fracAlg and interpolation info in dampingList
		  Info->PrintF( "", "\t\t\t using Gruenwald-Letnikov algorithm,\n");
		  if (interpol[0] == "no" )
			dampingList_[k] = FRACTIONAL_GL;
		  else {
			dampingList_[k] = FRACTIONAL_GL_INT;
			Info->PrintF("", 
						 "\t\t\t linear interpol. of single past values\n\n");
		  }
		}
		else if ( fracAlg[0] == "blank" ) {

		  Info->PrintF( "", "\t\t\t using Blanks algorithm,\n");
		  if (interpol[0] == "no" )
			dampingList_[k] = FRACTIONAL_BLANK;
		  else {
			dampingList_[k] = FRACTIONAL_BLANK_INT;
			Info->PrintF("", 
						 "\t\t\t linear interpol. of single past values\n\n");
		  }
		}

		// up to now take maximum of fracMemory
		if ( fracMem[0] > fracMemory_ )
		  fracMemory_ = fracMem[0];
	  }

	  // Determine, if regions have different types of damping
	  if ( k > 0 ) {
		if ( dampingList_[k] != dampingList_[k-1] )
		  identical = FALSE;
	  }
	}

	if ( dampingType_ == FRACTIONAL ) {

	  Info->PrintF(pdename_, "Memory size for fractional damping  is: %d\n",
				   fracMemory_ );

	  if ( dampingList_[firstFrac] == FRACTIONAL_GL ||
		   dampingList_[firstFrac] == FRACTIONAL_BLANK )
		inType_ = NOTUSED;
	  else
		inType_ = LIN1PT;
	}

	if ( identical==TRUE && dampingType_!=FRACTIONAL ) {
	  dampingType_ = dampingList_[0];
	}
	else if ( identical==TRUE && dampingType_==FRACTIONAL ) {
	  ;
	}
	else {
	  Info->PrintF(pdename_,
				   "Found different types of damping for regions!\n");
	}
  }

  void AcousticPDE::DefineIntegrators() {

    ENTER_FCN( "AcousticPDE::DefineIntegerators" );

    Boolean nonLin = FALSE;
    Double density, compressibility, c0, alpha, beta, BoverA;
    Double coeffmass, coeffdamp;

    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {

      density         = materialData_[actSD].GetDensity();
      compressibility = materialData_[actSD].GetCompressibility();
      c0 = sqrt(compressibility/density);
      alpha  = materialData_[actSD].GetDampingAlfa();
      beta   = materialData_[actSD].GetDampingBeta();
      BoverA = materialData_[actSD].GetBoverA();

      // if pde couples with mechanic, we have to multiply the density by -1
      if ( isMechCoupled_ == TRUE ) {
        density *= -1.0;
      }

      // ********************************************************************
      //   Linear wave equation
      // ********************************************************************

      // stiffness integrator
      BaseForm * bilinearStiff = new LaplaceInt(density,isaxi_);        
      IntegratorDescriptor * stiffIntDescr = 
        new IntegratorDescriptor(bilinearStiff, STIFFNESS);
      assemble_->AddIntegrator(stiffIntDescr, subdoms_[actSD]);

      // mass integrator
      coeffmass = density / (c0*c0);
#ifdef DEBUG
	  (*debug) << std::endl << "rho/c0^2 = "
			   << coeffmass << std::endl << std::endl;
#endif
      BaseForm * bilinearMass  = new MassInt(coeffmass, dofspernode_, isaxi_);
      IntegratorDescriptor * massIntDescr = 
        new IntegratorDescriptor(bilinearMass, MASS);
      assemble_->AddIntegrator(massIntDescr, subdoms_[actSD]);

      // ********************************************************************
      //   Additional terms for damping
      // ********************************************************************

      if ( !dampingList_.IsEmpty() ) {

		// We check, if damping has been specified for all regions.
		if ( dampingList_.GetSize() != subdoms_.GetSize() ) {
		  (*warning) << "Mismatch between dampingList_ and subdoms_!"
					 << "Size(dampingList_): " << dampingList_.GetSize()
					 << "Size(subdoms_): " << subdoms_.GetSize();
		  Warning(__FILE__, __LINE__);  
		}

        if (dampingList_[actSD] == RAYLEIGH) {
          // This works even after assemble_->AddIntegrator() is executed
          //   because of the pointers...

          // stiffness part
          stiffIntDescr->SetSecondaryMat(DAMPING, beta, analysistype_);
        
          // mass part
          massIntDescr->SetSecondaryMat(DAMPING, alpha, analysistype_);
        }
          
        else if ( dampingList_[actSD] == THERMOVISCOUS ) {
          coeffdamp  =  density * 2.0 * alpha * c0;
          BaseForm * bilinearStiff  = new LaplaceInt(coeffdamp, isaxi_);  
          IntegratorDescriptor * dampIntDescr = 
            new IntegratorDescriptor(bilinearStiff, DAMPING);
                
          assemble_->AddIntegrator(dampIntDescr, subdoms_[actSD]);
        }

        else if ( dampingList_[actSD] == FRACTIONAL_GL ||
				  dampingList_[actSD] == FRACTIONAL_GL_INT ) {
          coeffdamp = - density * 2.0 * alpha / c0 / sin((beta-1.0)*PI/2.0);
#ifdef DEBUG
		  (*debug) << std::endl << "-rho*2*alpha0/c0/sin((y-1)*0.5*pi) = "
				   << coeffdamp << std::endl << std::endl;
#endif
          BaseForm * bilinearDamp  = 
			new MassInt(coeffdamp, dofspernode_, isaxi_);
          bilinearDamp->SetFracDamping();
          // formulation using DAMPING matrix
          // adapt NewmarkFracDamp::Init and StdPDE::GetFracDampMatrixCoeff
          // IntegratorDescriptor * dampIntDescr = 
          //   new IntegratorDescriptor(bilinearDamp, DAMPING);

          // two matrices formulation
		  // added to STIFFNESS matrix because, because 
		  //   matrix_factors[STIFFNESS] = 1.0
          IntegratorDescriptor * dampIntDescr = 
            new IntegratorDescriptor(bilinearDamp, STIFFNESS);

          assemble_->AddIntegrator(dampIntDescr, subdoms_[actSD]);
        }

        else if  ( dampingList_[actSD] == FRACTIONAL_BLANK ||
				   dampingList_[actSD] == FRACTIONAL_BLANK_INT ) {
          coeffdamp =  - density * 2.0 * alpha / c0 / sin((beta-1.0)*PI/2.0);
          // prefactor of blank alg
          coeffdamp *= exp(-gammaln(1.0- (beta- 1.0)) ); 
          // weight factor of index 0
          coeffdamp *= 1.0/(1.0- (beta- 1.0));           
          BaseForm * bilinearDamp  = 
            new MassInt(coeffdamp, dofspernode_, isaxi_);
          bilinearDamp->SetFracDamping();
          // formulation using DAMPING matrix
          // adapt NewmarkFracDamp::Init and StdPDE::GetFracDampMatrixCoeff
          // IntegratorDescriptor * dampIntDescr = 
          //   new IntegratorDescriptor(bilinearDamp, DAMPING);

          // two matrices formulation
		  // added to STIFFNESS matrix because, because 
		  //   matrix_factors[STIFFNESS] = 1.0
          IntegratorDescriptor * dampIntDescr = 
            new IntegratorDescriptor(bilinearDamp, STIFFNESS);

          assemble_->AddIntegrator(dampIntDescr, subdoms_[actSD]);
        }
      }
    }

    // **********************************************************************
    //   surface-integration: Absorbing boundaries
    // **********************************************************************

    if ( absorbingBCs_ == TRUE && analysistype_ != HARMONIC ) {
      for (UInt actSD = 0; actSD < absBCs_.GetSize(); actSD++) {
        SurfForm * bilinear_damp = new AbsorbingBCsInt(isaxi_);
        bilinear_damp->SetFirstVoluInfo(pdename_, subdoms_, materialData_);

        // In the case of acou-mech coupling we have to multiply the 
        // abc-Integrator matrix with -1
        if ( isMechCoupled_ == TRUE ) {
          bilinear_damp->SetFactor(-1.0);
        }
          
        assemble_->AddSurfIntegrator(bilinear_damp,  absBCs_[actSD],
                                     DAMPING, nonLin);
      }
    }
  }

  void AcousticPDE::DefineSolveStep() {
    ENTER_FCN( "AcousticPDE::DefineSolveStep" );

    if( bubbleDynType_ != NOBUBBLETYPE ) {
      solveStep_ = new SolveStepAcousticBubble(*this, bubbleDynType_);
    }
    else {
      solveStep_ = new SolveStepAcoustic(*this);
    }
  }



  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================

  void AcousticPDE::InitTimeStepping() {
    ENTER_FCN( "AcousticPDE::InitTimeStepping" );

    UInt rhsSize = eqnData_->GetNumEQNs() *
      eqnData_->GetNumDofsPerEQN();

	// this includes rayleigh and thermoviscous damping
    if ( dampingType_!=FRACTIONAL ) { 
      if ( effectiveMass_ == FALSE ) {
        TS_alg_ = new Newmark( algsys_, rhsSize );
      }
      else {
        TS_alg_ = new NewmarkEffMass( algsys_, rhsSize );
      }
    }
    else {
      if ( effectiveMass_ == FALSE )
        TS_alg_ = new NewmarkFracDamp( algsys_, rhsSize, 
									   pdeId_, eqnData_, ptgrid_, this, 
									   subdoms_, dampingList_ );
      else
        Error("This needs to be implemented!",__FILE__,__LINE__);
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
    
    isIterCoupled_ = TRUE;
    ptCoupling_   = Coupling;
    
    // Intialize the memory of the coupling values
    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
      if (ptCoupling_->GetOutputQuantity(i) == ACOU_FORCE)    {
        ptCoupling_->CreateCouplingVector(i,isComplex_);
        
        // now since we need a incremental formulation, 
		//  initialize some necessary vectors
        isIncrFormulation_ = TRUE;
      }
    }
  }
  

  void AcousticPDE::CalcOutputCoupling() {

    ENTER_FCN( "AcousticPDE::CalcOutputCoupling" );

    UInt dof;
    SolutionType quantity;
    StdVector<Elem*> * couplingElems = NULL;
    StdVector<UInt> * couplingNodes = NULL;
    CFSVector * temp_values = NULL;
  
    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == FALSE)
      return;

    // loop over all output coupling quantities
    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
      quantity = ptCoupling_->GetOutputQuantity(i);
      ptCoupling_->GetOutputValues(i, temp_values);

      // hard coded cast, since coupling is only possible with
      // real valued entries
      Vector<Double> * values = dynamic_cast<Vector<Double>*>(temp_values);

      switch(ptCoupling_->GetOutputType(i)) {

      case NODE:
        if (quantity == ACOU_FORCE) {
          ptCoupling_->GetOutputElements(i, couplingElems);
          ptCoupling_->GetOutputNodes(i, couplingNodes);
          dof = ptCoupling_->GetOutputDof(i);
          

          CalcMechCouplingRHS(couplingElems, *couplingNodes,
                              *values, dof);                              
        }
        break;

      case ELEM:
        Error("No Element coupling output", __FILE__,__LINE__);
      }
    }
  }


  void AcousticPDE::
  CalcMechCouplingRHS( StdVector<Elem*> * couplingElems, 
                       StdVector<UInt> & couplingNodes,
                       Vector<Double>& elemCouplingSols,
                       UInt couplingdof ) {
    
    ENTER_FCN( "AcousticPDE::CalcMechCouplingRHS" );
    
    Double density = 0.0;
    Double sign = 0.0;
    Integer matIndex = -1;
    Elem * ptVolElem = NULL;
    Matrix<Double> ptCoord, elemMat;
    Vector<Double> sol, forceOnElem, normal;
    
    elemCouplingSols.Init(0.0);
    
    for (UInt actElem=0; actElem<couplingElems->GetSize(); actElem++) {
      
      // Perform cast from volume element to surface element, since
      // mech-acou coupling makes only sense on surface elements
      SurfElem * actCoupleElem = 
        dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);
      
      BaseFE * ptElem = actCoupleElem->ptElem;
      StdVector<UInt> & connecth = actCoupleElem->connect;
      GetElemCoords(connecth, ptCoord);
      
      // Try to find according region for first neighbouring volume
      // element of the surface element
      matIndex = subdoms_.Find(actCoupleElem->ptVolElem1->regionId);
      
      // If first volume element does not belong to acoustic PDE, try the
      // second one
      if ( matIndex == -1 ) {
        matIndex = subdoms_.Find(actCoupleElem->ptVolElem2->regionId);
        ptVolElem = actCoupleElem->ptVolElem2;
        sign = actCoupleElem->normalSign;
      } else {
        ptVolElem = actCoupleElem->ptVolElem1;
        sign = -1.0 * actCoupleElem->normalSign;
      }
      
      if ( matIndex == -1) {
        (*error) << "AcousticPDE::CalcMechCouplingRHS: The two volume "
                 << "element neighbours of surface element Nr. "
                 << actCoupleElem->elemNum << " do not belong to my regions!";
        Error( __FILE__, __LINE__ );
      }
      
      // Assign correct density
      density = materialData_[matIndex].GetDensity();
      
      
      
      BaseForm * bilinear_mass = new MassInt(ptElem, density, isaxi_);
      bilinear_mass->CalcElementMatrix(ptCoord, elemMat);
      delete bilinear_mass;     
      
      GetDerivSolVecOfElement(sol, connecth);
      
      forceOnElem = elemMat * sol;
      
      // force has to be added on RHS with negative sign
      forceOnElem *= - 1.0;
      
      // the normal vector points outwards of the MECHANICAL domain
      // (see. Kaltenbacher, "Num. Sim. of Mechatr. Act. & Sens." chapter 8.2)
      // ptgrid_->CalcSurfNormalOutOfVol(n, *actCoupleElem, 
      // *(*interfaceVolElems)[actElem]); 
      ptgrid_->CalcSurfNormal(normal, *actCoupleElem);
      normal *= sign;

      for (UInt actNode=0; actNode<ptCoord.GetSizeRow(); actNode++) {
        UInt nodePos = 0;
          
        while(connecth[actNode] != couplingNodes[nodePos] &&
              nodePos < couplingNodes.GetSize()) {
          nodePos++;
        }
          
        for (UInt actDof=0; actDof < couplingdof ; actDof++) {
          elemCouplingSols[nodePos*couplingdof +actDof] += 
            forceOnElem[actNode] * normal[actDof];
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

  void AcousticPDE::WriteResultsInFile(const UInt kstep,
                                       const Double asteptime,
                                       UInt stepOffset,
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
      
      
      Double actTime = asteptime + timeOffset;
      UInt actStep = kstep + stepOffset;
      
      if (analysistype_==HARMONIC) {
        solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);

        if (saveSol_)
          outFile_->WriteNodeSolutionHarmonic(*solHarmonic,  actStep, 
                                              actTime, complexFormat_);
        if (saveSolHist_)
          outFile_->WriteNodeHistoryHarmonic(*solHarmonic,  actStep, 
                                             actTime, complexFormat_);
      }
      else {  
        solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);

        if (saveSol_){
          outFile_->WriteNodeSolutionTransient(*solTransient,actStep,actTime);
         
        }
        if (saveSolHist_){
          outFile_->WriteNodeHistoryTransient(*solTransient, actStep,actTime);
        }
        
        if (saveDeriv1_) {
          solDeriv1_.SetAlgSysVector(getS1()); 
          outFile_->WriteNodeSolutionTransient(solDeriv1_, actStep, actTime);

        }
        
        if (saveDeriv1Hist_) {
          solDeriv1_.SetAlgSysVector(getS1()); 
          outFile_->WriteNodeHistoryTransient(solDeriv1_, actStep, actTime);
        }

        if (saveDeriv2_) {
          solDeriv2_.SetAlgSysVector(getS2());
          outFile_->WriteNodeSolutionTransient(solDeriv2_, actStep, actTime);
        }
        if (saveDeriv2Hist_){
          solDeriv2_.SetAlgSysVector(getS2());
          outFile_->WriteNodeHistoryTransient(solDeriv2_, actStep, actTime);
        }
        if (saveRHSval_){
          outFile_->WriteNodeSolutionTransient(rhs_, actStep, actTime);
        }
        if (saveRHSvalHist_){
          outFile_->WriteNodeHistoryTransient(rhs_, actStep, actTime); 
        }
      }

      // ------- for bubble results ----------------------
      if(bubbleDynType_ != NOBUBBLETYPE &&
         (saveSol_ == TRUE || saveDeriv1_ == TRUE || saveDeriv2_ == TRUE)) {

        StdVector<Double> radius = dynamic_cast<SolveStepAcousticBubble*>
          (solveStep_)->GetResultData("bubbleRadius");

        StdVector<Double> velocity = dynamic_cast<SolveStepAcousticBubble*>
          (solveStep_)->GetResultData("bubbleVelocity");

        ElemStoreSol<Double> bubbleResult;
    
        // Resize solution arrays
        bubbleResult.SetNumSolutions(1);
        bubbleResult.SetSolutionType(ELEC_FIELD_INTENSITY);
        bubbleResult.SetNumElems(numElems_);
          
        // dimension hard coded for .unv file!
        bubbleResult.SetNumDofs(3);  
        bubbleResult.SetPtrEQNData(eqnData_, ptgrid_);
        bubbleResult.Init();
          
        for (UInt el=0; el<numElems_; el++) {
          Vector<Double> result(3);
            
          result[0] = radius[el];
          result[1] = velocity[el];
          result[2] = (4.0/3.0)*PI*bubbleDensity_*radius[el]*radius[el]
			*radius[el];
          //        if (el == 90)
          //std::cerr<<actTime<<"   " <<el<< "   " << result[0] << "   " 
          // result[1] << "     " << result[2] << std::endl;
            
          bubbleResult.SetElemResult(el,result);
        }
          
        outFile_->WriteElemSolutionTransient(bubbleResult, actStep, actTime);
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
        (*error) << "It makes no sense to have a PDE in acoustic potential "
                 << "and retrieve \n node results in pressure. "
                 << "Try element results!";
        Error( __FILE__, __LINE__ );
      }

      // --- acoustic potential ---
      Enum2String(ACOU_POTENTIAL, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
        saveSol_ = TRUE;
        hasOutput_ = TRUE;
      }

      // --- acoustic rhsval ---
      Enum2String(ACOU_RHSVAL, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
	saveRHSval_ = TRUE;
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
        solDeriv1_.SetPtrEQNData(eqnData_, ptgrid_); 
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
        solDeriv2_.SetPtrEQNData(eqnData_, ptgrid_); 
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
        std::string warnMsg;
        warnMsg = "Due to the restrictions in the .unv file format, the ";
        warnMsg += "acoustic pressure is written as acoustic (fluid) ";
        warnMsg += "potential!";
        Warning(warnMsg.c_str(), __FILE__, __LINE__);

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
        (*error) <<  "It makes no sense to have a PDE in acoustic potential "
                 << " and retrieve \n history node results in pressure.";
        Error( __FILE__, __LINE__ );
      }

      // --- acoustic potential  ---
      Enum2String(ACOU_POTENTIAL, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
        saveSolHist_ = TRUE;
        hasOutput_ = TRUE;
        Info->PrintF( pdename_, "Saving acouPotential for Nodes:\n" );
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
        }

      }


      // --- acoustic RHS ---
      Enum2String(ACOU_RHSVAL, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
        saveRHSvalHist_ = TRUE;
        hasOutput_ = TRUE;
        Info->PrintF( pdename_, "Saving acouRHSval for Nodes:\n" );
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
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
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
        }
        solDeriv1_.SetNumSolutions(1);
        solDeriv1_.SetNumNodes(numPDENodes_);
        solDeriv1_.SetSolutionType(ACOU_POTENTIAL_DERIV_1);
        solDeriv1_.SetNumDofs(1);
        solDeriv1_.SetPtrEQNData(eqnData_, ptgrid_); 
        solDeriv1_.Init(0);
      }

      // --- acoustic potential, 2. Deriv ---
      Enum2String(ACOU_POTENTIAL_DERIV_2, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
        saveDeriv2Hist_ = TRUE;
        hasOutput_ = TRUE;
        Info->PrintF( pdename_, "Saving acouPotetentialD2 for Nodes:\n" );
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
        }
        solDeriv2_.SetNumSolutions(1);
        solDeriv2_.SetNumNodes(numPDENodes_);
        solDeriv2_.SetSolutionType(ACOU_POTENTIAL_DERIV_2);
        solDeriv2_.SetNumDofs(1);
        solDeriv2_.SetPtrEQNData(eqnData_, ptgrid_); 
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
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
        }
      }
    }
  }


} // end of namespace
