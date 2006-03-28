#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acousticXYZPDE.hh" 
#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "Forms/forms_header.hh"
#include "Forms/abcInt.hh"
#include "Forms/acouNeumannInt.hh"
#include "Estimator/spaceerror.hh"
#include "newmark.hh"
#include "newmarkFracDamp.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/scalarnodeEQN.hh"
#include "Utils/mathfunctions.hh"
#include "Utils/nodestoresol.hh"
#include "Driver/solveStepAcoustic.hh"
#include "CoupledPDE/pdecoupling.hh"

#ifdef TCL_INTERFACE
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

#ifdef PARALLEL
#include <mpi.h>
#endif

namespace CoupledField {



  // =========================================================================
  // set solution information
  // =========================================================================
  AcousticXYZPDE::AcousticXYZPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, 
                           WriteResults *aptOut)
    :SinglePDE(aptgrid,aptOut,aptTimeFunc) {

    ENTER_FCN( "AcousticXYZPDE::AcousticXYZPDE" );

    // ****************************
    // DETERMINE GEOMETRY
    // ****************************

    // Get problem geometry and PDE subtype
    std::string probGeo;
    params->Get( "type", probGeo, "geometry" );

    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry
    if ( probGeo == "3d" ) {
      dofspernode_ = 3;
      Info->PrintF("", "=== 3D PROBLEM\n");
    }
    else if (probGeo == "axi" ) {
      isaxi_ = TRUE;
      dofspernode_ = 2;
      Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
    }
    else if ( probGeo == "plane" ) {
      dofspernode_ = 2;
      Info->PrintF("", "=== PLANE PROBLEM\n");
    }

    solDofs_ = dofspernode_;
    pdename_          = "acousticXYZ";
    pdematerialclass_ = FLUID;

    solTypes_ = ACOU_PRESSUREXYZ;

    // timestepping formulation
    std::string str;
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

  }


  // *********************************************
  //   Check what type of damping should be used
  // *********************************************
  void AcousticXYZPDE::ReadDampingInformation( ) {

    ENTER_FCN( "AcousticXYZPDE::ReadDampingInformation" );

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    for (UInt k = 0; k < subdoms_.GetSize(); k++) {

      // Be aware, that the RegionId entries in subdoms_ need not be in order!
      // We follow the order in subdoms_, which has to be in acordance with
      //  AcousticXYZPDE::DefineIntegrators()

      RegionIdType actRegion = subdoms_[k];
      std::string actRegionName;
      actRegionName = ptgrid_->RegionIdToName( actRegion );
      keyVec = "acoustic" , "region" , "damping" , "type";
      attrVec= ""         , "name"   , "";
      valVec = ""         , actRegionName, "";
      StdVector<std::string> dampInfo;
      params->GetList( keyVec, attrVec, valVec, dampInfo);

      if ( dampInfo.IsEmpty() ) {
        dampingList_[actRegion] = NONE;
        Info->PrintF( pdename_, 
                      "No information specifying damping detected!\n" );
      }
      else if (dampInfo[0] == "no") {
        dampingList_[actRegion] = NONE;
        Info->PrintF( pdename_, 
                      "      * NO damping at all for region: %s\n", 
                      actRegionName.c_str() );;
      }
      else if (dampInfo[0] == "rayleigh") {
        dampingList_[actRegion] = RAYLEIGH;
        Info->PrintF( pdename_, 
                      "      * RAYLEIGH damping for region: %s\n",
                      actRegionName.c_str() );
      }
    }

  }


  void AcousticXYZPDE::DefineIntegrators() {

    ENTER_FCN( "AcousticXYZPDE::DefineIntegerators" );

    // help variables for parameter checking
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    Double density, compressibility, c0, alpha, beta, BoverA;
    Double coeffmass, coeffdamp, coeffstiff;

    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {

      // ********************************************************************
      //   Attention:
      //   In AcousticXYZPDE ALL integrators are multiplied with density!
      // ********************************************************************

      materialData_[actSD]->GetScalar(density,DENSITY,REAL);
      materialData_[actSD]->GetScalar(compressibility,ACOU_BULK_MODULUS,REAL);

      c0 = sqrt(compressibility/density);

      materialData_[actSD]->GetScalar(alpha,RAYLEIGH_ALPHA,REAL);
      materialData_[actSD]->GetScalar(beta,RAYLEIGH_BETA,REAL);
      materialData_[actSD]->GetScalar(BoverA,BOVERA,REAL);

      // ********************************************************************
      //   Linear wave equation
      // ********************************************************************


      //do the stiffness part;
      IntegratorDescriptor * stiffIntDescr;

      //check  for PML!
      RegionIdType actRegion = subdoms_[actSD];
      std::string actRegionName;
      actRegionName = ptgrid_->RegionIdToName( actRegion );
      keyVec = "acousticXYZ" , "region" , "pml" , "type";
      attrVec= ""         , "name"   , "";
      valVec = ""         , actRegionName, "";
      StdVector<std::string> pmlInfo;
      params->GetList( keyVec, attrVec, valVec, pmlInfo);

      if (pmlInfo.GetSize() > 0) {
	//read data for PML layer

	//type of PML damping
	std::string dampingTypePML;

	// inner / outer region
	Matrix<Double> inner;
	Matrix<Double> outer;

	//damping factor
	Double dampPML;

	ReadDataPML(dampingTypePML, inner, dampPML, actRegion);
	dampPML *= c0;

	GetPMLLayerData(inner, outer, actSD);

	//====================================================================
	//	 stiffness integrator for PML
	//====================================================================

	std::string formsType = "stiffness";
	coeffstiff = 1.0 / (c0*c0);
	BaseForm * bilinearStiffPML = 
	  new PMLTransXYZInt(formsType, coeffstiff, dampingTypePML, dampPML, 
			     dofspernode_, isaxi_);

	bilinearStiffPML->SetPosPML(inner,outer);

	IntegratorDescriptor * stiffIntDescrPML = 
	  new IntegratorDescriptor(bilinearStiffPML, STIFFNESS);

	stiffIntDescrPML->SetPDEIds(this, this);
	assemble_->AddIntegrator(stiffIntDescrPML, subdoms_[actSD]);

	//====================================================================
	//	 damping integrator for PML
	//====================================================================
	formsType = "damping";
	coeffdamp = 1.0 / (c0*c0);
	BaseForm * bilinearDampPML = 
	  new PMLTransXYZInt(formsType, coeffdamp, dampingTypePML, dampPML, 
			     dofspernode_, isaxi_);

	bilinearDampPML->SetPosPML(inner,outer);

	IntegratorDescriptor * dampIntDescrPML = 
	  new IntegratorDescriptor(bilinearDampPML, DAMPING);

	dampIntDescrPML->SetPDEIds(this, this);
	assemble_->AddIntegrator(dampIntDescrPML, subdoms_[actSD]);
      }

      //      else {

      // stiffness integrator 
      coeffstiff = 1.0;
      BaseForm * bilinearStiff = new LaplaceXYZInt(coeffstiff,dofspernode_,isaxi_);        
      stiffIntDescr = new IntegratorDescriptor(bilinearStiff, STIFFNESS);
      
      stiffIntDescr->SetPDEIds(this, this);
      
      // mass integrator
      coeffmass = 1.0 / (c0*c0);
      
      BaseForm * bilinearMass  = new MassInt(coeffmass, dofspernode_, isaxi_);
      IntegratorDescriptor * massIntDescr = 
	new IntegratorDescriptor(bilinearMass, MASS);
      
      massIntDescr->SetPDEIds(this, this);
      
      
      // ********************************************************************
      //   Additional terms for damping
      // ********************************************************************
      
      if ( dampingList_.size() > 0 ) {
	(*warning) << "NO Rayleigh damping supported!";
	Warning(__FILE__, __LINE__);  
      }
      
      // Finally add the stiffness/mass integrators
      assemble_->AddIntegrator(stiffIntDescr, subdoms_[actSD]);
      assemble_->AddIntegrator(massIntDescr, subdoms_[actSD]); 
      //      }
    }
  
  }

  void AcousticXYZPDE::DefineSolveStep() {
    ENTER_FCN( "AcousticXYZPDE::DefineSolveStep" );
   
    solveStep_ = new SolveStepAcoustic(*this);
   
  }


  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================

  void AcousticXYZPDE::InitTimeStepping() {
    ENTER_FCN( "AcousticXYZPDE::InitTimeStepping" );

    UInt rhsSize = eqnData_->GetNumEQNs() *
      eqnData_->GetNumDofsPerEQN();

    if ( effectiveMass_ == FALSE ) {
      TS_alg_ = new Newmark( algsys_, rhsSize );
    }
    else {
      TS_alg_ = new NewmarkEffMass( algsys_, rhsSize );
    }
    // Needed for fractional damping model, see Assemble::AssembleMatrices
    assemble_->SetPDEPointer(this);

  }



  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================

  void AcousticXYZPDE::PostProcess() {
    ENTER_FCN( "AcousticXYZPDE::PostProcess" );


    // Last but no least trigger postprocessing fromt within script-file
#ifdef TCL_INTERFACE
    StdVector<std::string> context;
    context.Push_back( pdename_ );
    context.Push_back( GenStr(solveStep_->GetActStep() ) );
    
    if ( analysistype_ == TRANSIENT ||
         analysistype_ == STATIC ) {
      context.Push_back( GenStr(solveStep_->GetActTime() ) );
    } else {
      context.Push_back( GenStr(solveStep_->GetActFreq() ) );
    }
    messenger->TriggerEvent( CFSMessenger::CFS_PostProcess, 
                             context );
#endif   
  }

  void AcousticXYZPDE::WriteResultsInFile(const UInt kstep,
                                       const Double asteptime,
                                       UInt stepOffset,
                                       Double timeOffset) {
    ENTER_FCN( "AcousticXYZPDE::WriteResultsInFile" );

#ifdef PARALLEL //only one thread should write the output
    int commrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
    if (!commrank) {
#endif
      NodeStoreSol<Double> solIm_mesh;
      NodeStoreSol<Double> * solTransient;
      NodeStoreSol<Complex> * solHarmonic;

      NodeStoreSol<Double> * solTransientPress;
      NodeStoreSol<Complex> * solHarmonicPress;

      
      Double actTime = asteptime + timeOffset;
      UInt actStep = kstep + stepOffset;
      
      if (analysistype_==HARMONIC) {
        solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);

        // Please remember, that in harmonic case actTime means indeed
        //  the actual frequency
        if (saveSol_)
          outFile_->WriteNodeSolutionHarmonic(*solHarmonic,  actStep, 
                                              actTime, complexFormat_);
      }
      else {  
        solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);

        if (saveSol_) {
          outFile_->WriteNodeSolutionTransient(*solTransient,actStep,actTime);
        }

	if (savePress_) {
	  //for pressure
	  solTransientPress = dynamic_cast<NodeStoreSol<Double>*>(press_);
	  Double val,valDof;
	  for (UInt node=0; node<numPDENodes_; node++) {
	    val= 0;
	    for (UInt dof=0; dof<dofspernode_; dof++) {
	      solTransient->Get(ACOU_PRESSUREXYZ, node,dof,valDof);
	      val += valDof;
	    }
	    solTransientPress->Set(ACOU_PRESSURE,node,0,val);
	  }
	  outFile_->WriteNodeSolutionTransient(*solTransientPress,actStep,actTime);
	}

        if (saveDeriv1_) {
          solDeriv1_.SetAlgSysVector(getS1()); 
          outFile_->WriteNodeSolutionTransient(solDeriv1_, actStep, actTime);
        }

        if (saveDeriv2_) {
          solDeriv2_.SetAlgSysVector(getS2());
          outFile_->WriteNodeSolutionTransient(solDeriv2_, actStep, actTime);
        }

      }
 
#ifdef PARALLEL
    }//!commrank
#endif
  }

  void AcousticXYZPDE::WriteHistoryInFile(const UInt kstep,
                                       const Double asteptime,
                                       UInt stepOffset,
                                       Double timeOffset) {
    ENTER_FCN( "AcousticXYZPDE::WriteHistoryInFile" );

#ifdef PARALLEL //only one thread should write the output
    int commrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
    if (!commrank) {
#endif
      NodeStoreSol<Double> solIm_mesh;
      NodeStoreSol<Double> * solTransient;
      NodeStoreSol<Complex> * solHarmonic;
      
      NodeStoreSol<Double> * solTransientPress;
      NodeStoreSol<Complex> * solHarmonicPress;

      Double actTime = asteptime + timeOffset;
      UInt actStep = kstep + stepOffset;
      
      if (analysistype_==HARMONIC) {
        solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);

        if (saveSolHist_)
          outFile_->WriteNodeHistoryHarmonic(*solHarmonic,  actStep, 
                                             actTime, complexFormat_);

        if (saveDeriv1Hist_) {
          // multiply solution with j * omega
        }
      }
      else {  
        solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);

        if (saveSolHist_){
          outFile_->WriteNodeHistoryTransient(*solTransient, actStep,actTime);
        }

	if (savePressHist_) {
	  //for pressure
	  solTransientPress = dynamic_cast<NodeStoreSol<Double>*>(press_);
	  Double val,valDof;
	  for (UInt node=0; node<numPDENodes_; node++) {
	    val= 0;
	    for (UInt dof=0; dof<dofspernode_; dof++) {
	      solTransient->Get(ACOU_PRESSUREXYZ, node,dof,valDof);
	      val += valDof;
	    }
	    solTransientPress->Set(ACOU_PRESSURE,node,0,val);
	  }
	  outFile_->WriteNodeHistoryTransient(*solTransientPress,actStep,actTime);
	}

        if (saveDeriv1Hist_) {
          solDeriv1_.SetAlgSysVector(getS1()); 
          outFile_->WriteNodeHistoryTransient(solDeriv1_, actStep, actTime);
        }

        if (saveDeriv2Hist_){
          solDeriv2_.SetAlgSysVector(getS2());
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
  void AcousticXYZPDE::ReadStoreResults() {  
    ENTER_FCN( "AcousticXYZPDE::ReadStoreResults" );
    
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

    // --- pressure ---
    Enum2String(ACOU_PRESSUREXYZ, quantity);
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    saveSol_ = FALSE;
    if (nodeValues.GetSize() > 0) {
      saveSol_ = TRUE;
      hasOutput_ = TRUE;
    }
  
    Enum2String(ACOU_PRESSURE, quantity);
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    savePress_ = FALSE;
    if (nodeValues.GetSize() > 0) {
      savePress_ = TRUE;
      if ( analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC ) {
	press_ = new NodeStoreSol<Complex>;
      }
      else {
	press_ = new NodeStoreSol<Double>;
      }

      press_->SetNumSolutions(1);
      press_->SetSolutionType(ACOU_PRESSURE);
      press_->SetNumNodes(numPDENodes_);
      press_->SetNumDofs(1);
      press_->SetPtrEQNData(eqnData_, ptgrid_);
      press_->Init();
    }
    
    // *****************************
    // Determine element results
    // *****************************
    StdVector<std::string> regionNames;
    keyVec  = pdename_, "storeResults", "elemResults", "region";
    attrVec = "", "", "type";

    // *****************************
    // Determine nodal history
    // *****************************
    StdVector<std::string> saveNodeHist; 
    keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
    attrVec = "", "", "type";

    // --- pressure ---
    Enum2String(ACOU_PRESSUREXYZ, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
    if (saveNodeHist.GetSize() > 0) {
      saveSolHist_ = TRUE;
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, "Saving acouPressureXYZ for Nodes:\n" );
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
	Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
      }
    }

    Enum2String(ACOU_PRESSURE, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
    if (saveNodeHist.GetSize() > 0) {
      savePressHist_ = TRUE;
      Info->PrintF( pdename_, "Saving acouPressure for Nodes:\n" );
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
	Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
      }
    }
    
  }

  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void AcousticXYZPDE::ReadDataPML(std::string& dampingTypePML, Matrix<Double>& inner, 
				Double& dampPML, RegionIdType actRegion) {
  
    ENTER_FCN( "AcousticXYZPDE::ReadDataPML" );

    inner.Resize(dim_,dim_);

    // help variables for parameter checking
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    StdVector<std::string> stringVal;

    // Construct vectors for restricted parameter search
    std::string actRegionName;
    actRegionName = ptgrid_->RegionIdToName( actRegion );
    keyVec = "acousticXYZ" , "region" , "pml" , "type";
    attrVec= ""         , "name"   , "";
    valVec = ""         , actRegionName, "";
    params->GetList( keyVec, attrVec, valVec, stringVal);
    dampingTypePML = stringVal[0];

    StdVector<Double> val;
 
    //xMin
    keyVec = "acousticXYZ" , "region" , "pml" , "xMin";
    params->GetList( keyVec, attrVec, valVec, val);
    inner[0][0] =  val[0];

    //xMax
    keyVec = "acousticXYZ" , "region" , "pml" , "xMax";
    params->GetList( keyVec, attrVec, valVec, val);
    inner[1][0] =  val[0];

    //yMin
    keyVec = "acousticXYZ" , "region" , "pml" , "yMin";
    params->GetList( keyVec, attrVec, valVec, val);
    inner[0][1] =  val[0];

    //yMax
    keyVec = "acousticXYZ" , "region" , "pml" , "yMax";
    params->GetList( keyVec, attrVec, valVec, val);
    inner[1][1] =  val[0];

    if ( dim_ == 3 ) {
      //zMin
      keyVec = "acousticXYZ" , "region" , "pml" , "zMin";
      params->GetList( keyVec, attrVec, valVec, val);
      inner[0][2] =  val[0];

      //zMax     
      keyVec = "acousticXYZ" , "region" , "pml" , "zMax";
      params->GetList( keyVec, attrVec, valVec, val);
      inner[1][2] =  val[0];
    }

    keyVec = "acousticXYZ" , "region" , "pml" , "dampFactor";
    params->GetList( keyVec, attrVec, valVec, val);
    dampPML = val[0];


  }

  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void AcousticXYZPDE::GetPMLLayerData(Matrix<Double>& inner, Matrix<Double>& outer,
				    UInt actSD)  {  

    ENTER_FCN( "AcousticXYZPDE::GetPMLLayerData" );

    outer.Resize(dim_,dim_);
    outer = inner;

    // inner/outer:   xmin  ymin  zmin
    //                xmax  ymax  zmax


    Double minXPML, minYPML, minZPML, maxXPML, maxYPML, maxZPML;

    StdVector<Elem*> elemssd;
    ptgrid_->GetElems(elemssd, subdoms_[actSD]);

    for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
      BaseFE * ptEl = elemssd[actEl]->ptElem;
      StdVector<UInt> connecth = elemssd[actEl]->connect;
                             
      Matrix<Double> ptCoord;
      GetElemCoords(connecth, ptCoord);
      for (UInt i=0; i< ptCoord.GetSizeCol(); i++) {
	//minXPML
	if ( ptCoord[0][i] < outer[0][0] )
	  outer[0][0] = ptCoord[0][i];

	//minYPML
	if ( ptCoord[1][i] < outer[0][1] )
	  outer[0][1] = ptCoord[1][i];

	if (inner.GetSizeCol() > 2 ) {
	  //minZPML
	  if ( ptCoord[2][i] < outer[0][2] )
	    outer[0][2] = ptCoord[2][i];
	}

	//maxXPML
	if ( ptCoord[0][i] > outer[1][0] )
	  outer[1][0] = ptCoord[0][i];

	//maxYPML
	if ( ptCoord[1][i] > outer[1][1] )
	  outer[1][1] = ptCoord[1][i];

	if (inner.GetSizeCol() > 2 ) {
	  //maxZPML
	  if ( ptCoord[2][i] > outer[1][2] )
	    outer[1][2] = ptCoord[2][i];
	}
      }

    }

    //   std::cout << "outer:\n" << outer << std::endl;

  }

} // end of namespace
