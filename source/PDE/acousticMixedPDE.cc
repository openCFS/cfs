#include "acousticMixedPDE.hh"
#include <sstream>
#include <iomanip>

#include "Forms/mixedInt.hh"
#include "Forms/linElastInt.hh"
#include "Forms/massInt.hh"
#include "Forms/linPressureInt.hh"
#include "Driver/assemble.hh"
#include "newmark.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"
#include "Driver/assemble.hh"
#include "trapezoidal.hh"
#include "bdf2.hh"
#include "StdPDE.hh"
#include "Driver/stdSolveStep.hh"
#include "PDE/mixedEqnMap.hh"
#include "Forms/PMLInt.hh"

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

namespace CoupledField
{

  AcousticMixedPDE::AcousticMixedPDE( Grid* aGrid, PtrParamNode paramNode )
    :SinglePDE(aGrid, paramNode )

  {

    pdename_          = "acousticMixed";
    pdematerialclass_ = FLUID;
    maxTimeDerivOrder_ = 1;

    // ****************************
    // DETERMINE GEOMETRY
    // ****************************

    // Get problem geometry and PDE subtype
    myParam_->GetValue( "subType", subType_ );
    std::string probGeo;
    param->Get("domain")->GetValue("geometryType", probGeo );

    approxType_ = myParam_->Get("type")->As<std::string>();
    if(approxType_ == "spectral"){
      std::cerr << std::endl << "Usng the mixed Equation Map" << std::endl;
      //eqnMap_->~EqnMap();
      eqnMap_ = shared_ptr<MixedEqnMap>(new MixedEqnMap( ptgrid_, pdeId_, usePenalty_ ));
    }


    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry
    UInt dofspernode = 0;
    if ( subType_ == "3d" && probGeo == "3d" ) {
      dofspernode = 4;
      Info->PrintF("", "=== 3D PROBLEM\n");
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = true;
      dofspernode = 3;
      Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
    }
    else if ( subType_ == "plane" && probGeo == "plane" ) {
      isaxi_ = false;
      dofspernode = 3;
      Info->PrintF("", "=== PLANE PROBLEM\n");
    }
    else
    {
      EXCEPTION("Subtype '" << subType_ << "' of PDE '" << pdename_ 
                <<"' does not fit to problem geometry '" << probGeo);
    }


    // timestepping formulation
    std::string str = "";
    myParam_->GetValue( "timeSteppingFormulation", str,  ParamNode::PASS );
    //myParam_->Get( "timeSteppingFormulation", str,  true );
    //str = "effMassMatrix";
    if ( str == "effMassMatrix" ) {
      effectiveMass_ = true;
      Info->PrintF( pdename_, 
                    "      * effective mass matrix timestepping\n");
    } 
    else if ( str == "diagMassMatrix" ) {
      diagMass_      = true;
      effectiveMass_ = true;
      Info->PrintF( pdename_, 
                    "      * diagonal mass matrix in explicit timestepping\n");
    }//just preparing for explicit schemes which are implemented later
    else if (str == "leapfrog"){
      diagMass_  = true;
      Info->PrintF( pdename_, 
                    "      * explicit timestepping using leapfrog second order scheme\n");
    }
    else if (str == "RK4"){
      diagMass_  = true;
      Info->PrintF( pdename_, 
                    "      * explicit timestepping using runge-kutta fourth order scheme\n");
    }
    else {
      effectiveMass_ = false;
      Info->PrintF( pdename_, 
                    "      * effective stiffness matrix timestepping\n");
    }
   }

  AcousticMixedPDE::~AcousticMixedPDE()
  {
  }

  void AcousticMixedPDE::DefineIntegrators()
  {

    // help variables for parameter checking
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    
    Double density, bulkModulus,c0;
    Double VelSurfIntFactor = 0;
    Double abcIntFactor = 0;
    Double MachNumber = 1.0;

    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( subdoms_[actSD] );
      RegionIdType actRegion = subdoms_[actSD];
      std::string actRegionName;
      
      actRegionName = ptgrid_->GetRegion().ToString( actRegion );
      BaseMaterial * actMat = materials_[subdoms_[actSD]];
      actMat->GetScalar(density, DENSITY,Global::REAL);
      actMat->GetScalar(bulkModulus, ACOU_BULK_MODULUS,Global::REAL);
      c0 = sqrt(bulkModulus/density);
      
      PtrParamNode actRegionNode = myParam_->Get("regionList")
          ->GetByVal( "region", "name", actRegionName );

      Info->PrintF( pdename_, "density = %e\n", density);
      Info->PrintF( pdename_, "bulk modulus = %e\n", bulkModulus);
      
      if( approxType_ == "spectral" ){
        VelSurfIntFactor = 1.0;
        abcIntFactor = density * c0;
        // Check for Perfectly matchec layers
        if ( dampingList_[actRegion] == PML ) {
          //read data for PML layer

          //type of PML damping
          std::string dampingTypePML;
          
          // inner / outer region
          Matrix<Double> inner;
          Matrix<Double> outer;
          std::string coordSysId;
          
          //damping factor
          Double dampPML;
          
          std::string id = actRegionNode->Get("dampingId")->As<std::string>();
          PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal   ("pml", "id", id);
          ReadDataPML(dampingTypePML, inner, dampPML, coordSysId, pmlNode );
          dampPML *= c0;
          
          GetPMLLayerData(inner, outer, actRegion, coordSysId);
          if ( analysistype_ == HARMONIC ) {        
            //====================================================================
            //	 mass integrator for PML PP
            //====================================================================
            
            std::string formsType = "mixedPMLMassPPInt";
            
            //set real part
            BaseForm * bilinearStiff_pml_PP =
              new PMLMixedInt(formsType, 1.0/bulkModulus, dampingTypePML, dampPML, isaxi_);
            
            bilinearStiff_pml_PP->SetPosPML(inner,outer,coordSysId);
            
            BiLinFormContext * stiffContextMPP =
              new BiLinFormContext( bilinearStiff_pml_PP, MASS );
            
            stiffContextMPP->SetPtPdes(this, this);
            stiffContextMPP->SetResults( results_[0], results_[0],
                                          actSDList, actSDList );
            // stiffContextReal->SetEntryType(matType);
            assemble_->AddBiLinearForm( stiffContextMPP);
            
            //====================================================================
            //	 mass integrator for PML VV
            //====================================================================
            
            formsType = "mixedPMLMassVVInt";
            
            //set real part
            BaseForm * bilinearStiff_pml_VV =
              new PMLMixedInt(formsType, density, dampingTypePML, dampPML, isaxi_);
            
            bilinearStiff_pml_VV->SetPosPML(inner,outer,coordSysId);
            
            BiLinFormContext * stiffContextMVV =
              new BiLinFormContext( bilinearStiff_pml_VV , MASS );
            
            stiffContextMVV->SetPtPdes(this, this);
            stiffContextMVV->SetResults( results_[1], results_[1],
                                          actSDList, actSDList );
            // stiffContextReal->SetEntryType(matType);
            assemble_->AddBiLinearForm( stiffContextMVV);
            
          } else{ // end of pml part
            //EXCEPTION("PML is only available in Frequency domain");
            //Define additional Integrators i.e.
            // 1. M_phi,phi Mass integrator for the PML Variable
            // 2. R_phi^sigma Integrator for \nabla v in sigma equation
            // 3. C_V,V Damped Stiffness integrator for velocity
            // 4. C_phi,phi damped Stiffness integrator for phi
            // 5. M_P,phi Integrator for Sigma in euler EQ

            std::string formsType = "NoInt";

            //================================================
            // 1. M_phi,phi Mass integrator for the PML Variable
            // ===============================================

              //formsType = "PMLAuxVecMass";

              //BaseForm * bilinearM_phi =
              //   new PMLMixedTimeInt(formsType, 1.0 , dampingTypePML, dampPML, isaxi_);

              //bilinearM_phi->SetPosPML(inner,outer);

              //BiLinFormContext * bilinearM_phiContext =
              //              new BiLinFormContext( bilinearM_phi, MASS );
              //bilinearM_phiContext->SetResults( results_[2], results_[2],
              //                               actSDList, actSDList );
              //bilinearM_phiContext->SetPtPdes(this, this);
              //assemble_->AddBiLinearForm( bilinearM_phiContext);
              
              BaseForm * bilinearMass_PhiPhi = new MassMixedInt_VV(1.0, dim_, isaxi_);

              BiLinFormContext * massContext_PhiPhi =  
                new BiLinFormContext( bilinearMass_PhiPhi, MASS );

              massContext_PhiPhi->SetPtPdes(this, this);
              massContext_PhiPhi->SetResults( results_[2], results_[2],
                     actSDList, actSDList );
              assemble_->AddBiLinearForm( massContext_PhiPhi );

            //================================================
            // 2. 2. R_Q^sigma Integrator 
            //================================================

              formsType = "PMLGradR_PhiSigma";

              BaseForm * PMLGradR_PhiSigma =
                 new PMLMixedTimeInt(formsType, -1.0 , dampingTypePML, dampPML, isaxi_);

              PMLGradR_PhiSigma->SetPosPML(inner,outer, coordSysId);

              BiLinFormContext * PMLGradR_PhiSigmaContext =
                            new BiLinFormContext( PMLGradR_PhiSigma, STIFFNESS );
              PMLGradR_PhiSigmaContext->SetResults( results_[2], results_[1],
                                             actSDList, actSDList );
              PMLGradR_PhiSigmaContext->SetPtPdes(this, this);
              assemble_->AddBiLinearForm( PMLGradR_PhiSigmaContext);

            //================================================
            // 3. C_VV Damped Stiffness integrator for velocity
            //===============================================

              formsType = "PMLVelStiff";

              BaseForm * bilinearC_VV =
                 new PMLMixedTimeInt(formsType, 1.0 * density , dampingTypePML, dampPML, isaxi_);

              bilinearC_VV->SetPosPML(inner,outer, coordSysId);

              BiLinFormContext * bilinearC_VVContext =
                            new BiLinFormContext( bilinearC_VV, STIFFNESS );
              bilinearC_VVContext->SetResults( results_[1], results_[1],
                                             actSDList, actSDList );
              bilinearC_VVContext->SetPtPdes(this, this);
              assemble_->AddBiLinearForm( bilinearC_VVContext);

            //================================================
            // 4. C_QQ damped Stiffness integrator for Q
            //================================================

              formsType = "PMLAccelStiff";

              BaseForm * bilinearC_phi =
                 new PMLMixedTimeInt(formsType, 1.0 , dampingTypePML, dampPML, isaxi_);

              bilinearC_phi->SetPosPML(inner,outer, coordSysId);

              BiLinFormContext * bilinearC_phiContext =
                            new BiLinFormContext( bilinearC_phi, STIFFNESS );
              bilinearC_phiContext->SetResults( results_[2], results_[2],
                                             actSDList, actSDList );
              bilinearC_phiContext->SetPtPdes(this, this);
              assemble_->AddBiLinearForm( bilinearC_phiContext);

            //================================================
            // 5. M_P,Q Integrator for Q in euler EQ
            //================================================

              formsType = "PMLStiffPhi";

              BaseForm * PMLMassPhi =
                   new PMLMixedTimeInt(formsType, 1.0 , dampingTypePML, dampPML, isaxi_);

              PMLMassPhi->SetPosPML(inner,outer,coordSysId);

              BiLinFormContext * PMLMassPhiContext =
                            new BiLinFormContext( PMLMassPhi, STIFFNESS );
              PMLMassPhiContext->SetResults( results_[0], results_[2],
                                             actSDList, actSDList );
              PMLMassPhiContext->SetPtPdes(this, this);
              assemble_->AddBiLinearForm( PMLMassPhiContext);

            // Give result to equation numbering class
            eqnMap_->AddResult( *results_[2], actSDList );


            //DEFINE STANDARD MASS INTEGRATORS FOR P AND V
            //==============  add MPP ======================================
            BaseForm * bilinearMass_PP = new MassMixedInt_PP(1.0/bulkModulus, isaxi_);


            BiLinFormContext * massContext_PP =  
            	new BiLinFormContext( bilinearMass_PP, MASS );

            massContext_PP->SetPtPdes(this, this);
            massContext_PP->SetResults( results_[0], results_[0],
			             actSDList, actSDList );
            assemble_->AddBiLinearForm( massContext_PP );


            //==============  add MVV ======================================  
            BaseForm * bilinearMass_VV = new MassPiolaMixedInt_VV(density, dim_, isaxi_);

            BiLinFormContext * massContext_VV =  
              new BiLinFormContext( bilinearMass_VV, MASS );

            massContext_VV->SetPtPdes(this, this);
            massContext_VV->SetResults( results_[1], results_[1],
                   actSDList, actSDList );
            assemble_->AddBiLinearForm( massContext_VV );
          }
        }else{
          //==============  add MPP ======================================
          BaseForm * bilinearMass_PP = new MassMixedInt_PP(1.0/bulkModulus, isaxi_);

          BiLinFormContext * massContext_PP =  
          	new BiLinFormContext( bilinearMass_PP, MASS );

          massContext_PP->SetPtPdes(this, this);
          massContext_PP->SetResults( results_[0], results_[0],
			           actSDList, actSDList );
          assemble_->AddBiLinearForm( massContext_PP );


          //==============  add MVV ======================================  
          BaseForm * bilinearMass_VV = new MassPiolaMixedInt_VV(density, dim_, isaxi_);

          BiLinFormContext * massContext_VV =  
            new BiLinFormContext( bilinearMass_VV, MASS );

          massContext_VV->SetPtPdes(this, this);
          massContext_VV->SetResults( results_[1], results_[1],
                 actSDList, actSDList );
          assemble_->AddBiLinearForm( massContext_VV );
        }
        Double kpvFactor = -1.0;
        Double kvpFactor = 1.0;

        // ==============  add stiffness ======================================
        BaseForm *bilinearStiff_KPV = new StiffPiolaMixedInt_KPV( kpvFactor , dim_, isaxi_);
        BiLinFormContext * stiffContext_KPV =
          new BiLinFormContext(bilinearStiff_KPV, STIFFNESS );

        stiffContext_KPV->SetPtPdes(this, this);
        stiffContext_KPV->SetResults( results_[0], results_[1],
              actSDList, actSDList );

        assemble_->AddBiLinearForm( stiffContext_KPV );

         // ==============  add KVP ======================================
        BaseForm *bilinearStiff_KVP = new StiffPiolaMixedInt_KVP( kvpFactor, dim_, isaxi_);
        BiLinFormContext * stiffContext_KVP =
          new BiLinFormContext(bilinearStiff_KVP, STIFFNESS );

        stiffContext_KVP->SetPtPdes(this, this);
        stiffContext_KVP->SetResults( results_[1], results_[0],
              actSDList, actSDList );

        assemble_->AddBiLinearForm( stiffContext_KVP );

      }
      else{
        VelSurfIntFactor = bulkModulus;
        abcIntFactor = c0;


	      // check for flow data

	      // get current region name and get grip of paramNode
	      RegionIdType actRegion = subdoms_[actSD];
	      std::string actRegionName;
	      actRegionName = ptgrid_->GetRegion().ToString(actRegion);
	      PtrParamNode actRegionNode =
	        myParam_->Get("regionList")->GetByVal( "region", "name", actRegionName );

              std::string id = actRegionNode->Get("flowId")->As<std::string>();
	      if ( id != "" ) {
	        // add the convective bilinear forms

	        PtrParamNode flowNode = myParam_->Get("flowList")->GetByVal("flow", "id", id);

	        Vector<Double> velVec;
	        Double MachNr;
	        ReadDataFlow(velVec, MachNr, flowNode );
	        
	        MachNumber = MachNr;
	        VelSurfIntFactor /= MachNumber;
	        abcIntFactor /= MachNumber;

	        // ==============  add convective PV ======================================
	        BaseForm *bilinearConv_KPP = new ConvectiveMixedInt_KPP( velVec, 
	      							   dim_, isaxi_);
	        BiLinFormContext * convContext_KPP =
	          new BiLinFormContext(bilinearConv_KPP, STIFFNESS );

	        convContext_KPP->SetPtPdes(this, this);
	        convContext_KPP->SetResults( results_[0], results_[0],
	      			       actSDList, actSDList );

	        assemble_->AddBiLinearForm( convContext_KPP );


	        // ==============  add convective VP ======================================
	        BaseForm *bilinearConv_KVV = new ConvectiveMixedInt_KVV( velVec, 
	      							   dim_, isaxi_);
	        BiLinFormContext * convContext_KVV =
	          new BiLinFormContext(bilinearConv_KVV, STIFFNESS );

	        convContext_KVV->SetPtPdes(this, this);
	        convContext_KVV->SetResults( results_[1], results_[1],
	      			       actSDList, actSDList );

	        assemble_->AddBiLinearForm( convContext_KVV );
	      }

        // ==============  add stiffness ======================================
        BaseForm *bilinearStiff_KPV = new StiffMixedInt_KPV( bulkModulus/MachNumber, 
			  				   dim_, isaxi_);
        BiLinFormContext * stiffContext_KPV =
          new BiLinFormContext(bilinearStiff_KPV, STIFFNESS );

        stiffContext_KPV->SetPtPdes(this, this);
        stiffContext_KPV->SetResults( results_[0], results_[1],
			  	    actSDList, actSDList );

        assemble_->AddBiLinearForm( stiffContext_KPV );

         // ==============  add KVP ======================================
        BaseForm *bilinearStiff_KVP = new StiffMixedInt_KVP( 1.0/(density*MachNumber), 
			  				   dim_, isaxi_);
        BiLinFormContext * stiffContext_KVP =
          new BiLinFormContext(bilinearStiff_KVP, STIFFNESS );

        stiffContext_KVP->SetPtPdes(this, this);
        stiffContext_KVP->SetResults( results_[1], results_[0],
			  	    actSDList, actSDList );

        assemble_->AddBiLinearForm( stiffContext_KVP );

        //==============  add MPP ======================================
        Double coeffMass = 1.0;
        BaseForm * bilinearMass_PP = new MassMixedInt_PP(coeffMass, isaxi_);

        BiLinFormContext * massContext_PP =  
        	new BiLinFormContext( bilinearMass_PP, MASS );

        massContext_PP->SetPtPdes(this, this);
        massContext_PP->SetResults( results_[0], results_[0],
			         actSDList, actSDList );
        assemble_->AddBiLinearForm( massContext_PP );


        //==============  add MVV ======================================	
        BaseForm * bilinearMass_VV = new MassMixedInt_VV(coeffMass, dim_, isaxi_);

        BiLinFormContext * massContext_VV =  
        	new BiLinFormContext( bilinearMass_VV, MASS );

        massContext_VV->SetPtPdes(this, this);
        massContext_VV->SetResults( results_[1], results_[1],
			         actSDList, actSDList );
        assemble_->AddBiLinearForm( massContext_VV );
      }

      // Give result to equation numbering class
      eqnMap_->AddResult( *results_[0], actSDList );
      eqnMap_->AddResult( *results_[1], actSDList );

    }
    
    // *******************************************************************************
    //   inhom. Neumann boundary condition: acutually given normal surface velocity!!
    // *******************************************************************************
    for( UInt iBc = 0; iBc < inBcs_.GetSize(); iBc++ ) {
      
      // get current Bc
      InhomNeumannBc const & actBc = *inBcs_[iBc];
      LinearSurfForm *neumannBC = new LinSurfVelocity( actBc.value, actBc.phase,
                                                     VelSurfIntFactor, isaxi_ );
      LinearFormContext * neumannContext = new LinearFormContext( neumannBC );
      neumannContext->SetPtPde( this );
      neumannContext->SetResult( results_[0], actBc.entities );
      assemble_->AddLinearForm( neumannContext );
      
      // Give result to equation numbering class
      eqnMap_->AddResult( *results_[0], actBc.entities );
    }


    // **********************************************************************
    //   surface-integration: Absorbing boundaries
    // **********************************************************************
    if ( absorbingBCs_ == true) { 
      for (UInt actSD = 0; actSD < absBCs_.GetSize(); actSD++) {
        ABC_MixedInt * bilinear_abc = new ABC_MixedInt(abcIntFactor,isaxi_);
        BiLinFormContext * abcContext = 
          new BiLinFormContext( bilinear_abc, STIFFNESS );
        abcContext->SetPtPdes(this, this);     
        abcContext->SetResults( results_[0], results_[0],
                                absBCs_[actSD], absBCs_[actSD] );
        assemble_->AddBiLinearForm( abcContext );

        // Give result to equation numbering class
        eqnMap_->AddResult( *results_[0], absBCs_[actSD] );
      }
    }

    // **********************************************************************
    //   surface-bilinear form for normal particle velocity
    // **********************************************************************
    if ( isSurfVelLHS_ == true) { 
      for (UInt actSD = 0; actSD < surfVelLHS_.GetSize(); actSD++) {
        SurfVel_MixedInt * bilinear_surf = new SurfVel_MixedInt( VelSurfIntFactor,
                 dim_, isaxi_);
        BiLinFormContext * surfContext = 
          new BiLinFormContext( bilinear_surf, STIFFNESS );
        surfContext->SetPtPdes(this, this);     
        surfContext->SetResults( results_[0], results_[1],
         surfVelLHS_[actSD], surfVelLHS_[actSD] );
        assemble_->AddBiLinearForm( surfContext );

        // Give result to equation numbering class
        eqnMap_->AddResult( *results_[0], surfVelLHS_[actSD] );
        eqnMap_->AddResult( *results_[1], surfVelLHS_[actSD] );
      }
    }

    // Add integrators for region loads
    VolForceInt * forceInt;
    std::map<RegionIdType, RegionLoad>::iterator loadIt = regionLoads_.begin();
    for( loadIt = regionLoads_.begin(); loadIt != regionLoads_.end(); loadIt++ ) {
      forceInt = (*loadIt).second.GetIntegrator();

      // Create new element list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( loadIt->first );
      LinearFormContext * forceContext = 
        new LinearFormContext( forceInt );
      forceContext->SetPtPde(this);
      forceContext->SetResult( results_[1], actSDList );
      assemble_->AddLinearForm( forceContext );

      (*loadIt).second.ToInfo(infoNode_->Get("regionLoad"));
    }
  }


  void AcousticMixedPDE::ReadSpecialBCs() {
    
    // read volume force definition
    ReadRegionLoads();

    // ***************************************************************
    //   Abbsorbing BCs
    // ***************************************************************
    StdVector<std::string> auxVec;
    absorbingBCs_ = false;
    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    if( bcNode ) {
      ParamNodeList abcNodes = bcNode->GetList( "absorbingBCs" );
      
      for( UInt i = 0; i < abcNodes.GetSize(); i++ ) {
        std::string regionName = abcNodes[i]->Get("name")->As<std::string>(); 
        absBCs_.Push_back( ptgrid_
          ->GetEntityList( EntityList::SURF_ELEM_LIST,
                           regionName, EntityList::REGION ) );
        absorbingBCs_ = true;
        Info->PrintF( pdename_, 
                      "Apply Absorbing Boundary Conditions on surfaceRegion '%s'\n",
                      regionName.c_str() );
      }
    }

    // ***************************************************************
    //   Surface bilinear form for particle velocity in normal direction
    // ***************************************************************
    isSurfVelLHS_ = false;
    PtrParamNode bcNodeVel = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    if( bcNodeVel ) {
      ParamNodeList surfNodes = bcNodeVel->GetList( "surfVelLHS" );
      
      for( UInt i = 0; i < surfNodes.GetSize(); i++ ) {
        std::string regionName = surfNodes[i]->Get("name")->As<std::string>(); 
        surfVelLHS_.Push_back( ptgrid_->GetEntityList( EntityList::SURF_ELEM_LIST,
                   regionName, EntityList::REGION ) );
        isSurfVelLHS_ = true;
        Info->PrintF( pdename_, 
                      "Apply Surface bilinear form for particle velocity in normal direction '%s'\n",
                      regionName.c_str() );
      }
    }

  }
  
  void AcousticMixedPDE::SetSpecialBCs(){
    return;
  }
  
  void AcousticMixedPDE::DefineSolveStep()
  {

    solveStep_ = new  StdSolveStep(*this);
  }


  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================


  void AcousticMixedPDE::InitTimeStepping()
  {
    PtrParamNode systemNode = FindLinearSystem(pdename_);
    if ( effectiveMass_ == true ) {
      TS_alg_ = new TrapezoidalEffMass( algsys_, systemNode );
    }
    else {
      TS_alg_ = new Trapezoidal( algsys_, systemNode );
      TS_alg_->SetTrapezoidalGamma(0.505);
      //      TS_alg_ = new Bdf2( algsys_ );
    }
  }



  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================


  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void AcousticMixedPDE::DefineAvailResults() {

    // =====================================================================
    // set solution information
    // =====================================================================

    // === PRESSURE ===
    shared_ptr<ResultInfo> res1(new ResultInfo);
    res1->resultType = ACOU_PRESSURE;
    res1->dofNames = "";
    res1->unit = "Pa";
    res1->entryType = ResultInfo::SCALAR;

    std::string approxType = myParam_->Get("type")->As<std::string>();
    if ( approxType == "lagrange" ) {
      shared_ptr<AnsatzFct> fct(new LagrangeFct);
      res1->fctType = fct;
      res1->definedOn = ResultInfo::NODE;
    } 
    else if( approxType == "taylorHood" ) {
      std::cerr << "Using taylorHood!\n";
      UInt order = myParam_->Get("order")->As<UInt>();
      shared_ptr<LegendreFct> fct(new LegendreFct);
      fct->SetIsoOrder( order );
      res1->fctType = fct;
      res1->definedOn = ResultInfo::PFEM;
    } 
    else if (  approxType == "spectral" ) {
      UInt order = myParam_->Get("order")->As<UInt>();
      shared_ptr<SpectralFct> fct(new SpectralFct);
      res1->definedOn = ResultInfo::PFEM;
      fct->SetOrder( order );
      res1->fctType = fct;
      res1->fctType->SetDiscontinuity(false);
    }
    else {
      EXCEPTION( "approximation type '" << approxType << "' not allowed");
    }
      
    results_.Push_back( res1 );
    availResults_.insert( res1 );

    // === VELOCITY ===
    shared_ptr<ResultInfo> res2(new ResultInfo); res2->resultType = ACOU_VELOCITY;
    if( subType_ == "3d" ) {
      res2->dofNames = "x", "y", "z";
    } 
    else if ( subType_ == "axi" ) {
      res2->dofNames = "r", "z";
    } 
    else if( subType_ == "plane") {
      res2->dofNames = "x", "y";
    } 
    else {
      EXCEPTION( "AcousticMixedPDE: subType '" << subType_ 
                   << "' not known'" );
    }

    res2->unit = "m/s";
    res2->entryType = ResultInfo::VECTOR;

    std::string approxTypeV = myParam_->Get("type")->As<std::string>();
    if ( approxTypeV == "lagrange" ) {
      shared_ptr<AnsatzFct> fctV(new LagrangeFct);
      res2->fctType = fctV;
      res2->definedOn = ResultInfo::NODE;
    } 
    else if( approxType == "taylorHood" ) {
      std::cerr << "Using taylorHood!\n";
      UInt order = myParam_->Get("order")->As<UInt>();
      shared_ptr<LegendreFct> fctV(new LegendreFct);
      fctV->SetIsoOrder( order+1 );
      res2->fctType = fctV;
      res2->definedOn = ResultInfo::PFEM;
    }
    else if(  approxType == "spectral" ) {
      UInt order = myParam_->Get("order")->As<UInt>();
      shared_ptr<SpectralFct> fct(new SpectralFct);
      res2->definedOn = ResultInfo::PFEM;
      fct->SetOrder( order );
      res2->fctType = fct;
      res2->fctType->SetDiscontinuity(true);
    }
    else {
      EXCEPTION( "approximation type '" << approxType << "' not allowd");
    }

    results_.Push_back( res2 );
    availResults_.insert( res2 );

    // for time domain PML we need auxillary variables
    // check for time domain PML
    if ( analysistype_ == TRANSIENT ) { 
      bool isPML = false;

      for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {
        // get current region name and get grip of paramNode
        RegionIdType actRegion = subdoms_[actSD];
        // Check for Perfectly matchec layers
        if ( dampingList_[actRegion] == PML )
          isPML = true;
      }

      if ( isPML ) {
        //======================================================================
        //AXILIARY VARIABLES FOR TIME DOMAIN PML
        //=====================================================================
        // === U ===
        shared_ptr<ResultInfo> res4(new ResultInfo); 
        res4->resultType = ACOU_ACCELERATION;
        if( subType_ == "3d" ) {
          res4->dofNames = "x", "y", "z";
        } 
        else if ( subType_ == "axi" ) {
          res4->dofNames = "r", "z";
        } 
        else if( subType_ == "plane") {
          res4->dofNames = "x", "y";
        } 
        else {
          EXCEPTION( "AcousticMixedPDE: subType '" << subType_ 
                       << "' not known'" );
        }

        res4->unit = "m/s^2";
        res4->entryType = ResultInfo::VECTOR;

        if ( approxType == "lagrange" ) {
          shared_ptr<AnsatzFct> fctU(new LagrangeFct);
          res4->fctType = fctU;
          res4->definedOn = ResultInfo::NODE;
        } 
        else if( approxType == "taylorHood" ) {
          std::cerr << "Using taylorHood!\n";
          UInt order = myParam_->Get("order")->As<UInt>();
          shared_ptr<LegendreFct> fctU(new LegendreFct);
          fctU->SetIsoOrder( order+1 );
          res4->fctType = fctU;
          res4->definedOn = ResultInfo::PFEM;
        }
        else if(  approxType == "spectral" ) {
          UInt order = myParam_->Get("order")->As<UInt>();
          shared_ptr<SpectralFct> fctU(new SpectralFct);
          res4->definedOn = ResultInfo::PFEM;
          fctU->SetOrder( order );
          res4->fctType = fctU;
          res4->fctType->SetDiscontinuity(false);
        }
        else {
          EXCEPTION( "approximation type '" << approxType << "' not allowed");
        }

        results_.Push_back( res4 );
        availResults_.insert( res4 );
      }
    }
    // ===  RHS LOAD ===
    //shared_ptr<ResultInfo> rhs(new ResultInfo);
    //rhs->resultType = ACOU_RHS_LOAD;
    //rhs->dofNames = dispDofNames;
    //rhs->unit = "N";
    //rhs->entryType = disp->entryType;
    //rhs->definedOn = disp->definedOn;
    //rhs->fctType = disp->fctType;
    //availResults_.insert( rhs );
    
  }


  void  AcousticMixedPDE::CalcResults( shared_ptr<BaseResult> result ) {
    
    switch (result->GetResultInfo()->resultType ) {
      
    case ACOU_PRESSURE:
      if( isComplex_ ) {
        ExtractResult<Complex>( result, sol_ );
      } else {
        ExtractResult<Double>( result, sol_ );
      }
      break;

    case ACOU_ACCELERATION:
      if( isComplex_ ) {
        ExtractResult<Complex>( result, sol_ );
      } else {
        ExtractResult<Double>( result, sol_ );
      }
      break;

    case ACOU_VELOCITY:
      if( isComplex_ ) {
        ExtractResult<Complex>( result, sol_ );
      } else {
        ExtractResult<Double>( result, sol_ );
      }
      break;

    case ACOU_RHS_LOAD:
      if( isComplex_ ) {
        ExtractRhsResult<Complex>( result, results_[0] );
      } else {
        ExtractRhsResult<Double>( result, results_[0] );
      }
      break;
      
    default:
      WARN( "Resulttype not computable by mechanic PDE" );
    }
  }

  // *********************************************
  //   Check what type of damping should be used
  // *********************************************
  void AcousticMixedPDE::ReadDampingInformation( ) {


    fracMemory_ = 0;
    std::map<std::string, DampingType> idDampType;
    std::map<std::string, Double>      idDampFreq;
    std::map<std::string, Double>      idDampRatioDeltaF;

    // try to get dampingList
    PtrParamNode dampListNode = myParam_->Get( "dampingList", ParamNode::PASS );
    if( dampListNode ) {

      // get specific damping nodes
      ParamNodeList dampNodes = dampListNode->GetChildren();

      for( UInt i = 0; i < dampNodes.GetSize(); i++ ) {

        std::string dampString = dampNodes[i]->GetName();
        std::string actId = dampNodes[i]->Get("id")->As<std::string>();

        // determine type of damping
        DampingType actType;
        String2Enum( dampString, actType );

        // store damping type string
        idDampType[actId] = actType;

      }
    }

    // Run over all region and set entry in "regionNonLinId"
    ParamNodeList regionNodes =
      myParam_->Get("regionList")->GetChildren();

    RegionIdType actRegionId;
    std::string actRegionName, actDampingId;

    if( regionNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Damping in following region(s)\n" );
    }

    for (UInt k = 0; k < regionNodes.GetSize(); k++) {
      regionNodes[k]->GetValue( "name", actRegionName );
      regionNodes[k]->GetValue( "dampingId", actDampingId );
      if( actDampingId == "" )
        continue;

      actRegionId = ptgrid_->GetRegion().Parse(actRegionName);

      // Check actDampingId was already registerd
      if( idDampType.count( actDampingId ) == 0 ) {
        EXCEPTION( "Damping with id '" << actDampingId
                   << "' was not defined in 'dampingList'" );
      }

      dampingList_[actRegionId] = idDampType[actDampingId];

      // Log to info file
      std::string dampString;
      Enum2String( dampingList_[actRegionId], dampString );

      Info->PrintF( pdename_, " %s: %s\n", actRegionName.c_str(),
                    dampString.c_str() );
    }
    Info->PrintF( pdename_, "\n" );
  }

  //   Obtain information about flow data
  // ***********************************************************************
  void AcousticMixedPDE::ReadDataFlow(Vector<Double>& velVec, Double& MachVal,
				      PtrParamNode actNode ) {


    velVec.Resize(dim_);
    // Check, if pml node has a child "propRegion"
    PtrParamNode velDataNode = actNode->Get( "velocityData", ParamNode::PASS );

    //Vx
    velDataNode->GetValue( "Vx", velVec[0] );

    //Vy
    velDataNode->GetValue( "Vy", velVec[1] );

    if ( dim_ == 3) {
      //Vz
      velDataNode->GetValue("Vz", velVec[2] );
    }

    //get Mach number
    actNode->GetValue( "Mach", MachVal );

  }


  // ************************************************************
  //   PostProcess
  // ************************************************************

  void AcousticMixedPDE::PostProcess() {
  }
} // end namespace CoupledField


