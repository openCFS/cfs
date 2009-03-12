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

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

namespace CoupledField
{

  AcousticMixedPDE::AcousticMixedPDE( Grid* aGrid, ParamNode* paramNode )
    :SinglePDE(aGrid, paramNode )

  {

    pdename_          = "acousticMixed";
    pdematerialclass_ = FLUID;
    maxTimeDerivOrder_ = 1;

    // ****************************
    // DETERMINE GEOMETRY
    // ****************************

    // Get problem geometry and PDE subtype
    myParam_->Get( "subType", subType_ );
    std::string probGeo;
    param->Get("domain")->Get("geometryType", probGeo );

    approxType_ = myParam_->Get("type")->AsString();
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
        std::string errmsg = "Subtype '" + subType_;
        errmsg += "' of PDE '" + pdename_ + "' does not fit to problem ";
        errmsg += "geometry '" + probGeo + "'\n";
        Info->Error( errmsg, __FILE__, __LINE__ );
      }


    // timestepping formulation
    std::string str = "";
    myParam_->Get( "timeSteppingFormulation", str,  false );
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

    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( subdoms_[actSD] );
      
      BaseMaterial * actMat = materials_[subdoms_[actSD]];
      actMat->GetScalar(density, DENSITY,Global::REAL);
      actMat->GetScalar(bulkModulus, ACOU_BULK_MODULUS,Global::REAL);
      c0 = sqrt(bulkModulus/density);
      
      Info->PrintF( pdename_, "density = %e\n", density);
      Info->PrintF( pdename_, "bulk modulus = %e\n", bulkModulus);
      
      std::cout << "dim: " << dim_ << std::endl;

      if( approxType_ == "spectral" ){
        VelSurfIntFactor = 1.0;
        abcIntFactor = density * c0;
        // ==============  add stiffness ======================================
        BaseForm *bilinearStiff_KPV = new StiffPiolaMixedInt_KPV( 1.0 , dim_, isaxi_);
        BiLinFormContext * stiffContext_KPV =
          new BiLinFormContext(bilinearStiff_KPV, STIFFNESS );

        stiffContext_KPV->SetPtPdes(this, this);
        stiffContext_KPV->SetResults( results_[0], results_[1],
              actSDList, actSDList );

        assemble_->AddBiLinearForm( stiffContext_KPV );

         // ==============  add KVP ======================================
        BaseForm *bilinearStiff_KVP = new StiffPiolaMixedInt_KVP( 1.0, dim_, isaxi_);
        BiLinFormContext * stiffContext_KVP =
          new BiLinFormContext(bilinearStiff_KVP, STIFFNESS );

        stiffContext_KVP->SetPtPdes(this, this);
        stiffContext_KVP->SetResults( results_[1], results_[0],
              actSDList, actSDList );

        assemble_->AddBiLinearForm( stiffContext_KVP );

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
      else{
        VelSurfIntFactor = bulkModulus;
        abcIntFactor = c0;
        // ==============  add stiffness ======================================
        BaseForm *bilinearStiff_KPV = new StiffMixedInt_KPV( bulkModulus, 
			  				   dim_, isaxi_);
        BiLinFormContext * stiffContext_KPV =
          new BiLinFormContext(bilinearStiff_KPV, STIFFNESS );

        stiffContext_KPV->SetPtPdes(this, this);
        stiffContext_KPV->SetResults( results_[0], results_[1],
			  	    actSDList, actSDList );

        assemble_->AddBiLinearForm( stiffContext_KPV );

         // ==============  add KVP ======================================
        BaseForm *bilinearStiff_KVP = new StiffMixedInt_KVP( 1.0/density, 
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
        surfContext->SetResults( results_[0], results_[0],
         surfVelLHS_[actSD], surfVelLHS_[actSD] );
        assemble_->AddBiLinearForm( surfContext );

        // Give result to equation numbering class
        eqnMap_->AddResult( *results_[0], surfVelLHS_[actSD] );
      }
    }

    // Add integrators for region loads
    VolForceInt * forceInt;
    std::map<RegionIdType, RegionLoad>::iterator loadIt = regionLoads_.begin();
    if (regionLoads_.size() != 0 ) {
      (*loadIt).second.Print(true, pdename_ );
    }
    for( loadIt = regionLoads_.begin(); loadIt != regionLoads_.end(); loadIt++ ) {
      forceInt = (*loadIt).second.GetIntegrator();

      // Create new element list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( loadIt->first );
      LinearFormContext * forceContext = 
        new LinearFormContext( forceInt );
      forceContext->SetPtPde(this);
      forceContext->SetResult( results_[0], actSDList );
      assemble_->AddLinearForm( forceContext );

      (*loadIt).second.Print(false, pdename_);
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
    ParamNode * bcNode = myParam_->Get( "bcsAndLoads", false );
    if( bcNode ) {
      StdVector<ParamNode*> abcNodes = bcNode->GetList( "absorbingBCs" );
      
      for( UInt i = 0; i < abcNodes.GetSize(); i++ ) {
        std::string regionName = abcNodes[i]->Get("name")->AsString(); 
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
    ParamNode * bcNodeVel = myParam_->Get( "bcsAndLoads", false );
    if( bcNodeVel ) {
      StdVector<ParamNode*> surfNodes = bcNodeVel->GetList( "surfVelLHS" );
      
      for( UInt i = 0; i < surfNodes.GetSize(); i++ ) {
        std::string regionName = surfNodes[i]->Get("name")->AsString(); 
        surfVelLHS_.Push_back( ptgrid_->GetEntityList( EntityList::SURF_ELEM_LIST,
                   regionName, EntityList::REGION ) );
        isSurfVelLHS_ = true;
        Info->PrintF( pdename_, 
                      "Apply Surface bilinear form for particle velocity in normal direction '%s'\n",
                      regionName.c_str() );
      }
    }

    //// ***************************************************************
    ////   Surface bilinear form for particle velocity in normal direction
    //// ***************************************************************
    //ParamNode * bcNodeDirichletVel = myParam_->Get( "bcsAndLoads", false );
    //if( bcNodeDirichletVel ) {
    //  StdVector<ParamNode*> velNodes = bcNodeVel->GetList( "dirichletVelCond" );
    //  for( UInt i = 0; i < velNodes.GetSize(); i++ ) {
    //    // read parameters
    //    dof = "";
    //    velNodes[i]->Get( "name", name );
    //    velNodes[i]->Get( "dof1", dof1, false );
    //    velNodes[i]->Get( "dof2", dof2, false );
    //    velNodes[i]->Get( "dof3", dof3, false );
    //    velNodes[i]->Get( "val1", val1, false );
    //    velNodes[i]->Get( "val2", val2, false );
    //    velNodes[i]->Get( "val3", val3, false );
    //    velNodes[i]->Get( "entityType", entType );

    //    surfVelLHS_.Push_back( ptgrid_->GetEntityList( EntityList::SURF_ELEM_LIST,
    //               regionName, EntityList::REGION ) );
    //    isSurfVelLHS_ = true;
    //    Info->PrintF( pdename_, 
    //                  "Apply Surface bilinear form for particle velocity in normal direction '%s'\n",
    //                  regionName.c_str() );
    //  }
    //}
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
    if ( effectiveMass_ == true ) {
      TS_alg_ = new TrapezoidalEffMass( algsys_ );
    }
    else {
      TS_alg_ = new Trapezoidal( algsys_ );
      TS_alg_->SetTrapezoidalGamma(0.51);
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

    std::string approxType = myParam_->Get("type")->AsString();
    if ( approxType == "lagrange" ) {
      shared_ptr<AnsatzFct> fct(new LagrangeFct);
      res1->fctType = fct;
      res1->definedOn = ResultInfo::NODE;
    } 
    else if( approxType == "taylorHood" ) {
      std::cerr << "Using taylorHood!\n";
      UInt order = myParam_->Get("order")->AsUInt();
      shared_ptr<LegendreFct> fct(new LegendreFct);
      fct->SetIsoOrder( order );
      res1->fctType = fct;
      res1->definedOn = ResultInfo::PFEM;
    } 
    else if (  approxType == "spectral" ) {
      UInt order = myParam_->Get("order")->AsUInt();
      shared_ptr<SpectralFct> fct(new SpectralFct);
      res1->definedOn = ResultInfo::PFEM;
      fct->SetOrder( order );
      res1->fctType = fct;
      res1->fctType->SetDiscontinuity(false);
    }
    else {
      EXCEPTION( "approximation type '" << approxType << "' not allowd");
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

    std::string approxTypeV = myParam_->Get("type")->AsString();
    if ( approxTypeV == "lagrange" ) {
      shared_ptr<AnsatzFct> fctV(new LagrangeFct);
      res2->fctType = fctV;
      res2->definedOn = ResultInfo::NODE;
    } 
    else if( approxType == "taylorHood" ) {
      std::cerr << "Using taylorHood!\n";
      UInt order = myParam_->Get("order")->AsUInt();
      shared_ptr<LegendreFct> fctV(new LegendreFct);
      fctV->SetIsoOrder( order+1 );
      res2->fctType = fctV;
      res2->definedOn = ResultInfo::PFEM;
    }
    else if(  approxType == "spectral" ) {
      UInt order = myParam_->Get("order")->AsUInt();
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
      Warning( "Resulttype not computable by mechanic PDE",
               __FILE__, __LINE__ );
    }
  }

  // ************************************************************
  //   PostProcess
  // ************************************************************

  void AcousticMixedPDE::PostProcess() {
  }
} // end namespace CoupledField


