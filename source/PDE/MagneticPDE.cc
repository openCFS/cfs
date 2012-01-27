// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MagneticPDE.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/Coil.hh"

// include fespaces
#include "FeBasis/H1/FeSpaceH1.hh"
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/DivOperator.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
#include "Forms/Operators/IdentityOperator.hh"

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "CoupledPDE/PDECoupling.hh"

namespace CoupledField {

DECLARE_LOG(magpde)
DEFINE_LOG(magpde, "magpde")

MagneticPDE::MagneticPDE(Grid * aptgrid, PtrParamNode paramNode )
    :SinglePDE( aptgrid, paramNode ) {

  // =====================================================================
  // set solution information
  // =====================================================================
  pdename_          = "magnetic";
  pdematerialclass_ = ELECTROMAGNETIC;
  maxTimeDerivOrder_ = 1;

  // check for use of mixed formulation, i.e. for transient / harmonic
}

  MagneticPDE::~MagneticPDE()
  {

  }


  void MagneticPDE::InitNonLin()
  {

    SinglePDE::InitNonLin();
  }


  void MagneticPDE::DefineIntegrators() {

    // determine tensor representation of the material parameters needed
    SubTensorType tensorType;
    if ( dim_ == 3 ) {
      tensorType = FULL;
    } else {
      if ( isaxi_ == true ) {
        tensorType = AXI;
      } else {
        // 2d: plane case
        tensorType = PLANE_STRAIN;
      }
    }
    
    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    
    
    //  Loop over all regions
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
      
      // Set current region and material
      actRegion = it->first;
      actMat = it->second;

      // Get current region name
      std::string regionName = ptgrid_->GetRegion().ToString(actRegion);

      // create new entity list and add it fefunction
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );
      myFct->AddEntityList( actSDList );
      
      // --------------------------
      //  Set region approximation
      // --------------------------
      // --- Set the approximation for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);
      
      
      // ====================================================================
      //  Standard Linear CASE (2D AND 3D)
      // ====================================================================
      if( !nonLin_ ) {
        BaseBDBInt * stiffInt = NULL;
        shared_ptr<CoefFunction> curCoef = 
            actMat->GetCoefFunction(MAG_RELUCTIVITY,tensorType,Global::REAL);
        if( dim_ == 2) {
          if( isaxi_ ) {
            // axisymmetric case
            stiffInt = new BDBInt<CurlOperatorAxi<Double> >(curCoef, 1.0);
          } else {
            // plane 2D case
            stiffInt = new BDBInt<CurlOperator<FeH1,2,Double> >(curCoef, 1.0);
          }
        } else {
          // 3D case
          stiffInt = new BDBInt<CurlOperator<FeH1,3,Double> >(curCoef, 1.0);
        }
        stiffInt->SetName("CurlCurlIntegrator");    
        stiffInt->SetFeSpace( mySpace);
        BiLinFormContext * stiffIntDescr =
            new BiLinFormContext(stiffInt, STIFFNESS );
        stiffIntDescr->SetEntities( actSDList, actSDList );
        stiffIntDescr->SetFeFunctions( myFct, myFct );
        
        assemble_->AddBiLinearForm( stiffIntDescr );
        
        // Important: Add bdb-integrator to global list, as we need them later
        // for calculation of postprocessing results
        bdbInts_[actRegion] = stiffInt;
        
        
        // ====================================================================
        //  3D CASE (additional stiffness integrator)
        // ====================================================================
        if( dim_ == 3 ) {
          BaseBDBInt * divInt = 
              new BDBInt<DivOperator<FeH1,3,Double> > (curCoef, 1.0);
          divInt->SetFeSpace( mySpace );
          divInt->SetName("DivDivIntegrator");
          BiLinFormContext * divIntDescr =
                      new BiLinFormContext(divInt, STIFFNESS );
          divIntDescr->SetEntities( actSDList, actSDList );
          divIntDescr->SetFeFunctions( myFct, myFct );
          assemble_->AddBiLinearForm( divIntDescr );
        }
        
        
        // ====================================================================
        //  3D MIXED CASE (coupling between vector and scalar potential)
        // ====================================================================
        if( isMixed_ ) {
          WARN("Mixed integrators have to be defined ....");
        }
        
      }
      
//      // ====================================================================
//      //  Standard Mass Integrator
//      // ====================================================================
//      Double density = 0.0;
//      actSDMat->GetScalar(density,DENSITY,Global::REAL);
//      shared_ptr<CoefFunction> densCoeff 
//          = CoefFunction::Generate(Global::REAL, lexical_cast<std::string>(density));
//      
//      BiLinearForm *massInt = NULL;
//      if( dim_ == 2 ) {
//        massInt = new BBInt<IdentityOperator<FeH1,2,2> >(densCoeff, 1.0);
//      } else {
//        massInt = new BBInt<IdentityOperator<FeH1,3,3> >(densCoeff, 1.0);
//      }
//      massInt->SetFeSpace( mySpace );
//
//      BiLinFormContext *massContext =  new BiLinFormContext( massInt, MASS );
//      massContext->SetEntities( actSDList, actSDList );
//      massContext->SetFeFunctions( myFct, myFct );
//      assemble_->AddBiLinearForm( massContext );

      
    }
  }
  
  
  void MagneticPDE::DefineRhsLoadIntegrators() {
    
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> mySpace = feFct->GetFeSpace();    
    
    // Loop over all coils
    for ( UInt coil = 0; coil < coilDef_.GetSize(); coil++ ) {

         // Set current region and material
         RegionIdType actRegion = coilRegionId_[coil];

         // Get current region name
         std::string regionName = ptgrid_->GetRegion().ToString(actRegion);

         // create new entity list
         shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
         actSDList->SetRegion( actRegion );

         LinearForm * curInt = NULL;
         std::string factor = coilDef_[coil]->value_ + "/" +
             lexical_cast<std::string>(coilDef_[coil]->windingCrossSection_);
         shared_ptr<CoefFunction> coef;
         if( dim_ == 3 ) {
           StdVector<std::string> currDensity(3);
           currDensity[0] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[0]);
           currDensity[1] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[1]);
           currDensity[2] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[2]);
           coef = CoefFunction::Generate(Global::REAL, currDensity);
           coef->SetCoordinateSystem(coilDef_[coil]->flowCoordSys_);
           curInt = new BUIntegrator<IdentityOperator<FeH1,3,3> >(1.0, coef);
         } else {
           StdVector<std::string> currDensity(1);
           currDensity[0] = factor;
           coef = CoefFunction::Generate(Global::REAL, currDensity);
           coef->SetCoordinateSystem(coilDef_[coil]->flowCoordSys_);
           curInt = new BUIntegrator<IdentityOperator<FeH1,2,1> >(1.0, coef);
         }
         
        
         LinearFormContext * coilContext =
             new LinearFormContext( curInt );
         coilContext->SetEntities( actSDList );
         coilContext->SetFeFunction( feFct );
         assemble_->AddLinearForm( coilContext );

    }
  }
  
  void MagneticPDE::ReadSpecialBCs() {


    // --------------------------------------------------------------------
    //   Get information about coils and open files for measurement coils
    // --------------------------------------------------------------------
    ReadCoils();

    // -----------------------------
    // Check for permanent magnets
    // -----------------------------
    ReadMagnets();

  }


 
  void MagneticPDE::DefineSolveStep()
  {
		  solveStep_ = new StdSolveStep(*this);
  }



  // ======================================================
  // COUPLING SECTION
  // ======================================================


  void MagneticPDE::InitCoupling(PDECoupling * Coupling)
  {

//    isIterCoupled_ = true;
//    ptCoupling_   = Coupling;
//
//    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
//      {
//        if (ptCoupling_->GetOutputQuantity(i) == MECH_DISPLACEMENT)
//          {
//            // Intialize the memory of the coupling values
//            ptCoupling_->CreateCouplingVector(i,isComplex_);
//          }
//
//        if (ptCoupling_->GetOutputQuantity(i) == MECH_VELOCITY)
//          {
//            // Intialize the memory of the coupling values
//            ptCoupling_->CreateCouplingVector(i,isComplex_);
//          }
//
//        if (ptCoupling_->GetOutputQuantity(i) == MECH_FORCE)
//          {
//            // Intialize the memory of the coupling values
//            ptCoupling_->CreateCouplingVector(i,isComplex_);
//
//            //now since we need a incremental formulation, initialize some necessary vectors
//            isIncrFormulation_ = true;
//          }
//      }

  }


  void MagneticPDE::CalcOutputCoupling()
  {
    REFACTOR;
  }

  bool MagneticPDE::HasOutput(SolutionType output)
  {
//
//    if (output == MECH_DISPLACEMENT || output == MECH_VELOCITY || output == MECH_FORCE)
//      return true;
//
    return false;
  }



  void MagneticPDE::CalcResults( shared_ptr<BaseResult> res ) {

    switch (res->GetResultInfo()->resultType ) {

      case MAG_POTENTIAL:
        feFunctions_[MAG_POTENTIAL]->ExtractResult( res );
        break;

      case MAG_RHS_LOAD:
        rhsFeFunctions_[MAG_POTENTIAL]->ExtractResult( res );
        break;

      case MAG_FLUX_DENSITY:
        resultFunctors_[MAG_FLUX_DENSITY]->EvalResult( res );
        break;

      default:
        WARN( "Result '" << 
              SolutionTypeEnum.ToString(res->GetResultInfo()->resultType)
              << "' type not computable by magnetic PDE" );
    }
  }

  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================
  void MagneticPDE::InitTimeStepping()
  {
//    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(TimeSchemeGLM::NEWMARK, 2) );
//    feFunctions_[MECH_DISPLACEMENT]->SetTimeScheme(myScheme);
  }


  // ******************************************************
  //   Query parameter object for information about coils
  // ******************************************************
  void MagneticPDE::ReadCoils() {

    // Check if the element "coils" is present at all.
    // Otherwise leave
    PtrParamNode coilNode = myParam_->Get( "coils", ParamNode::PASS );
    if ( !coilNode )
      return;

    // Get single coil nodes
    ParamNodeList coilNodes = coilNode->GetChildren();

    // Trigger reading in of definitions
    if( coilNodes.GetSize() > 0 ) {
      WARN("Adapt printing of coils to InfoNode");
      //Info->PrintF( pdename_, "Using the following coils:\n" );
      for( UInt i = 0; i < coilNodes.GetSize(); i++ ) {

        // get region name of actual coil
        std::string regionName = coilNodes[i]->Get("name")->As<std::string>();
        RegionIdType regionId = ptgrid_->GetRegion().Parse( regionName );

        coilRegionId_.Push_back( regionId );
        coilDef_.Push_back( shared_ptr<Coil>( new Coil( regionId,
                                                        coilNodes[i], ptgrid_) ) );
        //Info->PrintCoil( *coilDef_.Last(), analysistype_ );
      }
    }
  }


  // ********************************************************
  //   Query parameter object for information about magnets
  // ********************************************************
  void MagneticPDE::ReadMagnets() {

    // Check if the element "magnets" is present at all.
    // Otherwise leave
    PtrParamNode magnetNode = myParam_->Get( "magnets", ParamNode::PASS );
    if ( !magnetNode )
      return;

    // Get single magnet nodes
    ParamNodeList magnetNodes = magnetNode->GetChildren();

    // trigger definition of magnets
    if( magnetNodes.GetSize() > 0 ) {
      WARN("Adjust printing of permanent magnet definition to InfoNode");
//      Info->PrintF( pdename_,
//              "Found permanent magnets in the following regions:\n" );

      Double tmpDir;
      for( UInt i = 0; i < magnetNodes.GetSize(); i++ ) {

        // get region name of actual magnet
        std::string regionName = magnetNodes[i]->Get("name")->As<std::string>();
        RegionIdType regionId = ptgrid_->GetRegion().Parse( regionName );

        magnetsDomain_.Push_back( regionId );

        // read orientation
        magnetNodes[i]->GetValue( "orientX", tmpDir );
        magnetsOriX_.Push_back( tmpDir );

        magnetNodes[i]->GetValue( "orientY", tmpDir );
        magnetsOriY_.Push_back( tmpDir );

        magnetNodes[i]->GetValue( "orientZ", tmpDir );
        magnetsOriZ_.Push_back( tmpDir );

        // report name to logfile
        //Info->PrintF( pdename_, " %s\n", regionName.c_str());
      }
    }
  }


 
  void MagneticPDE::DefinePrimaryResults()
  {
    StdVector<std::string> vecComponents;
    if( dim_ == 3 ) {
      vecComponents = "x", "y", "z";
    }
    else if( isaxi_ ) {
      vecComponents = "phi";
    } 
    else {
      vecComponents = "z";
    }
    
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];
    
    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> res1(new ResultInfo);
    res1->resultType = MAG_POTENTIAL;
    res1->dofNames = vecComponents;
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::VECTOR;
    res1->unit = "Vs/m";
    res1->SetFeFunction(feFct);
    results_.Push_back( res1 );
    availResults_.insert( res1 );
    // define related interpolatory method
    shared_ptr<BaseFieldFunctor> vFunc;
    if( isComplex_ ) {
      if( dim_ == 2 ) {
      vFunc.reset(
          new FieldInterpolFunctor<IdentityOperator<FeH1,2,2,Complex>,
          Complex>(feFct, res1));
      } else {
        vFunc.reset(
            new FieldInterpolFunctor<IdentityOperator<FeH1,3,3,Complex>,
            Complex>(feFct, res1));
      }
    } else {
      if( dim_ == 2 ) {
      vFunc.reset(
          new FieldInterpolFunctor<IdentityOperator<FeH1,2,2>,
          Double>(feFct, res1));
      } else {
        vFunc.reset(
            new FieldInterpolFunctor<IdentityOperator<FeH1,3,3>,
            Double>(feFct, res1));
      }
    }
    resultFunctors_[MAG_POTENTIAL] = vFunc;
    fieldFunctors_[MAG_POTENTIAL] = vFunc;
    feFct->SetResultInfo(res1);
    //
//    if( subType_ == "3d" ) {
//      dispDofNames = "x", "y", "z";
//      stressDofNames = "xx", "yy", "zz", "yz", "xz", "xy";
//
//    } else if( subType_ == "planeStrain" ) {
//      dispDofNames = "x", "y";
//      stressDofNames = "xx", "yy", "xy";
//
//    } else if( subType_ == "planeStress" ) {
//      dispDofNames = "x", "y";
//      stressDofNames = "xx", "yy", "xy";
//
//    } else if( subType_ == "axi" ) {
//      dispDofNames = "r", "z";
//      stressDofNames = "rr", "zz", "rz", "phiphi";
//    }
//
//    // === MECHANIC DISPLACEMENT ===
//    shared_ptr<ResultInfo> disp(new ResultInfo);
//    disp->resultType = MECH_DISPLACEMENT;
//    disp->dofNames = dispDofNames;
//    disp->unit = "m";
//    disp->entryType = ResultInfo::VECTOR;
//    disp->SetFeFunction(feFunctions_[MECH_DISPLACEMENT]);
//    disp->definedOn = ResultInfo::NODE;
//    feFunctions_[MECH_DISPLACEMENT]->SetResultInfo(disp);
//    
//    // this defines the primary unknown
//    results_.Push_back( disp );
//    availResults_.insert( disp);
//    
//    // define functor for interpolation of the field
//    shared_ptr<BaseFeFunction> feFct = feFunctions_[MECH_DISPLACEMENT];
//    shared_ptr<BaseFieldFunctor> dFunc;
//       if( isComplex_ ) {
//         dFunc.reset(
//             new FieldInterpolFunctor<IdentityOperator<FeH1,2,2,Complex>,
//             Complex>(feFct, disp));
//       } else {
//         dFunc.reset(
//             new FieldInterpolFunctor<IdentityOperator<FeH1,2,2,Double>,
//             Double>(feFct, disp));
//       }
//       resultFunctors_[MECH_DISPLACEMENT] = dFunc;
//       fieldFunctors_[MECH_DISPLACEMENT] = dFunc;
//
//    if ( (analysistype_ == TRANSIENT) || (analysistype_ == HARMONIC))
//    {
//      // === MECHANIC VELOCITY ===
//      shared_ptr<ResultInfo> vel(new ResultInfo);
//      vel->resultType = MECH_VELOCITY;
//      vel->dofNames = dispDofNames;
//      vel->unit = "m/s";
//      vel->entryType = disp->entryType;
//      vel->definedOn = disp->definedOn;
//      availResults_.insert( vel );
//
//      // === MECHANIC ACCELERATION ===
//      shared_ptr<ResultInfo> acc(new ResultInfo);
//      acc->resultType = MECH_ACCELERATION;
//      acc->dofNames = dispDofNames;
//      acc->unit = "m/s^2";
//      acc->entryType = disp->entryType;
//      acc->definedOn = disp->definedOn;
//      availResults_.insert( acc );
//    }
//
//    // === MECHANIC RHS ===
//    shared_ptr<ResultInfo> rhs(new ResultInfo);
//    rhs->resultType = MECH_RHS_LOAD;
//    rhs->dofNames = dispDofNames;
//    rhs->unit = "N";
//    rhs->entryType = disp->entryType;
//    rhs->definedOn = disp->definedOn;
//    availResults_.insert( rhs );
//    postProcResults_[MECH_RHS_LOAD] = MECH_DISPLACEMENT;
//    //
    
  }
    
  void MagneticPDE::DefinePostProcResults() {
    StdVector<std::string> vecComponents;
       if( dim_ == 3 ) {
         vecComponents = "x", "y", "z";
       }
       else if( isaxi_ ) {
         vecComponents = "r", "phi";
       } 
       else {
         vecComponents = "x", "y";
       }
       shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];
             
       // === MAGNETIC FLUX DENSITY ===
       shared_ptr<ResultInfo> flux(new ResultInfo);
       flux->resultType = MAG_FLUX_DENSITY;
       flux->dofNames = vecComponents;
       flux->unit = "Vs/m^2";
       flux->definedOn = ResultInfo::ELEMENT;
       flux->entryType = ResultInfo::VECTOR;
       availResults_.insert( flux );
       postProcResults_[MAG_FLUX_DENSITY] = MAG_POTENTIAL;
       shared_ptr<BaseFieldFunctor> bFunc;
       if( isComplex_ ) {
         bFunc.reset(new DiffFieldFunctor<Complex>(feFct, flux));
       } else {
         bFunc.reset(new DiffFieldFunctor<Double>(feFct, flux));
       }
       resultFunctors_[MAG_FLUX_DENSITY] = bFunc;
       fieldFunctors_[MAG_FLUX_DENSITY] = bFunc;
//    
//    
//    StdVector<std::string > dispDofNames;
//    dispDofNames = feFunctions_[MECH_DISPLACEMENT]->GetResultInfo()->dofNames;
//    shared_ptr<BaseFeFunction> feFct = feFunctions_[MECH_DISPLACEMENT];
//    
//    // === MECHANIC STRESS ===
//    shared_ptr<ResultInfo> stress(new ResultInfo);
//    stress->resultType = MECH_STRESS;
//    stress->dofNames = stressComponents;
//    stress->unit =  "N/m^2";
//    stress->entryType = ResultInfo::TENSOR;
//    stress->definedOn = ResultInfo::ELEMENT;
//    availResults_.insert( stress );
//    postProcResults_[MECH_STRESS] = MECH_DISPLACEMENT;
//    shared_ptr<BaseFieldFunctor> sigmaFunc;
//    if( isComplex_ ) {
//      sigmaFunc.reset(new FluxFieldFunctor<Complex>(feFct, stress));
//    } else {
//      sigmaFunc.reset(new FluxFieldFunctor<Double>(feFct, stress));
//    }
//    resultFunctors_[MECH_STRESS] = sigmaFunc;
//    fieldFunctors_[MECH_STRESS] = sigmaFunc;
//       
//    
//    
//    // === MECHANIC STRAIN ===
//    shared_ptr<ResultInfo> strain(new ResultInfo);
//    strain->resultType = MECH_STRAIN;
//    strain->dofNames = stressComponents;
//    strain->unit =  "";
//    strain->entryType = ResultInfo::TENSOR;
//    strain->definedOn = ResultInfo::ELEMENT;
//    availResults_.insert( strain );
//    postProcResults_[MECH_STRAIN] = MECH_DISPLACEMENT;
//    shared_ptr<BaseFieldFunctor> strainFunc;
//    if( isComplex_ ) {
//      strainFunc.reset(new DiffFieldFunctor<Complex>(feFct, strain));
//    } else {
//      strainFunc.reset(new DiffFieldFunctor<Double>(feFct, strain));
//    }
//    resultFunctors_[MECH_STRAIN] = strainFunc;
//    fieldFunctors_[MECH_STRAIN] = strainFunc;
//
//    // ============================
//    // Initialize result functors:
//    // ============================
//    // 1) Loop over all BDB-integrators
//    std::map<RegionIdType, BaseBDBInt*>::iterator it = bdbInts_.begin();
//    for( ; it != bdbInts_.end(); ++it ) {
//      RegionIdType region = it->first;
//      BaseBDBInt* bdb = it->second;
//
//      // 2) pass integrators to functors
//      sigmaFunc->AddIntegrator(bdb, region);
//      strainFunc->AddIntegrator(bdb, region);
//    }
//    
       
       // ============================
       // Initialize result functors:
       // ============================
       // 1) Loop over all BDB-integrators
       std::map<RegionIdType, BaseBDBInt*>::iterator it = bdbInts_.begin();
       for( ; it != bdbInts_.end(); ++it ) {
         RegionIdType region = it->first;
         BaseBDBInt* bdb = it->second;

         // 2) pass integrators to functors
         bFunc->AddIntegrator(bdb, region);
       }

  }
  
  
  std::map<SolutionType, shared_ptr<FeSpace> >
   MagneticPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
     
    // and 3D analysis
    isMixed_ = false;
    if( analysistype_ != BasePDE::STATIC && dim_ == 3 ) {
      isMixed_ = true;
    }
    
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    
    if( formulation == "default" || formulation == "H1" ){
      
      // 1) create space for magnetic vector potential
      PtrParamNode potSpaceNode = infoNode->Get("magPotential");
      crSpaces[MAG_POTENTIAL] =
          FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1);
      crSpaces[MAG_POTENTIAL]->Init();
      
      // 1) create space for electric scalar potential
      if( isMixed_ ) {
        crSpaces[ELEC_POTENTIAL] =
            FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1);
        crSpaces[ELEC_POTENTIAL]->Init();
      }

    }else{
       EXCEPTION( "The formulation " << formulation 
                  << "of the mechanic PDE is not known!" );
     }
     return crSpaces;
   }

} // end namespace CoupledField
