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
#include "Forms/BiLinForms/ADBInt.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/DivOperator.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"


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
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // create new entity list and add it fefunction
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );
      myFct->AddEntityList( actSDList );
      
      // --------------------------
      //  Set region approximation
      // --------------------------
      // --- Set the approximation for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")
          ->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);
      
      
      
      // ====================================================================
      //  Standard Linear CASE (2D AND 3D)
      // ====================================================================
      if( !nonLin_ ) {
        BaseBDBInt * stiffInt = NULL;
        PtrCoefFct curCoef = 
            actMat->GetTensorCoefFnc( MAG_RELUCTIVITY, tensorType,
                                     Global::REAL );
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
        //  MASS - MATRIX (all dimensions)
        // ====================================================================
        // Note: here we have currently the problem, that we can only treat
        // constant / double valued conductivity values, as we must check for 
        // non-zero. However, how do we do this in case of non-constant parameters?
        Double conductivity;
        
        materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
        PtrCoefFct conducCoef =
            materials_[actRegion]->GetScalCoefFnc(MAG_CONDUCTIVITY,Global::REAL);
//                                         lexical_cast<std::string>(conductivity));
        {
          BiLinearForm *massInt = NULL;
          if( dim_ == 2 ) {
            massInt = new BBInt<IdentityOperator<FeH1,2,1> >(conducCoef,1.0);
          } else {
            massInt = new BBInt<IdentityOperator<FeH1,3,3> >(conducCoef,1.0);
          }
          massInt->SetName("MassIntegrator");
          BiLinFormContext * massContext = new BiLinFormContext(massInt, MASS );
          massContext->SetEntities( actSDList, actSDList );
          massContext->SetFeFunctions( myFct, myFct );
          assemble_->AddBiLinearForm( massContext );
        }
        
        // ====================================================================
        //  3D MIXED CASE (coupling between vector and scalar potential)
        // ====================================================================
        if( isMixed_ && conductivity > EPS ) {
          shared_ptr<BaseFeFunction> potFct = feFunctions_[ELEC_POTENTIAL];
          shared_ptr<FeSpace> potSpace = potFct->GetFeSpace();
          // Important: set region approximation also at space for scalar 
          // potential
          potSpace->SetRegionApproximation(actRegion, polyId,integId);
          
          // add region to feFunction for scalar potential
          potFct->AddEntityList( actSDList );
          
          // ------------------------------
          // diagonal mass matrix M_phiphi
          // ------------------------------
          {
            BiLinearForm * phiDivInt = 
                new BBInt<GradientOperator<FeH1,3,Double> > (conducCoef, 1.0);
            phiDivInt->SetName("MassIntegrator_PhiPhi");
            BiLinFormContext * massContext = 
                new BiLinFormContext(phiDivInt, MASS );
            massContext->SetEntities( actSDList, actSDList );
            massContext->SetFeFunctions( potFct, potFct );
            assemble_->AddBiLinearForm( massContext );
          }
          
          // ------------------------------
          // off-diagonal mass matrix M_A_phi
          // ------------------------------
          {
            BiLinearForm * cplInt = 
                new ABInt< IdentityOperator<FeH1,3,3,Double>,
                           GradientOperator<FeH1,3,Double> > (conducCoef, 1.0);
            cplInt->SetName("MassIntegrator_Coupling_Phi_A");
            BiLinFormContext * cplContext = 
                new BiLinFormContext(cplInt, MASS );
            cplContext->SetCounterPart(true);
            cplContext->SetEntities( actSDList, actSDList );
            cplContext->SetFeFunctions( myFct, potFct );
            assemble_->AddBiLinearForm( cplContext );
          }
        } // mixed
        
      } // linear part
    } // regions
  }
  
  
  void MagneticPDE::DefineRhsLoadIntegrators() {
    
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> mySpace = feFct->GetFeSpace();    
    
    // Loop over all coils
    for ( UInt coil = 0; coil < coilDef_.GetSize(); coil++ ) {

         // Set current region and material
         RegionIdType actRegion = coilRegionId_[coil];

         // Get current region name
         std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

         // create new entity list
         shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
         actSDList->SetRegion( actRegion );

         LinearForm * curInt = NULL;
         std::string factor = coilDef_[coil]->value_ + "/" +
             lexical_cast<std::string>(coilDef_[coil]->windingCrossSection_);
         PtrCoefFct coef;
         // ===========
         //  3D CASE
         // ===========
         if( dim_ == 3 ) {
           StdVector<std::string> currDensity(3);
           currDensity[0] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[0]);
           currDensity[1] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[1]);
           currDensity[2] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[2]);

           if( isComplex_ ) {
             StdVector<std::string> phaseVec(3);
             phaseVec.Init(coilDef_[coil]->phase_);
             coef = CoefFunction::Generate(Global::COMPLEX, currDensity, phaseVec );
             coef->SetCoordinateSystem(coilDef_[coil]->flowCoordSys_);
             curInt = new BUIntegrator<IdentityOperator<FeH1,3,3>, Complex >(1.0, coef);
           } else {
             coef = CoefFunction::Generate(Global::REAL, currDensity);
             coef->SetCoordinateSystem(coilDef_[coil]->flowCoordSys_);
             curInt = new BUIntegrator<IdentityOperator<FeH1,3,3>, Double >(1.0, coef);
           } // complex
         } else {
           // ===============
           //  2D / AXI CASE
           // ===============
           StdVector<std::string> currDensity(1);
           currDensity[0] = factor;
           
           if( isComplex_ ) {
             StdVector<std::string> phaseVec(1);
             phaseVec.Init(coilDef_[coil]->phase_);
             curInt = new BUIntegrator<IdentityOperator<FeH1,2,1>, Complex >(1.0, coef);
             coef = CoefFunction::Generate(Global::REAL, currDensity);
             coef->SetCoordinateSystem(coilDef_[coil]->flowCoordSys_);
           } else {
             curInt = new BUIntegrator<IdentityOperator<FeH1,2,1>, Double >(1.0, coef);
             coef = CoefFunction::Generate(Global::REAL, currDensity);
             coef->SetCoordinateSystem(coilDef_[coil]->flowCoordSys_);
           } // complex
         } // dimension
         
        
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


  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================
  void MagneticPDE::InitTimeStepping() {
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(TimeSchemeGLM::TRAPEZOIDAL, 0) );
    feFunctions_[MAG_POTENTIAL]->SetTimeScheme(myScheme);

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
        RegionIdType regionId = ptGrid_->GetRegion().Parse( regionName );

        coilRegionId_.Push_back( regionId );
        coilDef_.Push_back( shared_ptr<Coil>( new Coil( regionId,
                                                        coilNodes[i], ptGrid_) ) );
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
        RegionIdType regionId = ptGrid_->GetRegion().Parse( regionName );

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
    
    shared_ptr<BaseFeFunction> vecFct = feFunctions_[MAG_POTENTIAL];
    
    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> res1(new ResultInfo);
    res1->resultType = MAG_POTENTIAL;
    res1->dofNames = vecComponents;
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::VECTOR;
    res1->unit = "Vs/m";
    res1->SetFeFunction(vecFct);
    results_.Push_back( res1 );
    availResults_.insert( res1 );
    vecFct->SetResultInfo(res1);
    DefineFieldResult( vecFct, res1);
    
    
    // === MAGNETIC SCALAR POTENTIAL ===
    if (isMixed_) {
      shared_ptr<BaseFeFunction> scalFct = feFunctions_[ELEC_POTENTIAL];
      shared_ptr<ResultInfo> res2(new ResultInfo);
      res2->resultType = ELEC_POTENTIAL;
      res2->dofNames = "";
      res2->unit = "V";
      res2->definedOn = ResultInfo::NODE;
      res2->entryType = ResultInfo::SCALAR;
      results_.Push_back( res2 );
      availResults_.insert( res2 );
      scalFct->SetResultInfo(res2);
      DefineFieldResult( scalFct, res2 );
    }
    
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
    shared_ptr<CoefFunctionFormBased> bFunc;
    if( isComplex_ ) {
      bFunc.reset(new CoefFunctionBOp<Complex>(feFct, flux));
    } else {
      bFunc.reset(new CoefFunctionBOp<Double>(feFct, flux));
    }
    DefineFieldResult( bFunc, flux );


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
          FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[MAG_POTENTIAL]->Init(solStrat_);
      
      // 1) create space for electric scalar potential
      if( isMixed_ ) {
        crSpaces[ELEC_POTENTIAL] =
            FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
        crSpaces[ELEC_POTENTIAL]->Init(solStrat_);
      }

    }else{
       EXCEPTION( "The formulation " << formulation 
                  << "of the mechanic PDE is not known!" );
     }
     return crSpaces;
   }

} // end namespace CoupledField
