#include <fstream>

#include "MagEdgePDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"

#include "Driver/TimeSchemes/Trapezoidal.hh"
#include "CoupledPDE/PDECoupling.hh"
#include "Driver/Assemble.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "FeBasis/HCurl/FeSpaceHCurlHi.hh"
#include "FeBasis/HCurl/HCurlElems.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "Domain/CoefFunction/CoefFunctionExpression.hh"

#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
namespace CoupledField {

// declare class specific logging stream
DECLARE_LOG(magEdgePde)
DEFINE_LOG(magEdgePde, "magEdgePde")


  // **************
  //  Constructor
  // **************
  MagEdgePDE::MagEdgePDE( Grid * aptgrid, PtrParamNode paramNode )
    :SinglePDE( aptgrid, paramNode ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticEdge";
    pdematerialclass_ = ELECTROMAGNETIC;
    maxTimeDerivOrder_ = 1;

    // check if we have a 3d setup
    bool is3d = param->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
    if ( !is3d )
      EXCEPTION("MagEdgePDE is just implemented for 3D setups!");
  }


  // *************
  //  Destructor
  // *************
  MagEdgePDE::~MagEdgePDE() {
  }



  void MagEdgePDE::ReadSpecialBCs() {


    // --------------------------------------------------------------------
    //   Get information about coils and open files for measurement coils
    // --------------------------------------------------------------------
    ReadCoils();

    // -----------------------------
    // Check for permanent magnets
    // -----------------------------
    ReadMagnets();

  }


  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void MagEdgePDE::InitNonLin() {

    SinglePDE::InitNonLin();
  }


  // *****************************
  //   Definition of Integrators
  // *****************************
  void MagEdgePDE::DefineIntegrators() {

    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // flag for updatedLagrange formulation
    //bool upLagrangeForm = true;

    
    // initially, check for regularization factor
    Double regularizationFactor = 1e-6;
    myParam_->GetValue("penaltyFactor", regularizationFactor, ParamNode::PASS);
    
    //==============================================================
    //begin new implementation
    //==============================================================
    shared_ptr<FeSpace> mySpace = feFunctions_[MAG_POTENTIAL]->GetFeSpace();
    for(UInt iRegion = 0; iRegion < subdoms_.GetSize() ; iRegion ++){
      actRegion = subdoms_[iRegion];
      actMat    = materials_[actRegion];
      std::string regionName = ptgrid_->GetRegion().ToString(actRegion);

      shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_POTENTIAL];
      shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();

      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);
      
      //get possible nonlinearities defined in this region
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion]; 

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );
      // Switch, if region is linear / nonlinear
      if ( nonLinTypes.GetSize() > 0 ) {
        REFACTOR;
       // =================================
       //  Nonlinear Stiffness Integrator
       //// =================================
       //nLinCurlCurlEdgeInt* curlcurlNL =
       //    new nLinCurlCurlEdgeInt( actMat, upLagrangeForm );
       //curlcurlNL->SetNonLinMethod( nonLinMethod_ );
       //curlcurlNL->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
       //curlcurlNL->SetFeSpace( feSpace );
       ////BiLinFormContext * stiffContext = new BiLinFormContext( curlcurlNL, STIFFNESS );
       ////stiffContext->SetEntities( actSDList, actSDList );
       ////stiffContext->SetFeFunctions( feFunc, feFunc );
       ////assemble_->AddBiLinearForm( stiffContext );
       //
       //// we have to save the bilinear form for later usage
      //// nlinBilinForms_[actRegion] = curlcurlNL;
       //
       //// =================================
       ////  Nonlinear RHS-integrator
       //// =================================
       //nLinMagEdge_linFormInt* rhsSource
       //= new nLinMagEdge_linFormInt( actMat, upLagrangeForm);
       //rhsSource->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
       //rhsSource->SetFeSpace( feSpace );
       //LinearFormContext * rhsContext =
       //    new LinearFormContext( rhsSource );
       //rhsContext->SetEntities( actSDList );
       //rhsContext->SetFeFunction( feFunc );
       //assemble_->AddLinearForm( rhsContext );

      } else {

       // ***************************************
       // LINEAR PART
       // ***************************************

       // ===============================
       //  Standard Stiffness Integrator
       // ===============================
        shared_ptr<CoefFunction> curCoef = actMat->GetCoefFunction(MAG_RELUCTIVITY,FULL,Global::REAL);
        BDBInt< CurlOperator ,FeHCurl >* curlcurl;
        curlcurl = new BDBInt< CurlOperator ,FeHCurl >(curCoef,1.0) ;

       BiLinFormContext * stiffContext =
           new BiLinFormContext(curlcurl, STIFFNESS );
       stiffContext->SetEntities( actSDList, actSDList );
       stiffContext->SetFeFunctions( feFunc, feFunc );
       assemble_->AddBiLinearForm( stiffContext );
       // Important: Add bdb-integrator to global list, as we need them later
       // for calculation of postprocessing results
       bdbInts_[actRegion] = curlcurl;

       // === Additional RHS integrator in case of Non-linearity ===
       if ( nonLin_ == true ) {
         REFACTOR;
         //nLinMagEdge_linFormInt* rhsSource =
         //    new nLinMagEdge_linFormInt( actMat, upLagrangeForm );
         //rhsSource->SetFeSpace( feSpace );
         //rhsSource->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
         //
         //LinearFormContext * rhsContext =
         //    new LinearFormContext( rhsSource );
         //rhsContext->SetEntities( actSDList );
         //rhsContext->SetFeFunction( feFunc );
         //assemble_->AddLinearForm( rhsContext );
       }
     } // END OF NOLIN / LIN PART

      // ============================
      // Standard Mass Matrix
      // ============================
      Double conductivity = 0.0; // , maxPerm = 0.0; // TODO: Check if this is still needed
      bool scaleByEdgeSize = false;
      materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);

      if ( conductivity < 1e-10 || analysistype_ == STATIC ) {
        Matrix<Double> reluc; 
        // get tensor of permeability and determine max. value
        materials_[actRegion]->GetTensor( reluc, MAG_RELUCTIVITY, Global::REAL );
        scaleByEdgeSize = true;
        //regularizationFactor = 1e-6;
        conductivity =  regularizationFactor * reluc[0][0];
      }

      shared_ptr<CoefFunction> coeff =
          CoefFunction::Generate(Global::REAL, lexical_cast<std::string>(conductivity));
      BBIntMassEdge<ScaledByEdgeIdentityOperator,FeHCurl,Double> *massInt;
      massInt = new BBIntMassEdge<ScaledByEdgeIdentityOperator,FeHCurl,Double>(coeff,1.0);

      BiLinFormContext * massContext;
      if ( analysistype_ == STATIC) {
        // we have to guarantee, that we add some mass to
        // the curl-curl-matrix!!
        massContext =  new BiLinFormContext(massInt, STIFFNESS );
      } else {
        massContext =
            new BiLinFormContext(massInt, MASS );
      }
      massContext->SetEntities( actSDList, actSDList );
      massContext->SetFeFunctions( feFunc, feFunc );
      assemble_->AddBiLinearForm( massContext );

      feFunc->AddEntityList( actSDList );
      // ============================
      // COIL INTEGRATORS
      // ============================
      // If this subdomain is a coil we have to do special things
      for ( UInt coil = 0; coil < coilDef_.GetSize(); coil++ ) {
        if ( actRegion == coilRegionId_[coil] ) {
          BUIntegrator<IdentityOperator,FeHCurl,Double>* curInt;
          
          std::string factor = coilDef_[coil]->value_ + "/" +
            lexical_cast<std::string>(coilDef_[coil]->windingCrossSection_);

          StdVector<std::string> currDensity(3);
          currDensity[0] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[0]);
          currDensity[1] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[1]);
          currDensity[2] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[2]);
          //shared_ptr<CoefFunctionExpression<Double> > coef (new CoefFunctionExpression<Double>());
          //coef->SetVector( currDensity );
          shared_ptr<CoefFunction> coef(CoefFunction::Generate(Global::REAL, currDensity));
          coef->SetCoordinateSystem(coilDef_[coil]->flowCoordSys_);
          
          curInt = new BUIntegrator<IdentityOperator,FeHCurl,Double>(1.0, coef);
          LinearFormContext * coilContext =
              new LinearFormContext( curInt );
          coilContext->SetEntities( actSDList );
          coilContext->SetFeFunction( feFunc );
          assemble_->AddLinearForm( coilContext );
        }
      }
        
    }


        // ============================
        // PERMANENT MAGNETS
        // ============================
        //
//      // check, if this subdomain is a permanent magnet
//      for ( UInt perm = 0; perm < magnetsDomain_.GetSize(); perm++ ) {
//        if ( actRegion == magnetsDomain_[perm] ) {
//          EXCEPTION("Currently magnetic 3D with edge elements do not support permanent magnets");
//
//          Vector<Double> magnetization(dim_);
//          magnetization[0] = magnetsOriX_[perm];
//          magnetization[1] = magnetsOriY_[perm];
//          magnetization[2] = magnetsOriZ_[perm];
//
//          // Get reluctivity for this domain and perform consistency check
//          Double reluctivity;
//          actMat->GetScalar(reluctivity,MAG_RELUCTIVITY,Global::REAL);
//
//          std::string fncname = "none";
//          LinearForm *permSource =
//            new MagPerm3DInt(magnetization, reluctivity,
//                             isaxi_, upLagrangeForm );
//
//          LinearFormContext * permContext =
//            new LinearFormContext( permSource );
//          permContext->SetPtPde( this );
//          permContext->SetResult( results_[0], actSDList );
//          assemble_->AddLinearForm( permContext );
//        }
//      }
//    }
  }

  void MagEdgePDE::DefineSolveStep()
  {
    SolveStepMagEdge *magSolveStep = new SolveStepMagEdge(*this); 
    
    solveStep_ = magSolveStep; 
  }


  // ======================================================
  // TIME-STEPPING SECTION
  // ======================================================

  void MagEdgePDE::InitTimeStepping() {
    TS_alg_ = new Trapezoidal( algsys_, olasNode_ );
  }

  // ***************
  //   CalcResults
  // ***************
  void MagEdgePDE::CalcResults ( shared_ptr<BaseResult> res ) {

    switch (res->GetResultInfo()->resultType ) {
      case MAG_RHS_LOAD:
        feFunctions_[MAG_POTENTIAL]->ExtractResult( res );
        break;

      case MAG_FLUX_DENSITY:
        resultFunctors_[MAG_FLUX_DENSITY]->EvalResult(res);
//        }
        break;
        
      case MAG_POTENTIAL:
        resultFunctors_[MAG_POTENTIAL]->EvalResult(res);
//        if( isComplex_ ) {
//          CalcVecPotential<Complex>( res );
//        } else {
//          CalcVecPotential<Double>( res );
//        }
        break;
      
//      case MAG_POTENTIAL_DIV:
//        if( isComplex_ ) {
//          CalcMagPotentialDiv<Complex>( res );
//        } else {
//          CalcMagPotentialDiv<Double>( res );
//        }
//        break;
//      
//      case MAG_EDDY_CURRENT:
//         if( isComplex_ ) {
//           CalcEddyCurrent<Complex>( res );
//         } else {
//           CalcEddyCurrent<Double>( res );
//         }
//         break;
//         
      case MAG_ELEM_PERMEABILITY:
        if( isComplex_ ) {
          CalcPermeability<Complex>( res );
        } else {
          CalcPermeability<Double>( res );
        }
        break;
         
      case MAG_ENERGY:
        resultFunctors_[MAG_ENERGY]->EvalResult(res);
        break;

      default:
        WARN( "Resulttype not computable by magnetic PDE" );
    }

  }
  
//  void MagEdgePDE::CalcSpecialResults() {
//    
//
//    // fetch fluxlist
//    for( UInt iList = 0; iList < calcFlux_.GetSize(); ++iList ) {
//#ifndef USE_INTERPOLATION
//    EXCEPTION("Special Results can just be calculated with INTERPOLATION active");
//#else
//      FluxAtPoints & actList = calcFlux_[iList];
//
//      // open file
//      std::ofstream  out(actList.fileName.c_str(),std::ios::out );
//
//      out << "# Flux density\n";
//      out << "# #elem \t x\t y \t z \t x(Vs/Am) \ty(Vs/Am( \t z(Vs/Am)\txi \teta \tzeta\n";
//
//      // Loop over all points 
//      Vector<Double> locCoord, flux;
//      for( UInt iPoint = 0; iPoint < actList.points.GetSize(); iPoint++) { 
//        const Elem * ptElem = NULL;
//        const Vector<Double> & globCoord =actList.points[iPoint].coord;
//        if( globCoord.GetSize() == 0) continue;
//        ptElem = ptgrid_->GetElemAtGlobalCoord( globCoord, locCoord );
//        if( !ptElem ) {
//          std::string warnStr = "Cold not find element at position " 
//              + actList.points[iPoint].coord.ToString(); 
//          WARN( warnStr.c_str());
//        } else {
//          LocPoint lp(locCoord);
//          CalcFluxDensityAtIP(ptElem, lp, flux);
//
//          // write to file
//          out << ptElem->elemNum << "\t";
//          out << globCoord[0] << "\t" << globCoord[1] << "\t" << globCoord[2] << "\t";
//          out << flux[0] << "\t" << flux[1] << "\t" << flux[2] << "\t"
//              << locCoord[0] << "\t" << locCoord[1] << "\t" << locCoord[2] << std::endl;
//        }
//
//      }
//      // close file
//      out.close();
//#endif
//    }
//
//  }



  template<class TYPE>
  void MagEdgePDE::CalcEddyCurrent( shared_ptr<BaseResult> result ) {
    EXCEPTION("Not yet adapted to new implementation")

//    // fetch result object and convert data
//    Result<TYPE> &  actSol =
//      dynamic_cast<Result<TYPE>&>(*result);
//    EntityIterator it = actSol.GetEntityList()->GetIterator();
//    Vector<TYPE> & actVal = actSol.GetVector();
//    UInt jEddyDofs = 3;
//    actVal.Resize( actSol.GetEntityList()->GetSize() * jEddyDofs );
//
//    Vector<TYPE> jEddyElem;
//    for ( it.Begin(); !it.IsEnd(); it++ ) {
//      CalcEddyCurrentAtIP( it, 0, jEddyElem );
//      for( UInt iDof = 0; iDof < jEddyDofs; iDof++ ) {
//        actVal[it.GetPos()*jEddyDofs + iDof] = jEddyElem[iDof];
//      }
//    }
  }
  
    
    //    EXCEPTION("Not yet adapted to new implementation")
//    Matrix<Double> elemmat, ptCoord, massMat;
//
//    StdVector<UInt> connecth;
//    StdVector<Integer> Eqns;  
//    Vector<TYPE> help;
//    BaseForm *bilinear_stiff = NULL; 
//
//    // get result
//    Result<TYPE> &  actSol = 
//      dynamic_cast<Result<TYPE>&>(*result);      
//      EntityIterator regionIt = actSol.GetEntityList()->GetIterator();
//
//      // resize vector
//      Vector<TYPE> & actVal = actSol.GetVector();
//      actVal.Resize( actSol.GetEntityList()->GetSize() );
//
//      // Loop over regions
//      for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {
//
//
//        // Create stiffness integrator
//        if ( regionNonLinType_[regionIt.GetRegion()] != NO_NONLINEARITY ) {
//
////          bilinear_stiff = new nLinCurlCurlEdgeInt( materials_[regionIt.GetRegion()],
////                                                    true );
////
////          bilinear_stiff->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
////
////          // VERY IMPORTANT: Set nonlinear-method "fixpoint", as otherwise also
////          // the Frechet part of the stiffness is calculated!
////          bilinear_stiff->SetNonLinMethod( "fixPoint" );
//        } else {
//          bilinear_stiff = new CurlCurlEdgeInt( materials_[regionIt.GetRegion()], true);
//
//        }
//        ElemList actSDList(ptgrid_ );
//        actSDList.SetRegion( regionIt.GetRegion() );
//        EntityIterator elemIt = actSDList.GetIterator();
//
//        TYPE energy = 0.0;
//        Vector<TYPE> magvecpot;   
//        for ( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
//          sol_->GetElemSolution(magvecpot, elemIt);
//          bilinear_stiff->CalcElementMatrix(elemmat, elemIt, elemIt);
//          help = elemmat * magvecpot;
//          TYPE temp = 0.5 * (help * magvecpot);
//          energy += temp;
//        }
//        
//        actVal[regionIt.GetPos()] = energy;
//        delete bilinear_stiff;   
//      }



  // ******************************************************
  //   Query parameter object for information about coils
  // ******************************************************
  void MagEdgePDE::ReadCoils() {

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
  void MagEdgePDE::ReadMagnets() {

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


  void MagEdgePDE::DefinePrimaryResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> res1(new ResultInfo);
//    shared_ptr<AnsatzFct> fct(new NedelecFct);
    res1->resultType = MAG_POTENTIAL;
    res1->dofNames = "";
    res1->unit = "Vs/m";
    res1->definedOn = ResultInfo::EDGE;
    res1->entryType = ResultInfo::SCALAR;
//    res1->fctType = fct;
    
    results_.Push_back( res1 );
    availResults_.insert( res1 );
    
    feFunctions_[MAG_POTENTIAL]->SetResultInfo(res1);
    //feFunctions_[MAG_POTENTIAL]->SetPDE(shared_ptr<MagEdgePDE>(this));
   
  }
  
  void MagEdgePDE::DefinePostProcResults() {

      StdVector<std::string> vecComponents;
      vecComponents = "x", "y", "z";

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

      // === MAGNETIC VECTOR POTENTIAL ===
      shared_ptr<ResultInfo> magPot(new ResultInfo);
      magPot->resultType = MAG_POTENTIAL;
      magPot->dofNames = vecComponents;
      magPot->unit = "A/m^2";
      magPot->definedOn = ResultInfo::ELEMENT;
      magPot->entryType = ResultInfo::VECTOR;
      availResults_.insert( magPot );
      postProcResults_[MAG_POTENTIAL] = MAG_POTENTIAL;
      shared_ptr<BaseFieldFunctor> aFunc;
      if( isComplex_ ) {
        aFunc.reset(
            new FieldInterpolFunctor<IdentityOperator,
            FeHCurl,
            Complex>(feFct, magPot));
      } else {
        aFunc.reset(
            new FieldInterpolFunctor<IdentityOperator,
            FeHCurl,
            Double>(feFct, magPot));
      }
      resultFunctors_[MAG_POTENTIAL] = aFunc;
      fieldFunctors_[MAG_POTENTIAL] = aFunc;

      //resultFncs_[MAG_POTENTIAL] = &MagEdgePDE::CalcFluxDensityAtIP;
      
  //    // === EDDY CURRENT DENSITY ===
  //    shared_ptr<ResultInfo> eddy(new ResultInfo);
  //    eddy->resultType = MAG_EDDY_CURRENT;
  //    eddy->dofNames = vecComponents;
  //    eddy->unit = "A/m^2";
  //    eddy->definedOn = ResultInfo::ELEMENT;
  //    eddy->entryType = ResultInfo::VECTOR;
  //    eddy->fctType = shared_ptr<ConstFct>(new ConstFct() );
  //    availResults_.insert( eddy );
  //    
  //    // === DIVERGENCE OF MAGNETIC POTENTIAL ===
  //    shared_ptr<ResultInfo> div(new ResultInfo);
  //    div->resultType = MAG_POTENTIAL_DIV;
  //    div->dofNames = "";
  //    div->unit = "Vs/m^3";
  //    div->definedOn = ResultInfo::ELEMENT;
  //    div->entryType = ResultInfo::SCALAR;
  //    div->fctType = shared_ptr<ConstFct>(new ConstFct() );
  //    availResults_.insert( div );
      
      // === PERMEABILITY  ===
      shared_ptr<ResultInfo> perm(new ResultInfo);
      perm->resultType = MAG_ELEM_PERMEABILITY;
      perm->dofNames = "";
      perm->unit = "Vs/Am";
      perm->definedOn = ResultInfo::ELEMENT;
      perm->entryType = ResultInfo::SCALAR;
      availResults_.insert( perm );

      // === MAGNETIC ENERGY ===
      shared_ptr<ResultInfo> energy(new ResultInfo);
      energy->resultType = MAG_ENERGY;
      energy->dofNames = "";
      energy->unit = "Ws";
      energy->definedOn = ResultInfo::REGION;
      energy->entryType = ResultInfo::SCALAR;
      availResults_.insert( energy );
      postProcResults_[MAG_ENERGY] = MAG_POTENTIAL;
      shared_ptr<ResultFunctor> energyFunc;
      if( isComplex_ ) {
        energyFunc.reset(new EnergyResultFunctor<Complex>(feFct, energy));
      } else {
        energyFunc.reset(new EnergyResultFunctor<Double>(feFct, energy));
      }
      resultFunctors_[MAG_ENERGY] = energyFunc;
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
//        aFunc->AddIntegrator(bdb, region);
        energyFunc->AddIntegrator(bdb, region);
      }
    }

  // =======================================================================
  //   HELPER METHODS FOR CALCULATING AUXILIARY QUANTITIES
  // =======================================================================

  template<class TYPE>
  void MagEdgePDE::CalcFluxDensityAtIP( const Elem *el,
                                        LocPoint lp,
                                        Vector<TYPE>& field ) {
    Vector<TYPE> tempE, elemSol;
    field.Resize(dim_);
    LocPointMapped lpm;

    shared_ptr<ElemShapeMap> esm = ptgrid_->GetElemShapeMap( el, true );
    shared_ptr<BaseFeFunction> fct = GetFeFunction(MAG_POTENTIAL);

    CurlOperator<FeHCurl,TYPE> li;
    //CurlCurlEdgeInt * li = NULL;

    BaseFE*  ptFe = fct->GetFeSpace()->GetFe( el->elemNum );
    lpm.Set(lp, esm );
    ElemList elist(ptgrid_);
    elist.SetElement(el);
    fct->GetEntitySolution(elemSol, elist.GetIterator() );

    li.ApplyOp( field, lpm, ptFe, elemSol);
  }
  
  template<class TYPE>
  void MagEdgePDE::CalcVecPotentialAtIP( EntityIterator it,
                                         LocPoint lp,
                                         Vector<TYPE>& field ) {
      Vector<TYPE> elemSol;
      field.Resize(dim_);
      LocPointMapped lpm;
      const Elem * el = it.GetElem();
      shared_ptr<ElemShapeMap> esm = ptgrid_->GetElemShapeMap( el, true );
      shared_ptr<BaseFeFunction> fct = GetFeFunction(MAG_POTENTIAL);
                                                     
      FeHCurl*  ptFe = dynamic_cast<FeHCurl*>(fct->GetFeSpace()->GetFe( it ));
      lpm.Set(Elem::shapes[el->type].midPointCoord, esm );
      fct->GetEntitySolution(elemSol, it );
      
      Matrix<Double> shape;
      ptFe->GetShFnc(shape, lpm, el);
      field = shape * elemSol;
  }      
//  template<class TYPE>
//   void MagEdgePDE::CalcEddyCurrentAtIP( EntityIterator it,
//                                         UInt ip,
//                                         Vector<TYPE>& jEddy ) {
//
//     Vector<Double> lCoord;
//     Matrix<Double> shpFnc;
//     Vector<TYPE> magVecDeriv1Elem;
//     Double conductivity = 0.0;
//     
//      // Get the right material parameter for current element
//     RegionIdType actRegionId = it.GetElem()->regionId;
//     materials_[actRegionId]
//       ->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
//     BaseFE * ptEl = it.GetElem()->ptElem;
//
//     Matrix<Double> cornerCoords;
//     ptgrid_->GetElemNodesCoord( cornerCoords, it.GetElem()->connect, 
//                                 true );
//
//     jEddy.Resize(3);
//     
//     // case 1: dummy integration point
//     if( ip == 0 ) {
//       it.GetElem()->ptElem->GetCoordMidPoint( lCoord );
//       ptEl->CalcEdgeShapeFnc( shpFnc, lCoord, cornerCoords, it.GetElem() );
//     } 
//     else {
//       ptEl->CalcEdgeShapeFncAtIp( shpFnc, ip, cornerCoords, it.GetElem() );
//     }
//     
//     // 1) get part coming from vector potential
//     GetDerivSolVecOfElement(magVecDeriv1Elem,it,results_[0]);
//     jEddy.Init();
//     for( UInt iDof = 0; iDof < 3; iDof++ ) {
//       for( UInt i = 0; i < shpFnc.GetNumRows(); i++ ) {
//         jEddy[iDof] += shpFnc[i][iDof] * magVecDeriv1Elem[i];
//       }
//     }
//     jEddy *= -conductivity;   
//  }
//   

  template<class TYPE>
  void MagEdgePDE::CalcPermeability( shared_ptr<BaseResult> result ) {

    TYPE elemPerm;
    Vector<TYPE> elemFlux;
    Double bAbs, reluct;

    // fetch result object and convert data
    Result<TYPE> &  actSol = 
        dynamic_cast<Result<TYPE>&>(*result);
    Vector<TYPE> & actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize());

    // loop over elements
    EntityIterator it = actSol.GetEntityList()->GetIterator();
    for ( it.Begin(); !it.IsEnd(); it++ ) {

      // Determine regionId of element
      const Elem & actEl = *(it.GetElem());
      RegionIdType actRegion = actEl.regionId; 

      //get possible nonlinearities defined in this region
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion]; 

      // Check, if region is nonlinear
      if ( nonLinTypes[actRegion] == PERMEABILITY ) {

        // Obtain nonlinear approximation functional
        ApproxData * approx  = materials_[actRegion]->GetNonlinFncBH(MAG_PERMEABILITY);

        // Calculate flux density in element midpoint
        LocPoint lp;
        lp.coord = Elem::shapes[actEl.type].midPointCoord;
        CalcFluxDensityAtIP( it.GetElem(), lp, elemFlux );
        bAbs = elemFlux.NormL2();
        reluct = approx->EvaluateFuncNu(bAbs);
        elemPerm = 1.0 / reluct;
      } else {
        materials_[actRegion]->GetScalar( elemPerm, MAG_PERMEABILITY,Global::REAL); 
      }
      actVal[it.GetPos() ] = elemPerm;
    }
  }
//  template<class TYPE>
//  void MagEdgePDE::CalcMagPotentialDiv( shared_ptr<BaseResult> result ) {
//    // fetch result object and convert data
//    Result<TYPE> &  actSol = 
//      dynamic_cast<Result<TYPE>&>(*result);
//    EntityIterator it = actSol.GetEntityList()->GetIterator();
//    Vector<TYPE> & actVal = actSol.GetVector();
//    actVal.Resize( actSol.GetEntityList()->GetSize() );
//
//    StdVector<Matrix<Double> >  globDeriv;
//    Matrix<Double> cornerCoords;
//    Vector<Double> elemMidPoint;
//    Vector<TYPE> elemSol;
//    TYPE elemDiv;
//
//    for( it.Begin(); !it.IsEnd(); it++ ) {
//      const Elem * ptElem = it.GetElem();
//
//      // get coordinates of element
//      ptgrid_->GetElemNodesCoord( cornerCoords, ptElem->connect, true );
//
//      // fetch global gradient
//      ptElem->ptElem->GetCoordMidPoint( elemMidPoint );
//      ptElem->ptElem->GetEdgeGlobalDerivShapeFnc( globDeriv, elemMidPoint, 
//                                                  cornerCoords, ptElem );
//
//      // getch solution of element
//      sol_->GetElemSolution(elemSol, it);
//
//      // loop over shape functions
//      elemDiv = 0.0;
//      for( UInt i = 0; i < elemSol.GetSize(); i++ )  {
//
//        elemDiv += elemSol[i] * ( globDeriv[i][0][0] 
//                                + globDeriv[i][1][1] 
//                                + globDeriv[i][2][2] );
//      }
//      actVal[it.GetPos()] = elemDiv;
//    }
//  }

  
  // ======================================================
  // COUPLING SECTION
  // ======================================================


  void MagEdgePDE::InitCoupling(PDECoupling * Coupling) {

    isIterCoupled_ = true;
    ptCoupling_   = Coupling;

    // Enable update of geometry
    const UInt numCouplings = ptCoupling_->GetNumOutputCouplings();


    cplNodeNumPos_.Resize( numCouplings );

    for ( UInt actCoupling = 0; actCoupling < numCouplings; actCoupling++ ){

      if (ptCoupling_->GetOutputQuantity(actCoupling) == MAG_FORCE_LORENTZ) {

        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector(actCoupling,isComplex_);

        //get the element-node to coupling node matching
        StdVector<std::string> couplRegions;
        StdVector<RegionIdType> regionIds;
        ptCoupling_->GetOutputRegions(actCoupling, couplRegions);
        ptgrid_->GetRegion().Parse( couplRegions , regionIds );

        // Check, that every coupling region is part of
        // the magnetic pde itself
        for( UInt iReg = 0; iReg < couplRegions.GetSize(); iReg++ ) {
          if( subdoms_.Find(regionIds[iReg]) == -1 ) {
            EXCEPTION( "Coupling region '" << couplRegions[iReg]
                       << "' is not contained in regions for"
                       << " magnetic PDE" );
          }
        }

        StdVector<UInt> * couplingnodes = NULL;
        ptCoupling_->GetOutputNodes(actCoupling, couplingnodes);

        if ( couplingnodes == NULL ) {
          EXCEPTION( "Pointer couplingnodes = NULL!" );
        }

        for( UInt iNode = 0; iNode < couplingnodes->GetSize(); iNode++ ) {
          UInt actNode = (*couplingnodes)[iNode];
          cplNodeNumPos_[actCoupling][actNode] = iNode;
        }
      }
    }
  }

  void MagEdgePDE::CalcOutputCoupling() {

    SolutionType quantity;
    StdVector<UInt> * couplingNodes = NULL;
    SingleVector * values = NULL;
    UInt forcesCount = 0;

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == false)
      return;

    // loop over all output coupling quantities
    for ( UInt actCoupling = 0;
          actCoupling < ptCoupling_->GetNumOutputCouplings();
          actCoupling++ ) {

      quantity = ptCoupling_->GetOutputQuantity(actCoupling);
      ptCoupling_->GetOutputValues(actCoupling, values);
      Vector<Double> *temp = dynamic_cast<Vector<Double> *>(values);

      switch(ptCoupling_->GetOutputType(actCoupling)) {

      case NODE:
        ptCoupling_->GetOutputNodes(actCoupling, couplingNodes);

        if (quantity == MAG_FORCE_LORENTZ) {
          CalcNodeForceLorentz( *temp, cplNodeNumPos_[forcesCount],
                                actCoupling, couplingNodes->GetSize());
          forcesCount++;
        }

        break;

      case ELEM:
        EXCEPTION( "No Element coupling output" );
      }
    }
  }

  void MagEdgePDE::
  CalcNodeForceLorentz( Vector<Double> & force, 
                        std::map<UInt, UInt> & cplNodeNumPos,
                        UInt actCoupling, UInt numCouplingNodes) {
    EXCEPTION("Not adjusted to new implementation");
    
//    //get the coupling regions
//    StdVector<std::string> couplRegions;
//    StdVector<RegionIdType> regionIds;
//    ptCoupling_->GetOutputRegions(actCoupling, couplRegions);
//    ptgrid_->GetRegion().Parse(couplRegions, regionIds);
//
//   
//    force.Init();
//    Vector<Double>  fluxAtIp(dim_), elemForce, fAtIp, jAtIp;
//    Matrix<Double> ptCoord;
//
//    for (UInt reg=0; reg<couplRegions.GetSize(); reg++) {
//
//      //find subdomain index
//      Integer sdIndex = subdoms_.Find( regionIds[reg] );
//      if( sdIndex == -1 ) {
//        EXCEPTION( "The region coupling region '" <<
//                   ptgrid_->GetRegion().ToString( regionIds[reg] )
//                   << "' was not found in magneticPDE" );
//      }
//      
//      RegionIdType actRegionId = subdoms_[sdIndex] ;
//      Double conductivity;
//      materials_[actRegionId]->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
//      
//      // Check if this region is a coil
//      Integer coilIndex = coilRegionId_.Find(actRegionId);
//      
//      ElemList actSDList(ptgrid_ );
//      actSDList.SetRegion(actRegionId );
//      
//       EntityIterator it = actSDList.GetIterator();
//       UInt actEl = 0;
//       // iterate over all elements of regions
//       for ( it.Begin(); !it.IsEnd(); it++, actEl++ ) {
//         
//         BaseFE * ptElem = it.GetElem()->ptElem;
//
//         //need info for ansatz functions of mechanical FE!
//         shared_ptr<AnsatzFct> fct(new LagrangeFct);
//         ptElem->SetAnsatzFct( fct );
//         UInt numFncs = ptElem->GetNumFncs(fct);
////          ptElem->SetAnsatzFct( results_[0]->fctType );
////          UInt numFncs = ptElem->GetNumFncs(results_[0]->fctType);
//
//         const UInt nrIntPts = ptElem->GetNumIntPoints();
//         const Vector<Double> & intWeights = ptElem->GetIntWeights();  
//         
//         Vector<Double> shpFncAtIp;         
//         ptgrid_->GetElemNodesCoord( ptCoord, 
//                                     it.GetElem()->connect,
//                                     true );
//         
//         // iterate over all integration points
//         fAtIp.Resize( dim_ * numFncs );
//         elemForce.Resize( dim_ * numFncs );
//         elemForce.Init();
//         for( UInt ip = 1; ip <= nrIntPts; ip++ ) {
//           ptElem->GetShFncAtIp(shpFncAtIp, ip, it.GetElem() );
//
//           // CHECK: If this region is a current coil, we simply take the
//           // prescribed current density value
//           if( coilIndex != -1 ) {
//             MathParser * mParser =  domain->GetMathParser();
//             std::string factor = coilDef_[coilIndex]->value_ + "/" 
//               + lexical_cast<std::string>(coilDef_[coilIndex]->windingCrossSection_ );
//             mParser->SetExpr( mHandle_, factor );
//             // TODO: Check if this is still needed
//             // Double currDens = mParser->Eval(mHandle_);
//           } 
//           else {
//             CalcEddyCurrentAtIP( it, 0, jAtIp );
//           }
//           
//           CalcFluxDensityAtIP( it, 0, fluxAtIp );
//           
//           // calculate cross product 
//           fAtIp.Init();
//           Vector<Double> tempCross;
//           jAtIp.CrossProduct( fluxAtIp, tempCross );
//           for (UInt iFnc=0; iFnc<numFncs; iFnc++) {
//             fAtIp[iFnc*3+0] = -tempCross[0] * shpFncAtIp[iFnc];
//             fAtIp[iFnc*3+1] = -tempCross[1] * shpFncAtIp[iFnc];
//             fAtIp[iFnc*3+2] = -tempCross[2] * shpFncAtIp[iFnc];
//           }
//
//           Double jacDet = ptElem->CalcJacobianDetAtIp(ip, ptCoord, 
//                                                       it.GetElem() );
//           elemForce += -fAtIp * (jacDet * intWeights[ip-1] );
//
//         }
//         StdVector<UInt> const & connecth = it.GetElem()->connect;
//
//         // Add the element force to the according coupling node
//         for (UInt ielemnode=0; ielemnode<connecth.GetSize(); ielemnode++) {
//           UInt pos = cplNodeNumPos[connecth[ielemnode]];
//           for( UInt idim=0; idim<dim_; idim++) {
//             force[pos*dim_+idim] += elemForce[ielemnode*dim_+idim];
//           }
//         }
//       }
//       
// 
//     }
  }
  
  std::map<SolutionType, shared_ptr<FeSpace> > 
  MagEdgePDE::CreateFeSpaces(const std::string& formulation,
                             PtrParamNode infoNode ) {
    //ok default case so we create grid based approximation H1 elements
    //and standard Gauss integration
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H_CURL"){
      PtrParamNode potSpaceNode = infoNode->Get("magPotential");
      crSpaces[MAG_POTENTIAL] = 
          FeSpace::CreateInstance(myParam_, potSpaceNode, FeSpace::HCURL );
      crSpaces[MAG_POTENTIAL]->Init();
    }else{
      EXCEPTION("The formulation " << formulation 
                << "of magnetic edge PDE is not known!");
    }
    return crSpaces;
  }

  bool MagEdgePDE::HasOutput( SolutionType output ) {
    if (output == MAG_FORCE_LORENTZ) {
      return true;
    }
    return false;
  }

  LinearFormContext* MagEdgePDE::CreateRhsLinearForm(SolutionType rhsType,shared_ptr<CoefFunction > rhsCoef){

      Exception("Right hand side quantity not known for MagEdgePDE");

      LinearFormContext * mContext = new LinearFormContext(NULL);
      return mContext;
    }

} // end of namespace

