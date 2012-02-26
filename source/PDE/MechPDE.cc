// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MechPDE.hh"

#include <sstream>
#include <iomanip>
#include <set>

#include "General/defs.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Driver/Assemble.hh"


// include fespaces
#include "FeBasis/H1/FeSpaceH1.hh"
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/StrainOperator.hh"
#include "Forms/Operators/IdentityOperatorNormalTrans.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
#include "Forms/Operators/IdentityOperator.hh"

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "CoupledPDE/PDECoupling.hh"

namespace CoupledField {

DECLARE_LOG(mechpde)
DEFINE_LOG(mechpde, "mechpde")

MechPDE::MechPDE(Grid * aptgrid, PtrParamNode paramNode )
    :SinglePDE( aptgrid, paramNode ) {

    pdename_          = "mechanic";
    pdematerialclass_ = MECHANIC;
    maxTimeDerivOrder_ = 2;

    fracDamping_   = false;
    nonLin_        = false;
    nonLinMaterial_= false;

    needSolPrev_ = true;

    // ****************************
    // DETERMINE GEOMETRY
    // ****************************

    // Get problem geometry and PDE subtype
    myParam_->GetValue("subType", subType_ );

    std::string probGeo;
    param->Get("domain")->GetValue("geometryType", probGeo );

    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry
    if ( subType_ == "3d" && probGeo == "3d" ) {
      stressDim_ = 6;
      tensorType_ = FULL;
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = true;
      stressDim_ = 4;
      tensorType_ = AXI;
    }
    else if ( subType_ == "planeStrain" && probGeo == "plane" ) {
      stressDim_ = 3;
      tensorType_ = PLANE_STRAIN;
    }

    else if ( subType_ == "planeStress" && probGeo == "plane" ) {
      stressDim_ = 3;
      tensorType_ = PLANE_STRESS;
    }

    else if ( subType_ == "flatShell" ) {
      stressDim_ = 3;
    }
    else {
      EXCEPTION( "Subtype '" <<  subType_ << "' of PDE '"
                 <<  pdename_ <<  "' does not fit to problem  geometry '"
                 << probGeo << "'"; );
    }
    
    // Sanity check: 3D can only be computed if 3D elements are present
    if( subType_ == "3d" && ptgrid_->GetNumElemOfDim(3) == 0 ) {
      EXCEPTION("Can not calculate 3D mechanics without 3D elements in the grid!");
    }
}


  MechPDE::~MechPDE()
  {

  }

  void MechPDE::ReadDampingInformation( )
  {
    REFACTOR
  }

  void MechPDE::ReadSoftening() {

    REFACTOR;
  }

  void MechPDE::InitNonLin()
  {

    nonLin_ = false;
    REFACTOR;
  }


  void MechPDE::DefineIntegrators() {

    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;

    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[MECH_DISPLACEMENT];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    
    
    //  Loop over all regions
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
      
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;

      // Get current region name
      std::string regionName = ptgrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );
      
      // --------------------------
      //  Set region approximation
      // --------------------------
      // --- Set the approximation for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);
      
      // ====================================================================
      //  Standard Linear Stiffness
      // ====================================================================
      if( !nonLin_ ) {
        BaseBDBInt * stiffInt = 
            GetStiffIntegrator( actSDMat, actRegion, complexMatData_[actRegion] );
        stiffInt->SetName("LinElastInt");
        stiffInt->SetFeSpace( mySpace);
        
        BiLinFormContext * stiffIntDescr =
            new BiLinFormContext(stiffInt, STIFFNESS );
        
        stiffIntDescr->SetEntities( actSDList, actSDList );
        stiffIntDescr->SetFeFunctions( myFct, myFct );
        
        assemble_->AddBiLinearForm( stiffIntDescr );
        
        // Important: Add bdb-integrator to global list, as we need them later
        // for calculation of postprocessing results
        bdbInts_[actRegion] = stiffInt;
        
      }
      
      // ====================================================================
      //  Standard Mass Integrator
      // ====================================================================
      Double density = 0.0;
      actSDMat->GetScalar(density,DENSITY,Global::REAL);
      shared_ptr<CoefFunction> densCoeff 
          = CoefFunction::Generate(Global::REAL, lexical_cast<std::string>(density));
      
      BiLinearForm *massInt = NULL;
      if( dim_ == 2 ) {
        massInt = new BBInt<IdentityOperator<FeH1,2,2> >(densCoeff, 1.0);
      } else {
        massInt = new BBInt<IdentityOperator<FeH1,3,3> >(densCoeff, 1.0);
      }
      massInt->SetName("MassInt");
      massInt->SetFeSpace( mySpace );

      BiLinFormContext *massContext =  new BiLinFormContext( massInt, MASS );
      massContext->SetEntities( actSDList, actSDList );
      massContext->SetFeFunctions( myFct, myFct );
      assemble_->AddBiLinearForm( massContext );

      
      // in the end, at the region to the feFunction      // to be implemented
      myFct->AddEntityList( actSDList );
    }
  }
  
  void MechPDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(mechpde) << "Defining rhs load integrators for mechanic PDE";
    
    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[MECH_DISPLACEMENT];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<shared_ptr<CoefFunction> > coef;
    LinearForm * lin = NULL;
    StdVector<std::string> dispDofNames = myFct->GetResultInfo()->dofNames;
    
    // ========================
    //  FORCES (volume, nodal)
    // ========================
    LOG_DBG(mechpde) << "Reading forces";
    ReadRhsExcitation( "force", dispDofNames, ResultInfo::VECTOR, isComplex_, ent, coef );
    
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        
        // --------------
        //  Nodal Forces 
        // --------------
        // In case of a nodal force, we divide the coeffunction by the number
        // of nodes to ensure, that the total force is prescribed. This works
        // only for constant coefficient functions
        if( coef[i]->GetDependency() == CoefFunction::GENERAL ) {
          EXCEPTION("Nodal forces must not be spatial dependent");
        }
        
        UInt numNodes = ent[i]->GetSize();
        if( numNodes > 1 ) {
          // Here we would divide the nodal force by the number of nodes
          // in the list, in order to ensure that the whole force corresponds
          // to the prescribed value. However, this requires modification of 
          // the expressions of the coefficient functions, which depends on real/harm
          // and the specific type (const, timefreq, variable).
          WARN("The force value will not be divided by the number of nodes and thus "
              << "depends on the number of nodes" );
        }
        
        lin = new SingleEntryInt(coef[i]);
        lin->SetName("NodalForceInt");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
        
      } else {
        
        // -------------------------
        //  Surface / Volume Forces 
        // -------------------------
        EXCEPTION("Not yet implemented")
    
        // Same issue here as above: We need to "divide" the total force by the
        // area / volume to get the force density.
      }
    } // for
    
    // ===============
    //  FORCE DENSITY 
    // ===============
    LOG_DBG(mechpde) << "Reading force densities";
    
    ReadRhsExcitation( "forceDensity", dispDofNames, ResultInfo::VECTOR, isComplex_, ent, coef );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Force density must be defined on elements")
      }
      
      // determine dimension
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != dim_ ) {
        EXCEPTION("Foce density can only be defined on surface elements");
      }
      
      if( dim_ == 2) {
        if(isComplex_) {
          lin = new BUIntegrator<IdentityOperator<FeH1,2,2>, Complex>(Complex(1.0), coef[i]);
        } else {
          lin = new BUIntegrator<IdentityOperator<FeH1,2,2>, Double>(1.0, coef[i]);
        }
      } else  {
        if(isComplex_) {
          lin = new BUIntegrator<IdentityOperator<FeH1,3,3>, Complex>(Complex(1.0), coef[i]);
        } else {
          lin = new BUIntegrator<IdentityOperator<FeH1,3,3>, Double>(1.0, coef[i]);
        }
      }
      lin->SetName("ForceDensityInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    } // for
    
    
    // ===============
    //  PRESSURE 
    // ===============
    LOG_DBG(mechpde) << "Reading mechanical pressure";
    StdVector<std::string> empty;
    ReadRhsExcitation( "pressure", empty, ResultInfo::VECTOR, isComplex_, ent, coef );
    std::set<RegionIdType> volRegions (subdoms_.Begin(), subdoms_.End() );
    
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Mechanical pressure must be defined on elements")
      }
            
      if( dim_ == 2) {
        if(isComplex_) {
          lin = new BUIntegrator<IdentityOperatorNormalTrans<FeH1,2>, 
                                 Complex,true>(Complex(1.0), coef[i], volRegions);
        } else {
          lin = new BUIntegrator<IdentityOperatorNormalTrans<FeH1,2>, 
                                 Double,true>(1.0, coef[i], volRegions);
        }
      } else  {
        if(isComplex_) {
          lin = new BUIntegrator<IdentityOperatorNormalTrans<FeH1,3>, 
                                 Complex, true>(Complex(1.0), coef[i], volRegions);
        } else {
          lin = new BUIntegrator<IdentityOperatorNormalTrans<FeH1,3>, 
                                 Double, true>(1.0, coef[i], volRegions);
        }
      }
      lin->SetName("PressureInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    } // for

    
    // ==================
    //  SURFACE TRACTION  
    // ==================
    LOG_DBG(mechpde) << "Reading surface tractions";
      
      ReadRhsExcitation( "traction", dispDofNames, ResultInfo::VECTOR, isComplex_, ent, coef );
      for( UInt i = 0; i < ent.GetSize(); ++i ) {
        // check type of entitylist
        if (ent[i]->GetType() == EntityList::NODE_LIST) {
          EXCEPTION("Surface traction must be defined on elements")
        }
        // ensure that list contains only surface elements
        EntityIterator it = ent[i]->GetIterator();
        UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
        if( elemDim != (dim_-1) ) {
          EXCEPTION("Surface traction can only be defined on surface elements");
        }
        
        if( dim_ == 2) {
          if(isComplex_) {
            lin = new BUIntegrator<IdentityOperator<FeH1,2,2>, Complex>(Complex(1.0), coef[i]);
          } else {
            lin = new BUIntegrator<IdentityOperator<FeH1,2,2>, Double>(1.0, coef[i]);
          }
        } else  {
          if(isComplex_) {
            lin = new BUIntegrator<IdentityOperator<FeH1,3,3>, Complex>(Complex(1.0), coef[i]);
          } else {
            lin = new BUIntegrator<IdentityOperator<FeH1,3,3>, Double>(1.0, coef[i]);
          }
        }
        lin->SetName("TractionIntegrator");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
        myFct->AddEntityList(ent[i]);
      } // for
    
    
  }

  BaseBDBInt *
  MechPDE::GetStiffIntegrator( BaseMaterial* actSDMat,
                               RegionIdType regionId,
                               bool isComplex )
  {

    // Get region name
    std::string regionName = ptgrid_->GetRegion().ToString( regionId );

    // ------------------------
    //  Obtain linear material
    // ------------------------
    shared_ptr<CoefFunction > curCoef;
    if( isComplex ) {
      curCoef = actSDMat->GetCoefFunction(MECH_STIFFNESS_TENSOR,
                                          tensorType_, Global::COMPLEX, false);
    } else {
      curCoef = actSDMat->GetCoefFunction(MECH_STIFFNESS_TENSOR,
                                          tensorType_, Global::REAL, false);
    }
    
    // store coefficient function for later use (e.g. in boundary integrators)
    regionStiffness_[regionId] = curCoef;
    
    // ----------------------------------------
    //  Determine correct stiffness integrator 
    // ----------------------------------------
    BaseBDBInt * integ = NULL;
    if( subType_ == "axi" ) {
      integ = new BDBInt<StrainOperatorAxi<FeH1> >(curCoef, 1.0);
    } else if( subType_ == "planeStrain" ) {
      integ = new BDBInt<StrainOperator2D<FeH1> >(curCoef, 1.0);
    } else if( subType_ == "planeStress" ) {
      EXCEPTION("Not implemented");
    } else if( subType_ == "3d") {
      integ = new BDBInt<StrainOperator3D<FeH1> >(curCoef, 1.0);
    } else {
      EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
    }
    return integ;
  }

  void MechPDE::DefineSolveStep()
  {
		  solveStep_ = new StdSolveStep(*this);
  }



  // ======================================================
  // COUPLING SECTION
  // ======================================================


  void MechPDE::InitCoupling(PDECoupling * Coupling)
  {

    isIterCoupled_ = true;
    ptCoupling_   = Coupling;

    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
      {
        if (ptCoupling_->GetOutputQuantity(i) == MECH_DISPLACEMENT)
          {
            // Intialize the memory of the coupling values
            ptCoupling_->CreateCouplingVector(i,isComplex_);
          }

        if (ptCoupling_->GetOutputQuantity(i) == MECH_VELOCITY)
          {
            // Intialize the memory of the coupling values
            ptCoupling_->CreateCouplingVector(i,isComplex_);
          }

        if (ptCoupling_->GetOutputQuantity(i) == MECH_FORCE)
          {
            // Intialize the memory of the coupling values
            ptCoupling_->CreateCouplingVector(i,isComplex_);

            //now since we need a incremental formulation, initialize some necessary vectors
            isIncrFormulation_ = true;
          }
      }

  }


  void MechPDE::CalcOutputCoupling()
  {
    REFACTOR;
  }

  bool MechPDE::HasOutput(SolutionType output)
  {

    if (output == MECH_DISPLACEMENT || output == MECH_VELOCITY || output == MECH_FORCE)
      return true;

    return false;
  }



  void MechPDE::CalcResults( shared_ptr<BaseResult> res ) {

    switch (res->GetResultInfo()->resultType ) {

      case MECH_DISPLACEMENT:
        feFunctions_[MECH_DISPLACEMENT]->ExtractResult( res );
        break;

      case MECH_VELOCITY:
        if ( (analysistype_ == TRANSIENT) || (analysistype_ == HARMONIC)){
          this->timeDerivFeFunctions_[MECH_VELOCITY]->ExtractResult( res );
        }
        break;

      case MECH_ACCELERATION:
        if ( (analysistype_ == TRANSIENT) || (analysistype_ == HARMONIC)){
          this->timeDerivFeFunctions_[MECH_ACCELERATION]->ExtractResult( res );
        }
        break;

      case MECH_RHS_LOAD:
        rhsFeFunctions_[MECH_DISPLACEMENT]->ExtractResult( res );
        break;

      case MECH_STRESS:
        resultFunctors_[MECH_STRESS]->EvalResult( res );
        break;

      case MECH_STRAIN:
        resultFunctors_[MECH_STRAIN]->EvalResult( res );
        break;

      default:
        WARN( "Result '" << 
              SolutionTypeEnum.ToString(res->GetResultInfo()->resultType)
              << "' type not computable by mechanic PDE" );
        break;
    }
  }

  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================
  void MechPDE::InitTimeStepping()
  {
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(TimeSchemeGLM::NEWMARK, 0) );
    feFunctions_[MECH_DISPLACEMENT]->SetTimeScheme(myScheme);
  }

  //Define time derivative feFunctions
  //the idea is the folowing:
  // in the case of a transient analysis, we can directly relate
  // the coefficient vector of the feFunction to the one of the
  // time stepping scheme.
  // in case of harmonics we create a new coefVector for the functions
  // which has to be updated within the harmonic driver
  void MechPDE::DefineTimeDerivFeFunctions(){

    if ( (analysistype_ != TRANSIENT) && (analysistype_ != HARMONIC)){
      //ok for the static case we do not need to bother
      return;
    }
    //first we find the result infos
    shared_ptr<ResultInfo> velInfo;
    shared_ptr<ResultInfo> accelInfo;
    ResultSet::const_iterator it = availResults_.begin();
    for( ; it != availResults_.end(); ++it ) {
      if( (*it)->resultType == MECH_VELOCITY ) {
        velInfo = *it;
      }else if( (*it)->resultType == MECH_ACCELERATION ){
        accelInfo = *it;
      }
    }
    if(!velInfo || !accelInfo){
      EXCEPTION("Could not determine the result info for time derivative FeFunctions.")
    }
    //TODO: This is not very nice so we need to implement a copy method to fefunction.

    if (!isComplex_)
    {
      this->timeDerivFeFunctions_[MECH_VELOCITY].reset(new FeFunction<Double>());
      this->timeDerivFeFunctions_[MECH_ACCELERATION].reset(new FeFunction<Double>());

    }else
    {
      this->timeDerivFeFunctions_[MECH_VELOCITY].reset(new FeFunction<Complex>());
      this->timeDerivFeFunctions_[MECH_ACCELERATION].reset(new FeFunction<Complex>());
    }
    shared_ptr<BaseFeFunction> dispFe = feFunctions_[MECH_DISPLACEMENT];
    shared_ptr<BaseFeFunction> vFct = timeDerivFeFunctions_[MECH_VELOCITY];
    shared_ptr<BaseFeFunction> aFct = timeDerivFeFunctions_[MECH_ACCELERATION];

    //lets create now a copy and add entity lists and stuff
    //TODO: Make there a copy function to do such work....
    StdVector< shared_ptr<EntityList> > entList = dispFe->GetEntityList();
    for(UInt i =0;i<entList.GetSize();i++){
      vFct->AddEntityList(entList[i]);
      aFct->AddEntityList(entList[i]);
    }

    //Copy info for velocity
    vFct->SetResultInfo(velInfo);
    vFct->SetFeSpace(dispFe->GetFeSpace());
    vFct->SetFctId(dispFe->GetFctId());
    vFct->Finalize();

    aFct->SetResultInfo(accelInfo);
    aFct->SetFeSpace(dispFe->GetFeSpace());
    aFct->SetFctId(dispFe->GetFctId());
    aFct->Finalize();

    //the only thing left to do is the association of corresponding coefVectors in the transient case
    if ( analysistype_ == TRANSIENT)
    {
      dispFe->GetTimeScheme()->SetTimeDerivVector(1,vFct->GetSingleVector());
      dispFe->GetTimeScheme()->SetTimeDerivVector(2,aFct->GetSingleVector());
    }

    //in any case we create a Field Functor for our functions

    shared_ptr<BaseFieldFunctor> vFunc;
    shared_ptr<BaseFieldFunctor> aFunc;
    if( isComplex_ ) {
        vFunc.reset(
            new FieldInterpolFunctor<IdentityOperator<FeH1,2,2,Complex>,
            Complex>(vFct, velInfo));
        aFunc.reset(
            new FieldInterpolFunctor<IdentityOperator<FeH1,2,2,Complex>,
           Complex>(aFct, accelInfo));

    } else {
        vFunc.reset(
            new FieldInterpolFunctor<IdentityOperator<FeH1,2,2,Double>,
            Double>(vFct, velInfo));
        aFunc.reset(
            new FieldInterpolFunctor<IdentityOperator<FeH1,2,2,Double>,
            Double>(aFct, accelInfo));
    }
    resultFunctors_[MECH_VELOCITY] = vFunc;
    fieldFunctors_[MECH_VELOCITY] = vFunc;
    resultFunctors_[MECH_ACCELERATION] = aFunc;
    fieldFunctors_[MECH_ACCELERATION] = aFunc;
  }




 
  void MechPDE::DefinePrimaryResults()
  {
    // Check for subType
    StdVector<std::string> dispDofNames, stressDofNames;

    if( subType_ == "3d" ) {
      dispDofNames = "x", "y", "z";
      stressDofNames = "xx", "yy", "zz", "yz", "xz", "xy";

    } else if( subType_ == "planeStrain" ) {
      dispDofNames = "x", "y";
      stressDofNames = "xx", "yy", "xy";

    } else if( subType_ == "planeStress" ) {
      dispDofNames = "x", "y";
      stressDofNames = "xx", "yy", "xy";

    } else if( subType_ == "axi" ) {
      dispDofNames = "r", "z";
      stressDofNames = "rr", "zz", "rz", "phiphi";
    }

    // === MECHANIC DISPLACEMENT ===
    shared_ptr<ResultInfo> disp(new ResultInfo);
    disp->resultType = MECH_DISPLACEMENT;
    disp->dofNames = dispDofNames;
    disp->unit = "m";
    disp->entryType = ResultInfo::VECTOR;
    disp->SetFeFunction(feFunctions_[MECH_DISPLACEMENT]);
    disp->definedOn = ResultInfo::NODE;
    feFunctions_[MECH_DISPLACEMENT]->SetResultInfo(disp);
    
    // this defines the primary unknown
    results_.Push_back( disp );
    availResults_.insert( disp);
    
    // define functor for interpolation of the field
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MECH_DISPLACEMENT];
    shared_ptr<BaseFieldFunctor> dFunc;
       if( isComplex_ ) {
         dFunc.reset(
             new FieldInterpolFunctor<IdentityOperator<FeH1,2,2,Complex>,
             Complex>(feFct, disp));
       } else {
         dFunc.reset(
             new FieldInterpolFunctor<IdentityOperator<FeH1,2,2,Double>,
             Double>(feFct, disp));
       }
       resultFunctors_[MECH_DISPLACEMENT] = dFunc;
       fieldFunctors_[MECH_DISPLACEMENT] = dFunc;

    if ( (analysistype_ == TRANSIENT) || (analysistype_ == HARMONIC))
    {
      // === MECHANIC VELOCITY ===
      shared_ptr<ResultInfo> vel(new ResultInfo);
      vel->resultType = MECH_VELOCITY;
      vel->dofNames = dispDofNames;
      vel->unit = "m/s";
      vel->entryType = disp->entryType;
      vel->definedOn = disp->definedOn;
      availResults_.insert( vel );

      // === MECHANIC ACCELERATION ===
      shared_ptr<ResultInfo> acc(new ResultInfo);
      acc->resultType = MECH_ACCELERATION;
      acc->dofNames = dispDofNames;
      acc->unit = "m/s^2";
      acc->entryType = disp->entryType;
      acc->definedOn = disp->definedOn;
      availResults_.insert( acc );
    }

    // === MECHANIC RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = MECH_RHS_LOAD;
    rhs->dofNames = dispDofNames;
    rhs->unit = "N";
    rhs->entryType = disp->entryType;
    rhs->definedOn = disp->definedOn;
    availResults_.insert( rhs );
    postProcResults_[MECH_RHS_LOAD] = MECH_DISPLACEMENT;
    //
    
  }
    
  void MechPDE::DefinePostProcResults() {
    StdVector<std::string> stressComponents;
    if( subType_ == "3d" ) {
      stressComponents = "xx", "yy", "zz", "yz", "xz", "xy";
    } else if( subType_ == "planeStrain" ) {
      stressComponents = "xx", "yy", "xy";
    } else if( subType_ == "planeStress" ) {
      stressComponents = "xx", "yy", "xy";
    } else if( subType_ == "axi" ) {
      stressComponents = "rr", "zz", "rz", "phiphi";
    }
    
    
    StdVector<std::string > dispDofNames;
    dispDofNames = feFunctions_[MECH_DISPLACEMENT]->GetResultInfo()->dofNames;
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MECH_DISPLACEMENT];
    
    // === MECHANIC STRESS ===
    shared_ptr<ResultInfo> stress(new ResultInfo);
    stress->resultType = MECH_STRESS;
    stress->dofNames = stressComponents;
    stress->unit =  "N/m^2";
    stress->entryType = ResultInfo::TENSOR;
    stress->definedOn = ResultInfo::ELEMENT;
    availResults_.insert( stress );
    postProcResults_[MECH_STRESS] = MECH_DISPLACEMENT;
    shared_ptr<BaseFieldFunctor> sigmaFunc;
    if( isComplex_ ) {
      sigmaFunc.reset(new FluxFieldFunctor<Complex>(feFct, stress));
    } else {
      sigmaFunc.reset(new FluxFieldFunctor<Double>(feFct, stress));
    }
    resultFunctors_[MECH_STRESS] = sigmaFunc;
    fieldFunctors_[MECH_STRESS] = sigmaFunc;
       
    
    
    // === MECHANIC STRAIN ===
    shared_ptr<ResultInfo> strain(new ResultInfo);
    strain->resultType = MECH_STRAIN;
    strain->dofNames = stressComponents;
    strain->unit =  "";
    strain->entryType = ResultInfo::TENSOR;
    strain->definedOn = ResultInfo::ELEMENT;
    availResults_.insert( strain );
    postProcResults_[MECH_STRAIN] = MECH_DISPLACEMENT;
    shared_ptr<BaseFieldFunctor> strainFunc;
    if( isComplex_ ) {
      strainFunc.reset(new DiffFieldFunctor<Complex>(feFct, strain));
    } else {
      strainFunc.reset(new DiffFieldFunctor<Double>(feFct, strain));
    }
    resultFunctors_[MECH_STRAIN] = strainFunc;
    fieldFunctors_[MECH_STRAIN] = strainFunc;

    // ============================
    // Initialize result functors:
    // ============================
    // 1) Loop over all BDB-integrators
    std::map<RegionIdType, BaseBDBInt*>::iterator it = bdbInts_.begin();
    for( ; it != bdbInts_.end(); ++it ) {
      RegionIdType region = it->first;
      BaseBDBInt* bdb = it->second;

      // 2) pass integrators to functors
      sigmaFunc->AddIntegrator(bdb, region);
      strainFunc->AddIntegrator(bdb, region);
    }
    
//
//    // === MECHANIC PRESSURE ===
//    shared_ptr<ResultInfo> pressure(new ResultInfo);
//    pressure->resultType = MECH_PRESSURE;
//    pressure->dofNames = "";
//    pressure->unit = "N/m^2";
//    pressure->entryType = ResultInfo::SCALAR;
//    pressure->definedOn = ResultInfo::SURF_ELEM;
//    pressure->fctType = shared_ptr<ConstFct>( new ConstFct() );
//    availResults_.insert( pressure );
//
//    // === MECHANIC VON MISES STRESS (yield criterion, a scalar value)===
//    shared_ptr<ResultInfo> vonMises(new ResultInfo);
//    vonMises->resultType = VON_MISES_STRESS;
//    vonMises->dofNames = "";
//    vonMises->unit =  "";
//    vonMises->entryType = ResultInfo::SCALAR;
//    vonMises->definedOn = ResultInfo::ELEMENT;
//    vonMises->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert( vonMises );
//
//
//    // === MECHANIC VON MISES STRAIN (yield criterion, a scalar value)===
//    shared_ptr<ResultInfo> vonMisesStrain(new ResultInfo);
//    vonMisesStrain->resultType = VON_MISES_STRAIN;
//    vonMisesStrain->dofNames = "";
//    vonMisesStrain->unit =  "";
//    vonMisesStrain->entryType = ResultInfo::SCALAR;
//    vonMisesStrain->definedOn = ResultInfo::ELEMENT;
//    vonMisesStrain->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert( vonMisesStrain );
//
//
//

//
//    // === MECHANIC ENERGY ===
//    shared_ptr<ResultInfo> energy(new ResultInfo);
//    energy->resultType = MECH_ENERGY;
//    energy->dofNames = "";
//    energy->unit = "Ws";
//    energy->entryType = ResultInfo::SCALAR;
//    energy->definedOn = ResultInfo::REGION;
//    energy->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert( energy );
//
//    // === LUMPED_MECHANIC DISPLACEMENT ===
//    shared_ptr<ResultInfo> elem_disp(new ResultInfo);
//    elem_disp->resultType = LUMPED_MECH_DISPLACEMENT;
//    elem_disp->dofNames = dispDofNames;
//    elem_disp->unit = "m";
//    if(subType_ != "flatShell") {
//      elem_disp->entryType = ResultInfo::VECTOR;
//    } else {
//      elem_disp->entryType = ResultInfo::TENSOR;
//    }
//    elem_disp->definedOn = ResultInfo::ELEMENT;
//    elem_disp->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert(elem_disp);
//
//    // === PSEUDO DENSITY for SIMP ===
//    shared_ptr<ResultInfo> mechPD(new ResultInfo);
//    mechPD->resultType = MECH_PSEUDO_DENSITY;
//    mechPD->dofNames = "";
//    mechPD->unit = "";
//    mechPD->entryType = ResultInfo::SCALAR;
//    mechPD->definedOn = ResultInfo::ELEMENT;
//    mechPD->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert( mechPD );
//
//    shared_ptr<ResultInfo> pysicalPD(new ResultInfo);
//    pysicalPD->resultType = PHYSICAL_PSEUDO_DENSITY;
//    pysicalPD->dofNames = "";
//    pysicalPD->unit = "";
//    pysicalPD->entryType = ResultInfo::SCALAR;
//    pysicalPD->definedOn = ResultInfo::ELEMENT;
//    pysicalPD->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert( pysicalPD );
//
//
//    // === SHAPE offset ===
//    shared_ptr<ResultInfo> shape(new ResultInfo);
//    shape->resultType = MECH_SHAPE;
//    shape->dofNames = dispDofNames;
//    shape->unit = "m";
//    shape->entryType = disp->entryType;
//    shape->definedOn = disp->definedOn;
//    shape->fctType = disp->fctType;
//    availResults_.insert( shape );
//
//    // === HOMOGENIZED_TENSOR ===
//    // calculated as a special optimization result
//    shared_ptr<ResultInfo> homogenTensor(new ResultInfo);
//    homogenTensor->resultType = HOMOGENIZED_TENSOR;
//    homogenTensor->dofNames = "";
//    homogenTensor->unit = "";
//    homogenTensor->entryType = ResultInfo::TENSOR;
//    homogenTensor->definedOn = ResultInfo::REGION;
//    homogenTensor->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert(homogenTensor);
//
//
//    // gradDisplacement as a nodal property
//    // TODO a method would reduce the copy and paste!
//    shared_ptr<ResultInfo> grad_x_displ(new ResultInfo);
//    grad_x_displ->resultType = GRAD_X_DISPLACEMENT;
//    grad_x_displ->SetVectorDOFs(dim_, isaxi_);
//    grad_x_displ->unit = "m/m";
//    grad_x_displ->entryType = ResultInfo::VECTOR;
//    grad_x_displ->definedOn = ResultInfo::NODE;
//    grad_x_displ->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert(grad_x_displ);
//
//    shared_ptr<ResultInfo> grad_y_displ(new ResultInfo);
//    grad_y_displ->resultType = GRAD_Y_DISPLACEMENT;
//    grad_y_displ->SetVectorDOFs(dim_, isaxi_);
//    grad_y_displ->unit = "m/m";
//    grad_y_displ->entryType = ResultInfo::VECTOR;
//    grad_y_displ->definedOn = ResultInfo::NODE;
//    grad_y_displ->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert(grad_y_displ);
//
//    shared_ptr<ResultInfo> grad_z_displ(new ResultInfo);
//    grad_z_displ->resultType = GRAD_Z_DISPLACEMENT;
//    grad_z_displ->SetVectorDOFs(dim_, isaxi_);
//    grad_z_displ->unit = "m/m";
//    grad_z_displ->entryType = ResultInfo::VECTOR;
//    grad_z_displ->definedOn = ResultInfo::NODE;
//    grad_z_displ->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert(grad_z_displ);
  }
  
  
  std::map<SolutionType, shared_ptr<FeSpace> >
   MechPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
     
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    
    if( formulation == "default" || formulation == "H1" ){
      PtrParamNode potSpaceNode = infoNode->Get("mechDisplacement");
      crSpaces[MECH_DISPLACEMENT] =
          FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1);
      crSpaces[MECH_DISPLACEMENT]->Init(solStrat_);

    }else{
       EXCEPTION( "The formulation " << formulation 
                  << "of the mechanic PDE is not known!" );
     }
     return crSpaces;
   }

} // end namespace CoupledField
