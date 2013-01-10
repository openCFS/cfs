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


// include elements
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormalTrans.hh"
#include "Forms/Operators/StrainOperator.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"

#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
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
    if( subType_ == "3d" && ptGrid_->GetNumElemOfDim(3) == 0 ) {
      EXCEPTION("Can not calculate 3D mechanics without 3D elements in the grid!");
    }
}


  MechPDE::~MechPDE()
  {

  }

  void MechPDE::ReadDampingInformation( ) {

//    bool identical = true; // i.e. same type of damping for all regions

    std::map<std::string, DampingType> idDampType;
    std::map<std::string, shared_ptr<RaylDampingData> > idRaylData;

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

        if( actType == RAYLEIGH ) {
          // set data for Rayleigh damping
          shared_ptr<RaylDampingData> actRaylDamp(new RaylDampingData());
          actRaylDamp->alpha = "0.0";
          actRaylDamp->beta = "0.0";
          actRaylDamp->adjustDamping = true;
          actRaylDamp->ratioDeltaF = 0.01;
          actRaylDamp->freq = 0.0;

          dampNodes[i]->GetValue( "freq", actRaylDamp->freq, ParamNode::PASS);
          dampNodes[i]->GetValue( "ratioDeltaF", actRaylDamp->ratioDeltaF, ParamNode::PASS );
          dampNodes[i]->GetValue( "adjustDamping", actRaylDamp->adjustDamping, ParamNode::PASS );
          idRaylData[actId] = actRaylDamp;
        }

        // store damping type string
        idDampType[actId] = actType;
      }
    }

    // Run over all region and set entry in "regionNonLinId"
    ParamNodeList regionNodes =
      myParam_->Get("regionList")->GetChildren();

    RegionIdType actRegionId;
    std::string actRegionName, actDampingId;

//    if( regionNodes.GetSize() > 0 ) {
//      Info->PrintF( pdename_, "Damping in following region(s)\n" );
//    }

    for (UInt k = 0; k < regionNodes.GetSize(); k++) {
      regionNodes[k]->GetValue( "name", actRegionName );
      regionNodes[k]->GetValue( "dampingId", actDampingId );
      if( actDampingId == "" )
        continue;

      actRegionId = ptGrid_->GetRegion().Parse( actRegionName );

      // Check actDampingId was already registerd
      if( idDampType.count( actDampingId ) == 0 ) {
        EXCEPTION( "Damping with id '" << actDampingId
                   << "' was not defined in 'dampingList'" );
      }

      dampingList_[actRegionId] = idDampType[actDampingId];
      if ( dampingList_[actRegionId] == RAYLEIGH ){
        RaylDampingData actRayl = *(idRaylData[actDampingId]); 
        Double dampFreq;
        
        if( actRayl.freq == 0.0 ) {
          materials_[actRegionId]->GetScalar(dampFreq,RAYLEIGH_FREQUENCY,Global::REAL);
        } else { 
          dampFreq = actRayl.freq;
        }
        
        // Check if analysis is harmonic
        bool isHarmonic = analysistype_ == BasePDE::HARMONIC;

        // Compute Rayleigh damping parameters
        materials_[actRegionId]->
         ComputeRayleighDamping( actRayl.alpha, actRayl.beta,
                                 dampFreq, actRayl.ratioDeltaF, 
                                 actRayl.adjustDamping, isHarmonic );
        regionRaylDamping_[actRegionId] = actRayl;

        PtrParamNode in = infoNode_->Get(ParamNode::HEADER)->GetByVal("region", "name", domain->GetGrid()->GetRegion().ToString(actRegionId));
        in->Get("alpha_M")->SetValue(actRayl.alpha);
        in->Get("alpha_K")->SetValue(actRayl.beta);
      }
    }

    // Check, if all entries are identical
    for ( UInt i = 1; i < dampingList_.size(); i++ ) {
      if ( dampingList_[regions_[i-1]] != dampingList_[regions_[i]] ) {
//        identical = false;
        break;
      }
    }
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
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
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
        
        //check for damping
        if ( dampingList_[actRegion] == RAYLEIGH ) {
          RaylDampingData & actDamp = (regionRaylDamping_[actRegion]);
          stiffIntDescr->SetSecDestMat(DAMPING, actDamp.beta );
        }
        
        assemble_->AddBiLinearForm( stiffIntDescr );
        
        // Important: Add bdb-integrator to global list, as we need them later
        // for calculation of postprocessing results
        bdbInts_[actRegion] = stiffInt;
        
      }
      
      // ====================================================================
      //  Standard Mass Integrator
      // ====================================================================
      PtrCoefFct densCoeff = actSDMat->GetScalCoefFnc( DENSITY,Global::REAL );
      
      BaseBDBInt *massInt = NULL;
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
      
      // Check for damping (mass part)
      if ( dampingList_[actRegion] == RAYLEIGH ) {
        RaylDampingData & actDamp = regionRaylDamping_[actRegion];
        massContext->SetSecDestMat( DAMPING, actDamp.alpha );
      }
      
      // Important: Add mass-integrator to global list, as we need them later
      // for calculation of postprocessing results
      massInts_[actRegion] = massInt;
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
    StdVector<PtrCoefFct > coef;
    LinearForm * lin = NULL;
    StdVector<std::string> dispDofNames = myFct->GetResultInfo()->dofNames;
    
    // ========================
    //  FORCES (volume, nodal)
    // ========================
    LOG_DBG(mechpde) << "Reading forces";
    ReadRhsExcitation( "force", dispDofNames, ResultInfo::VECTOR, isComplex_,
                       ent, coef );
    
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      
      // In case of a total force, we can not have a spatial dependency
      if( coef[i]->GetDependency() == CoefFunction::GENERAL ) {
        EXCEPTION("Total forces must not be spatial dependent");
      }
              
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        
        // --------------
        //  Nodal Forces 
        // --------------
        UInt numNodes = ent[i]->GetSize();
        // If more than one node is defined, we divide the total force by the number
        // of nodes to ensure that the total force is applied, independent of the 
        // number of nodes
        if( numNodes > 1 ) {
          Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;  
          coef[i] = CoefFunction::Generate(part, CoefXprVecScalOp(coef[i], 
                    boost::lexical_cast<std::string>(numNodes), CoefXpr::OP_DIV) );
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
    
    ReadRhsExcitation( "forceDensity", dispDofNames, ResultInfo::VECTOR, isComplex_, 
                        ent, coef );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Force density must be defined on elements")
      }
      
      // determine dimension
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != dim_-1 ) {
        EXCEPTION("Force density can only be defined on surface elements");
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
    ReadRhsExcitation( "pressure", empty, ResultInfo::VECTOR, isComplex_, 
                        ent, coef );
    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );
    
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
      
      ReadRhsExcitation( "traction", dispDofNames, ResultInfo::VECTOR, isComplex_, 
                          ent, coef );
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
    std::string regionName = ptGrid_->GetRegion().ToString( regionId );

    // ------------------------
    //  Obtain linear material
    // ------------------------
    shared_ptr<CoefFunction > curCoef;
    if( isComplex ) {
      curCoef = actSDMat->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR,
                                          tensorType_, Global::COMPLEX );
    } else {
      curCoef = actSDMat->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR,
                                          tensorType_, Global::REAL );
    }
    
    // store coefficient function for later use (e.g. in boundary integrators)
    regionStiffness_[regionId] = curCoef;
    
    // ----------------------------------------
    //  Determine correct stiffness integrator 
    // ----------------------------------------
    BaseBDBInt * integ = NULL;
    if( isComplex ) {
      if( subType_ == "axi" ) {
        integ = new BDBInt<StrainOperatorAxi<FeH1,Complex>,Complex,Complex>(curCoef, 1.0);
      } else if( subType_ == "planeStrain" ) {
        integ = new BDBInt<StrainOperator2D<FeH1,Complex>,Complex,Complex>(curCoef, 1.0);
      } else if( subType_ == "planeStress" ) {
        integ = new BDBInt<StrainOperator2D<FeH1,Complex>,Complex,Complex>(curCoef, 1.0);
      } else if( subType_ == "3d") {
        integ = new BDBInt<StrainOperator3D<FeH1,Complex>,Complex,Complex>(curCoef, 1.0);
      } else {
        EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
      }
    }
    else {
      if( subType_ == "axi" ) {
        integ = new BDBInt<StrainOperatorAxi<FeH1> >(curCoef, 1.0);
      } else if( subType_ == "planeStrain" ) {
        integ = new BDBInt<StrainOperator2D<FeH1> >(curCoef, 1.0);
      } else if( subType_ == "planeStress" ) {
        integ = new BDBInt<StrainOperator2D<FeH1> >(curCoef, 1.0);
      } else if( subType_ == "3d") {
        integ = new BDBInt<StrainOperator3D<FeH1> >(curCoef, 1.0);
      } else {
        EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
      }
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

  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================
  void MechPDE::InitTimeStepping()
  {
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(TimeSchemeGLM::NEWMARK, 0) );
    feFunctions_[MECH_DISPLACEMENT]->SetTimeScheme(myScheme);
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
    
    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MECH_DISPLACEMENT] = "fix";
    idbcSolNameMap_[MECH_DISPLACEMENT] = "displacement";
    
    // this defines the primary unknown
    results_.Push_back( disp );
    
    // define functor for interpolation of the field
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MECH_DISPLACEMENT];
    DefineFieldResult( feFct, disp );
  }
    
  void MechPDE::DefinePostProcResults() {

    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
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
    shared_ptr<BaseFeFunction> vFct;
    if ( analysistype_ != STATIC ) {
      // === MECHANIC VELOCITY ===
      shared_ptr<ResultInfo> vel(new ResultInfo);
      vel->resultType = MECH_VELOCITY;
      vel->dofNames = dispDofNames;
      vel->unit = "m/s";
      vel->entryType = ResultInfo::VECTOR;
      vel->definedOn = ResultInfo::NODE;
      availResults_.insert( vel );
      DefineTimeDerivResult( MECH_VELOCITY, 1, MECH_DISPLACEMENT );
      vFct = timeDerivFeFunctions_[MECH_VELOCITY];
      
      // === MECHANIC ACCELERATION ===
      shared_ptr<ResultInfo> acc(new ResultInfo);
      acc->resultType = MECH_ACCELERATION;
      acc->dofNames = dispDofNames;
      acc->unit = "m/s^2";
      acc->entryType = ResultInfo::VECTOR;
      acc->definedOn = ResultInfo::NODE;
      availResults_.insert( acc );
      DefineTimeDerivResult( MECH_ACCELERATION, 2, MECH_DISPLACEMENT );
    }

    // === MECHANIC RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = MECH_RHS_LOAD;
    rhs->dofNames = dispDofNames;
    rhs->unit = "N";
    rhs->entryType = ResultInfo::VECTOR;
    rhs->definedOn = ResultInfo::NODE;
    DefineFieldResult( rhsFeFunctions_[MECH_DISPLACEMENT], rhs );
    
    // === MECHANIC STRESS ===
    shared_ptr<ResultInfo> stress(new ResultInfo);
    stress->resultType = MECH_STRESS;
    stress->dofNames = stressComponents;
    stress->unit =  "N/m^2";
    stress->entryType = ResultInfo::TENSOR;
    stress->definedOn = ResultInfo::ELEMENT;
    shared_ptr<CoefFunctionFormBased> sigmaFunc;
    if( isComplex_ ) {
      sigmaFunc.reset(new CoefFunctionFlux<Complex>(feFct, stress));
    } else {
      sigmaFunc.reset(new CoefFunctionFlux<Double>(feFct, stress));
    }
    DefineFieldResult( sigmaFunc, stress );
    
    // === MECHANIC STRAIN ===
    shared_ptr<ResultInfo> strain(new ResultInfo);
    strain->resultType = MECH_STRAIN;
    strain->dofNames = stressComponents;
    strain->unit =  "";
    strain->entryType = ResultInfo::TENSOR;
    strain->definedOn = ResultInfo::ELEMENT;
    shared_ptr<CoefFunctionFormBased> strainFunc;
    if( isComplex_ ) {
      strainFunc.reset(new CoefFunctionBOp<Complex>(feFct, strain));
    } else {
      strainFunc.reset(new CoefFunctionBOp<Double>(feFct, strain));
    }
    DefineFieldResult( strainFunc, strain );

    
    PtrCoefFct intensFct;
    shared_ptr<CoefFunctionFormBased> kedFunc;
    shared_ptr<ResultFunctor> keFunc;
    shared_ptr<CoefFunctionSurf> sNormStructIntens;
    shared_ptr<ResultFunctor> powerFunc;
    if ( analysistype_ != STATIC ) {
      // === MECHANIC STRUCTURAL INTENSTIY ===
      shared_ptr<ResultInfo> intens(new ResultInfo);
      intens->resultType = MECH_STRUCT_INTENSTIY;
      intens->dofNames = dispDofNames ;
      intens->unit =  "N/ms";
      intens->entryType = ResultInfo::VECTOR;
      intens->definedOn = ResultInfo::ELEMENT;

      // The mechanic structural intensity calculates as
      // I = -[stress] * conj(v)
      PtrCoefFct velFnc = this->GetCoefFct( MECH_VELOCITY );
      // define temporary function, without the -1 sign
      PtrCoefFct intensTmp = 
          CoefFunction::Generate(part,
                                 CoefXprBinOp( sigmaFunc, velFnc,
                                               CoefXpr::OP_MULT_VOIGT_TENSOR_VEC_CONJ ) );
      intensFct = 
          CoefFunction::Generate(part,
                                 CoefXprBinOp( "-1.0", intensTmp , CoefXpr::OP_MULT ));
      DefineFieldResult( intensFct, intens );
      
      
      // === MECHANIC KINETIC ENERGY DENSITY ===
      shared_ptr<ResultInfo> kinEnergyDens(new ResultInfo);
      kinEnergyDens->resultType = MECH_KIN_ENERGY_DENS;
      kinEnergyDens->dofNames = "";
      kinEnergyDens->unit = "Ws/m^3";
      kinEnergyDens->entryType = ResultInfo::SCALAR;
      kinEnergyDens->definedOn = ResultInfo::ELEMENT;
      if( isComplex_ ) {
        kedFunc.reset(new CoefFunctionBdBKernel<Complex>(vFct, 0.5));
      } else {
        kedFunc.reset(new CoefFunctionBdBKernel<Double>(vFct, 0.5));
      }
      DefineFieldResult( kedFunc, kinEnergyDens );
      
      // === MECHANIC KINETIC ENERGY ===
      shared_ptr<ResultInfo> kinEnergy(new ResultInfo);
      kinEnergy->resultType = MECH_KIN_ENERGY;
      kinEnergy->dofNames = "";
      kinEnergy->unit = "Ws";
      kinEnergy->entryType = ResultInfo::SCALAR;
      kinEnergy->definedOn = ResultInfo::REGION;
      availResults_.insert ( kinEnergy );
      if( isComplex_ ) {
        keFunc.reset(new EnergyResultFunctor<Complex>(vFct, kinEnergy));
      } else {
        keFunc.reset(new EnergyResultFunctor<Double>(vFct, kinEnergy));
      }
      
      resultFunctors_[MECH_KIN_ENERGY] = keFunc;
      
      // === MECHANIC POWER ===
      shared_ptr<ResultInfo> power(new ResultInfo);
      power->resultType = MECH_POWER;
      power->dofNames = "";
      power->unit = "W";
      power->entryType = ResultInfo::SCALAR;
      power->definedOn = ResultInfo::SURF_REGION;
      availResults_.insert( power );
      sNormStructIntens.reset(new CoefFunctionSurf(true));
      // then, integrate values
      if( isComplex_ ) {
        powerFunc.reset(new ResultFunctorIntegrate<Complex>(sNormStructIntens, 
                                                            feFct, power ) );
      } else {
        powerFunc.reset(new ResultFunctorIntegrate<Double>(sNormStructIntens, 
                                                           feFct, power ) );
      }
      resultFunctors_[MECH_POWER] = powerFunc;
    }


    // === MECHANIC DEFORMATION ENERGY DENSITY ===
    shared_ptr<ResultInfo> defEnergyDens(new ResultInfo);
    defEnergyDens->resultType = MECH_DEFORM_ENERGY_DENS;
    defEnergyDens->dofNames = "";
    defEnergyDens->unit = "Ws/m^3";
    defEnergyDens->entryType = ResultInfo::SCALAR;
    defEnergyDens->definedOn = ResultInfo::ELEMENT;
    shared_ptr<CoefFunctionFormBased> dedFunc;
    if( isComplex_ ) {
      dedFunc.reset(new CoefFunctionBdBKernel<Complex>(feFct, 0.5));
    } else {
      dedFunc.reset(new CoefFunctionBdBKernel<Double>(feFct, 0.5));
    }
    DefineFieldResult( dedFunc, defEnergyDens );

    // === MECHANIC DEFORMATION ENERGY ===
    shared_ptr<ResultInfo> defEnergy(new ResultInfo);
    defEnergy->resultType = MECH_DEFORM_ENERGY;
    defEnergy->dofNames = "";
    defEnergy->unit = "Ws";
    defEnergy->entryType = ResultInfo::SCALAR;
    defEnergy->definedOn = ResultInfo::REGION;
    availResults_.insert( defEnergy );
    shared_ptr<ResultFunctor> deFunc;
    if( isComplex_ ) {
      deFunc.reset(new EnergyResultFunctor<Complex>(feFct, defEnergy));
    } else {
      deFunc.reset(new EnergyResultFunctor<Double>(feFct, defEnergy));
    }
    resultFunctors_[MECH_DEFORM_ENERGY] = deFunc;

    // === MECHANIC TOTAL ENERGY DENSITY ===
    shared_ptr<ResultInfo> totEnergyDens(new ResultInfo);
    totEnergyDens->resultType = MECH_TOTAL_ENERGY_DENS;
    totEnergyDens->dofNames = "";
    totEnergyDens->unit = "Ws";
    totEnergyDens->entryType = ResultInfo::SCALAR;
    totEnergyDens->definedOn = ResultInfo::ELEMENT;
    shared_ptr<CoefFunction> tedFunc;
    // in static analysis, the total energy density equals the deformation one
    if (analysistype_ == STATIC ) {
      tedFunc = dedFunc;
    } else {
      tedFunc = CoefFunction::Generate(part, 
                                       CoefXprBinOp( dedFunc, kedFunc, CoefXpr::OP_ADD) );
    }
    DefineFieldResult(tedFunc, totEnergyDens ); 

    // === MECHANIC TOTALENERGY ===
    shared_ptr<ResultInfo> tEnergy(new ResultInfo);
    tEnergy->resultType = MECH_TOTAL_ENERGY;
    tEnergy->dofNames = "";
    tEnergy->unit = "Ws";
    tEnergy->entryType = ResultInfo::SCALAR;
    tEnergy->definedOn = ResultInfo::REGION;
    availResults_.insert( tEnergy );
    shared_ptr<ResultFunctor> teFunc;
    if( isComplex_ ) {
      teFunc.reset(new ResultFunctorIntegrate<Complex>(tedFunc, feFct, tEnergy ) );
    } else {
      teFunc.reset(new ResultFunctorIntegrate<Double>(tedFunc, feFct, tEnergy ) );
    }
    resultFunctors_[MECH_TOTAL_ENERGY] = teFunc;


    // === MECHANIC DISPLACED SURFACE VOLUME ===
     // ... to be implemented
    
    // ============================
    // Initialize result functors:
    // ============================
     
    // 1) Loop over all BDB-integrators
    std::map<RegionIdType, BaseBDBInt*>::iterator it = bdbInts_.begin();
    for( ; it != bdbInts_.end(); ++it ) {
      RegionIdType region = it->first;
      BaseBDBInt* bdb = it->second;

      sigmaFunc->AddIntegrator(bdb, region);
      strainFunc->AddIntegrator(bdb, region);
      dedFunc->AddIntegrator(bdb, region);
      deFunc->AddIntegrator(bdb, region);
      if( analysistype_ != STATIC ) {
        sNormStructIntens->SetVolumeCoef( region, intensFct );
      }
    }
    
    // 2) Loop over all mass integrators
    if( analysistype_ != STATIC ) {
      it = massInts_.begin();
      for( ; it != massInts_.end(); ++it ) {
        RegionIdType region = it->first;
        BaseBDBInt* massInt = it->second;
        kedFunc->AddIntegrator(massInt, region);
        keFunc->AddIntegrator(massInt, region);
      }
    }
  }
  
  
  std::map<SolutionType, shared_ptr<FeSpace> >
   MechPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
     
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    
    if( formulation == "default" || formulation == "H1" ){
      PtrParamNode potSpaceNode = infoNode->Get("mechDisplacement");
      crSpaces[MECH_DISPLACEMENT] =
          FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[MECH_DISPLACEMENT]->Init(solStrat_);

    }else{
       EXCEPTION( "The formulation " << formulation 
                  << "of the mechanic PDE is not known!" );
     }
     return crSpaces;
   }

} // end namespace CoupledField
