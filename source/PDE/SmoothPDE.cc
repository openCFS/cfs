// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "SmoothPDE.hh"

#include <sstream>
#include <iomanip>
#include <set>

#include "General/defs.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Driver/Assemble.hh"

#include "DataInOut/SimState.hh"
#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"

// include elements
#include "FeBasis/FeFunctions.hh"
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormalTrans.hh"
#include "Forms/Operators/SurfaceNormalStressOperator.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionContactForceDensity.hh"
#include "Domain/CoefFunction/CoefFunctionGeom.hh"

#include "Domain/CoefFunction/CoefXpr.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

namespace CoupledField {

  DEFINE_LOG(smoothpde, "smoothpde")
          
  
  SmoothPDE::SmoothPDE(Grid * aptgrid, PtrParamNode paramNode,PtrParamNode infoNode,
          shared_ptr<SimState> simState, Domain* domain )
  :SinglePDE( aptgrid, paramNode, infoNode, simState, domain ) {
    pdename_          = "smooth";
    pdematerialclass_ = SMOOTH;
    
    nonLin_        = false;
    nonLinMaterial_= false;
    
    //! Always use total Lagrangian formulation 
    updatedGeo_        = false;
    
    // ****************************
    // DETERMINE GEOMETRY
    // ****************************
    
    // Get problem geometry and PDE subtype
    myParam_->GetValue("subType", subType_ );
    
    std::string probGeo = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>();
    
    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry
    if ( subType_ == "3d" && probGeo == "3d" ) {
      stressDim_ = 6;
      tensorType_ = FULL;
      dofNames_ = "x", "y", "z";
    }
    else if (subType_ == "2.5d" && probGeo == "plane")
    {
      stressDim_ = 6;
      tensorType_ = FULL;
      dofNames_ = "x", "y", "z";
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = true;
      stressDim_ = 4;
      tensorType_ = AXI;
      dofNames_ = "r", "z";
    }
    else if ( subType_ == "planeStrain" && probGeo == "plane" ) {
      stressDim_ = 3;
      tensorType_ = PLANE_STRAIN;
      dofNames_ = "x", "y";
    }
    
    else if ( subType_ == "planeStress" && probGeo == "plane" ) {
      stressDim_ = 3;
      tensorType_ = PLANE_STRESS;
      dofNames_ = "x", "y";
    }
    
    else if ( subType_ == "flatShell" ) {
      stressDim_ = 3;
      dofNames_ = "x", "y";
    }
    else
      EXCEPTION("Subtype '" <<  subType_ << "' of PDE '" <<  pdename_ <<  "' does not fit to problem  geometry '" << probGeo << "'");
    
    // Sanity check: 3D can only be computed if 3D elements are present/
    if(subType_ == "3d" && ptGrid_->GetNumElemOfDim(3) == 0)
      EXCEPTION("Can not calculate 3D smooth without 3D elements in the grid!");
    
  }
  
  
  SmoothPDE::~SmoothPDE()
  {
    
  }
  
  void SmoothPDE::InitNonLin()
  {
    SinglePDE::InitNonLin();
  }
  
  
  void SmoothPDE::DefineIntegrators() {
    
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    // Get FESpace and FeFunction of smooth displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[SMOOTH_DISPLACEMENT];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    
    unsigned int rows = 0;
    unsigned int cols = 0;
    materials_.begin()->second->GetTensorCoefFnc(SMOOTH_STIFFNESS_TENSOR, tensorType_, Global::REAL)->GetTensorSize(rows, cols);
    assert(rows > 0 && cols > 0);
    
    //  Loop over all regions
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion++){
      actRegion = regions_[iRegion];
      actSDMat    = materials_[actRegion];
      
      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
      
      // complex material or bloch mode with complex B-matrices
      bool isComplex = complexMatData_[actRegion];

      // get list of nonlinearities
      //StdVector<NonLinType> & nonLinTypes = regionNonLinTypes_[actRegion];
      
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
      //  Standard Linear Stiffness (with potential non-linear scaling)
      // ====================================================================
      
        
      PtrCoefFct factor = CoefFunction::Generate( mp_, Global::REAL, "1.0");
      
      BaseBDBInt* stiffInt;
      
      stiffInt =  GetStiffIntegrator(actSDMat, actRegion, myFct, isComplex);
      
      stiffInt->SetName("LinElastInt");
      stiffInt->SetFeSpace( mySpace);
      
      BiLinFormContext * stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS );
      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions( myFct, myFct );
      
      
      assemble_->AddBiLinearForm( stiffIntDescr );
      
      // Important: Add bdb-integrator to global list, as we need them later
      // for calculation of postprocessing results
      bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) );
      LOG_TRACE(smoothpde) << "Add Lin BDB" << std::endl;
      
      // write to info-xml
      PtrParamNode form = infoNode_->Get("header")->Get("integrators")->Get("matrixBiLinearForms")->GetByVal("bilinearForm","integrator","LinElastInt",ParamNode::APPEND);
      PtrParamNode coef = form->Get("coef", ParamNode::APPEND);
      coef->Get("value")->SetValue(stiffInt->GetCoef()->ToString());
      
      
      // in the end, at the region to the feFunction      // to be implemented
      myFct->AddEntityList( actSDList );
    }
  }

  
  void SmoothPDE::DefineSurfaceIntegrators( ){

    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    if( bcNode ) {

      //========================================================================================
      // Normal X : assumes a normal traction proportional to the normal Y
      // X/Y = Stiffness/Displacement or Damping/Velocity or Mass/Acceleration
      // Essentially the smooth boundary traction arising in the weak form u'*t_n
      // is replaced by u'*k u_n = u'*(k u*n n) = k u'*n u*n
      // where k is the parameter read from the input file, and u represents the unknown Y
      //========================================================================================
      typedef std::pair<std::string, FEMatrixType> normalBCtype;
      std::vector<normalBCtype> normalBCs;
      normalBCs.push_back(std::make_pair("normalStiffness", STIFFNESS));
      normalBCs.push_back(std::make_pair("normalDamping", DAMPING));
      normalBCs.push_back(std::make_pair("normalMass", MASS));
      for( std::vector<normalBCtype>::iterator it = normalBCs.begin() ; it != normalBCs.end(); ++it ){
        std::string xmlName = it->first;
        FEMatrixType feMat = it->second;
        LOG_DBG(smoothpde) << "Reading '" << xmlName << "' definition";
        StdVector<shared_ptr<EntityList> > ent;
        StdVector<PtrCoefFct > kCoef;
        StdVector<std::string> volumeRegions;
        ReadRhsExcitation( xmlName , feFunctions_[SMOOTH_DISPLACEMENT]->GetResultInfo()->dofNames, ResultInfo::SCALAR, ent, kCoef, updatedGeo_, volumeRegions);
        for( UInt i = 0; i < ent.GetSize(); ++i ) {
          // get the volume region for defining the correct normal direction
          RegionIdType aRegion = ptGrid_->GetRegion().Parse(volumeRegions[i]);
          std::set<RegionIdType> volRegion;
          volRegion.insert(aRegion);
          // check type of entitylist
          if (ent[i]->GetType() == EntityList::NODE_LIST) {
            EXCEPTION( xmlName << " must be defined on (surface) elements")
          }
          if ( kCoef[i]->IsComplex() && !(isComplex_) ) {
            EXCEPTION( xmlName << " is defied as complex but PDE is not")
          }
          // setup the integrator for: u'*t_n = u'*k u_n = u'*(k u*n n) = k u'*n u*n
          BiLinearForm * tangInt = NULL;
          if( kCoef[i]->IsComplex() ) {
            if (dim_ == 2){
              tangInt = new SurfaceBBInt<Complex,Double>(new IdentityOperatorNormalTrans<FeH1,2,2>(), kCoef[i], Complex(1.0,0), volRegion, updatedGeo_ );
            } else {
              tangInt = new SurfaceBBInt<Complex,Double>(new IdentityOperatorNormalTrans<FeH1,3,3>(), kCoef[i], Complex(1.0,0), volRegion, updatedGeo_ );
            }
          } else {
            if (dim_ == 2){
              tangInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,2,2>(), kCoef[i], 1.0, volRegion, updatedGeo_ );
            } else {
              tangInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,3,3>(), kCoef[i], 1.0, volRegion, updatedGeo_ );
            }
          }
          tangInt->SetName(xmlName + "Integrator");
          BiLinFormContext *tangContext = new BiLinFormContext(tangInt, feMat );
          tangContext->SetEntities( ent[i], ent[i]);
          tangContext->SetFeFunctions( feFunctions_[SMOOTH_DISPLACEMENT], feFunctions_[SMOOTH_DISPLACEMENT]);
          feFunctions_[SMOOTH_DISPLACEMENT]->AddEntityList( ent[i] );
          assemble_->AddBiLinearForm( tangContext );
        }
      }
    }
  }
  
  
  
  
  void SmoothPDE::DefineNcIntegrators() {
    EXCEPTION("NC integrators not defined for smooth PDE!");
  }
  
  void SmoothPDE::DefineRhsLoadIntegrators(PtrParamNode input)
  {
    LOG_TRACE(smoothpde) << "Defining rhs load integrators for smooth PDE";
    
    // Get FESpace and FeFunction of smooth displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[SMOOTH_DISPLACEMENT];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    LinearForm * lin = NULL;
    StdVector<std::string> dispDofNames = myFct->GetResultInfo()->dofNames;
    
    /* Flag, if coefficient function lives on updated geometry (updated 
     * Lagrangian formulation). For analytically prescribed values
     * (pressure, force, stress), this is in general not the case. 
     * In the coupled case, the magnetic force density however could live
     * on an updated geometry. In this case, the ReadRhsExcitation()
     * method will set the coefUpdateGeo flag to true, so the RHS-
     * integrator also has to use this value.
     * 
     */
    bool coefUpdateGeo = false;
    
    // ========================
    //  FORCES (volume, nodal)
    // ========================
    LOG_DBG(smoothpde) << "Reading forces";
    ReadRhsExcitation("force", dispDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo, input);
    
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      
      // In case of a total force, we can not have a spatial dependency
      if(coef[i]->IsSpacialDependent()) {
        EXCEPTION("Check changes in smoothPDE for spatial dependent forces");
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
        if(numNodes > 1  && coef[i]->DoNormalize()) {
          Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;  
          coef[i] = CoefFunction::Generate(mp_, part, CoefXprVecScalOp(mp_, coef[i],
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
    //  PRESSURE 
    // ===============
    LOG_DBG(smoothpde) << "Reading smooth pressure";
    StdVector<std::string> empty;
    ReadRhsExcitation("pressure", empty, ResultInfo::SCALAR, isComplex_, ent, coef, coefUpdateGeo, input);
    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );
    
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Smooth pressure must be defined on elements")
      }
      
      // Factor for the pressure:
      // The pressure is by definition in the opposite direction as the 
      // normal stress, i.e. a positive pressure means a compressive stress
      // (<0). Thus we have to take the minus sign into account
      const Double presFac = -1.0;
      
      if( dim_ == 2) {
        if(isComplex_) {
          lin = new BUIntegrator<Complex, true> ( new IdentityOperatorNormalTrans<FeH1,2>(),
                  Complex(presFac), coef[i],
                  volRegions, coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double,true> ( new IdentityOperatorNormalTrans<FeH1,2>(),
                  presFac, coef[i], volRegions,
                  coefUpdateGeo);
        }
      } else  {
        if(isComplex_) {
          lin = new BUIntegrator<Complex, true> ( new IdentityOperatorNormalTrans<FeH1,3>(),
                  Complex(presFac), coef[i],
                  volRegions, coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double, true> ( new IdentityOperatorNormalTrans<FeH1,3>(),
                  presFac, coef[i],
                  volRegions, coefUpdateGeo);
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
    LOG_DBG(smoothpde) << "Reading surface tractions";
    
    ReadRhsExcitation("traction", dispDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo, input);
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
          lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,2,2>(),
                  Complex(1.0), coef[i], coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,2,2>(),
                  1.0, coef[i],coefUpdateGeo);
        }
      } else  {
        if(isComplex_) {
          lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,3,3>(),
                  Complex(1.0), coef[i], coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,3,3>(),
                  1.0, coef[i], coefUpdateGeo);
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
  
  BaseBDBInt* SmoothPDE::GetStiffIntegrator(BaseMaterial* actSDMat, RegionIdType regionId, shared_ptr<BaseFeFunction> feFct, bool isComplex)
  {
    // Get region name
    std::string regionName = ptGrid_->GetRegion().ToString( regionId );

    // get list of nonlinearities
    StdVector<NonLinType> & nonLinTypes = regionNonLinTypes_[regionId];

    PtrCoefFct preFac = CoefFunction::Generate( mp_, Global::REAL, "1.0");;

    // for non-lin scaling factor used for strain stiffening (modified preFac from above)
    if ( nonLinTypes.Find(STRAIN_STIFFENING) != -1 ) {
      StdVector<std::string> stressComponents;
      if( subType_ == "3d" || subType_ == "2.5d") {
        stressComponents = "xx", "yy", "zz", "yz", "xz", "xy";
      } else if( subType_ == "planeStrain" ) {
        stressComponents = "xx", "yy", "xy";
      } else if( subType_ == "planeStress" ) {
        stressComponents = "xx", "yy", "xy";
      } else if( subType_ == "axi" ) {
        stressComponents = "rr", "zz", "rz", "phiphi";
      }
      shared_ptr<CoefFunctionFormBased> energyFunc;
      shared_ptr<ResultInfo> strain(new ResultInfo);
      strain->resultType = MECH_STRAIN;
      strain->dofNames = stressComponents;
      strain->unit =  "";
      strain->entryType = ResultInfo::TENSOR;
      strain->definedOn = ResultInfo::ELEMENT;
      shared_ptr<CoefFunctionFormBased> strainFunc;
      // scaling fac to steer the influence of the strain stiffening
      PtrCoefFct scalingFac = CoefFunction::Generate( mp_, Global::REAL, std::to_string(nonLinScalarFac_));

      if( isComplex_ ) {
        strainFunc.reset(new CoefFunctionBOp<Complex>(feFct, strain));
        PtrCoefFct squaredStrain = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, strainFunc, strainFunc, CoefXpr::OP_MULT));
        PtrCoefFct rootFunc = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprUnaryOp(mp_, squaredStrain, CoefXpr::OP_SQRT));
        PtrCoefFct scaledRootFunc = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, scalingFac, rootFunc, CoefXpr::OP_MULT));
        preFac = CoefFunction::Generate(mp_, Global::COMPLEX, 
          CoefXprBinOp(mp_, CoefFunction::Generate( mp_, Global::REAL, "1.0"), scaledRootFunc, CoefXpr::OP_ADD));
      } else {
        strainFunc.reset(new CoefFunctionBOp<Double>(feFct, strain));
        PtrCoefFct squaredStrain = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, strainFunc, strainFunc, CoefXpr::OP_MULT));
        PtrCoefFct rootFunc = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, squaredStrain, CoefXpr::OP_SQRT));
        PtrCoefFct scaledRootFunc = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, scalingFac, rootFunc, CoefXpr::OP_MULT));
        preFac = CoefFunction::Generate(mp_, Global::REAL, 
          CoefXprBinOp(mp_, CoefFunction::Generate( mp_, Global::REAL, "1.0"), scaledRootFunc, CoefXpr::OP_ADD));
      }
      stiffFormCoefs_.insert(strainFunc);
    } // strain stiffening end
    
    // ------------------------
    //  Obtain linear material
    // ------------------------
    shared_ptr<CoefFunction > curCoef;
    PtrCoefFct stiffCoef;
    if ( isComplex ) {
      curCoef = actSDMat->GetTensorCoefFnc(SMOOTH_STIFFNESS_TENSOR, tensorType_, Global::COMPLEX);
      stiffCoef = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprTensScalOp(mp_, curCoef, preFac, CoefXpr::OP_MULT));
    } else {
      curCoef = actSDMat->GetTensorCoefFnc(SMOOTH_STIFFNESS_TENSOR, tensorType_, Global::REAL);
      stiffCoef = CoefFunction::Generate(mp_, Global::REAL, CoefXprTensScalOp(mp_, curCoef, preFac, CoefXpr::OP_MULT));
    }
    
    
    
    // store coefficient function for later use (e.g. in boundary integrators)
    regionStiffness_[regionId] = stiffCoef;
    
    // ----------------------------------------
    //  Determine correct stiffness integrator 
    // ----------------------------------------
    
    BaseBDBInt* integ = NULL;
    BaseBOperator* bOp = GetStrainOperator(isComplex, false);
    
    // ====================
    //  Standard Stiffness
    // ====================
    if (isComplex )
      integ = new BDBInt<Complex>(bOp, stiffCoef, 1.0, updatedGeo_);
    else
      integ = new BDBInt<Double>(bOp, stiffCoef, 1.0, updatedGeo_);
    
    return integ;
  }
  
  BaseBOperator* SmoothPDE::GetStrainOperator(bool isComplex, bool icModes)
  {
    BaseBOperator* bOp = NULL;
    // determine if we do bloch eigenfrequency analysis

    if(isComplex)
    {
      if( subType_ == "planeStrain" )
      {
        bOp = new StrainOperator2D<FeH1,Complex>(icModes);
      }
      if( subType_ == "axi" )
        bOp = new StrainOperatorAxi<FeH1,Complex>(icModes);
      if(subType_ == "planeStress")
        bOp = new StrainOperator2D<FeH1,Complex>(icModes);
      if(subType_ == "3d")
      {
        bOp = new StrainOperator3D<FeH1,Complex>(icModes);
      }
      if(subType_ == "2.5d")
        bOp = new StrainOperator2p5D<FeH1,Complex>(icModes);
    }
    else // not complex
    {
      if(subType_ == "axi")
        bOp = new StrainOperatorAxi<FeH1,Double>(icModes);
      if(subType_ == "planeStrain")
        bOp = new StrainOperator2D<FeH1,Double>(icModes);
      if(subType_ == "planeStress")
        bOp = new StrainOperator2D<FeH1,Double>(icModes);
      if(subType_ == "3d")
        bOp = new StrainOperator3D<FeH1,Double>(icModes);
      if(subType_ == "2.5d")
        bOp = new StrainOperator2p5D<FeH1,Double>(icModes);
    }

    if(bOp == NULL)
      EXCEPTION("strain operator not implemented for analysis type");

    return bOp;
  }

  void SmoothPDE::DefineSolveStep()
  {
    solveStep_ = new StdSolveStep(*this);
  }
  
  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================
  void SmoothPDE::InitTimeStepping()  
  {
    // we do not require a special time stepping scheme for the non-linear strain stiffening, so we set it to NONE all the time
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(GLMScheme::BDF2, TimeSchemeGLM::NONE) );
    
    feFunctions_[SMOOTH_DISPLACEMENT]->SetTimeScheme(myScheme);
  }
  
  void SmoothPDE::DefinePrimaryResults()
  {
    // Check for subType
    StdVector<std::string> stressDofNames;
    
    if(subType_ == "3d" || subType_ == "2.5d")
      stressDofNames = "xx", "yy", "zz", "yz", "xz", "xy";
    else if(subType_ == "planeStrain")
      stressDofNames = "xx", "yy", "xy";
    else if(subType_ == "planeStress")
      stressDofNames = "xx", "yy", "xy";
    else if(subType_ == "axi")
      stressDofNames = "rr", "zz", "rz", "phiphi";
    else if(subType_ == "2.5d")
      stressDofNames = "xx", "yy", "zz", "yz", "xz", "xy";
    
    // === SMOOTH DISPLACEMENT ===
    shared_ptr<ResultInfo> disp(new ResultInfo);
    disp->resultType = SMOOTH_DISPLACEMENT;
    disp->dofNames = dofNames_;
    disp->unit = "m";
    disp->entryType = ResultInfo::VECTOR;
    disp->SetFeFunction(feFunctions_[SMOOTH_DISPLACEMENT]);
    disp->definedOn = ResultInfo::NODE;
    feFunctions_[SMOOTH_DISPLACEMENT]->SetResultInfo(disp);

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[SMOOTH_DISPLACEMENT] = "fix";
    idbcSolNameMap_[SMOOTH_DISPLACEMENT] = "displacement";
    idbcSolNameMapD1_[SMOOTH_DISPLACEMENT] = "velocity";
    idbcSolNameMapD2_[SMOOTH_DISPLACEMENT] = "acceleration";
    
    // this defines the primary unknown
    results_.Push_back( disp );
    
    // define functor for interpolation of the field
    shared_ptr<BaseFeFunction> feFct = feFunctions_[SMOOTH_DISPLACEMENT];
    DefineFieldResult( feFct, disp );
  }
  
  void SmoothPDE::DefinePostProcResults() {
    StdVector<std::string> stressComponents;
    if( subType_ == "3d" || subType_ == "2.5d") {
      stressComponents = "xx", "yy", "zz", "yz", "xz", "xy";
    } else if( subType_ == "planeStrain" ) {
      stressComponents = "xx", "yy", "xy";
    } else if( subType_ == "planeStress" ) {
      stressComponents = "xx", "yy", "xy";
    } else if( subType_ == "axi" ) {
      stressComponents = "rr", "zz", "rz", "phiphi";
    }
    StdVector<std::string > dispDofNames;
    dispDofNames = feFunctions_[SMOOTH_DISPLACEMENT]->GetResultInfo()->dofNames;
    shared_ptr<BaseFeFunction> feFct = feFunctions_[SMOOTH_DISPLACEMENT];
    shared_ptr<BaseFeFunction> vFct;
    shared_ptr<BaseFeFunction> aFct;
    if ( analysistype_ != STATIC && analysistype_ != BUCKLING) {
      // === GRID VELOCITY ===
      shared_ptr<ResultInfo> vel(new ResultInfo);
      vel->resultType = SMOOTH_VELOCITY;
      vel->dofNames = dispDofNames;
      vel->unit = "m/s";
      vel->entryType = ResultInfo::VECTOR;
      vel->definedOn = ResultInfo::NODE;
      availResults_.insert( vel );
      DefineTimeDerivResult( SMOOTH_VELOCITY, 1, SMOOTH_DISPLACEMENT );
      vFct = timeDerivFeFunctions_[SMOOTH_VELOCITY];
      //feFunctions_[SMOOTH_VELOCITY] = vFct;

      // === GRID ACCELERATION ===

      // shared_ptr<ResultInfo> acc(new ResultInfo);
      // acc->resultType = SMOOTH_ACCELERATION;
      // acc->dofNames = dispDofNames;
      // acc->unit = "m^2/s";
      // acc->entryType = ResultInfo::VECTOR;
      // acc->definedOn = ResultInfo::NODE;
      // availResults_.insert( acc );
      // DefineTimeDerivResult( SMOOTH_ACCELERATION, 2, SMOOTH_DISPLACEMENT );

      // the following two lines cause problems - include whenever needed
      // aFct = timeDerivFeFunctions_[SMOOTH_ACCELERATION];
      // feFunctions_[SMOOTH_ACCELERATION] = aFct;
    }
    

    // === SMOOTH ZERO STRESS ===
    // This is a dummy result in order to get the iterative coupling to recognize
    // both PDEs as iterative coupled although it's just a forward coupling
    shared_ptr<ResultInfo> stressZero(new ResultInfo);
    stressZero->resultType = SMOOTH_ZERO_PRESSURE;
    stressZero->dofNames = "";
    stressZero->unit = MapSolTypeToUnit(SMOOTH_ZERO_PRESSURE);
    stressZero->entryType = ResultInfo::SCALAR;
    stressZero->definedOn = ResultInfo::NODE;
    availResults_.insert( stressZero );
    PtrCoefFct constZero = CoefFunction::Generate( mp_, Global::REAL, "0.0");

    DefineFieldResult( constZero, stressZero );

    // === SMOOTH STRAIN ===
    shared_ptr<ResultInfo> strain(new ResultInfo);
    strain->resultType = SMOOTH_STRAIN;
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
    stiffFormCoefs_.insert(strainFunc);

    // === SMOOTH CONTACT FORCE DENSITY ===
    shared_ptr<ResultInfo> contactForceDensity(new ResultInfo);
    contactForceDensity->resultType = SMOOTH_CONTACT_FORCE_DENSITY;
    contactForceDensity->dofNames = dispDofNames;
    contactForceDensity->unit = MapSolTypeToUnit(SMOOTH_CONTACT_FORCE_DENSITY);
    contactForceDensity->entryType = ResultInfo::VECTOR;
    contactForceDensity->definedOn = ResultInfo::SURF_ELEM;
    availResults_.insert( contactForceDensity );

    StdVector<std::string> surfList1;
    StdVector<std::string> surfList2;
    StdVector<std::string> volumeList;
    StdVector<std::string> contactLawList;
    StdVector<bool> useSurfaceMidpointsList;
    ReadContact(surfList1, surfList2, volumeList, contactLawList, useSurfaceMidpointsList);

    shared_ptr<CoefFunctionContactForceDensity> contactForceDensityFunc;
    contactForceDensityFunc.reset(new CoefFunctionContactForceDensity(feFct, surfList1, surfList2, volumeList, contactLawList, useSurfaceMidpointsList));
    
    DefineFieldResult( contactForceDensityFunc, contactForceDensity );


    // === SMOOTH CONTACT FORCE (= integral of contact force density) ===
    shared_ptr<ResultInfo> contactForce;
    contactForce.reset(new ResultInfo);
    contactForce->resultType = SMOOTH_CONTACT_FORCE;
    contactForce->dofNames = dispDofNames;
    contactForce->unit = MapSolTypeToUnit(SMOOTH_CONTACT_FORCE);
    contactForce->entryType = ResultInfo::VECTOR;
    contactForce->definedOn = ResultInfo::SURF_REGION;
    // Integrate surface traction
    shared_ptr<ResultFunctor> contactForceFct;
    if(isComplex_)
        contactForceFct.reset(new ResultFunctorIntegrate<Complex>(contactForceDensityFunc, feFct, contactForce));
    else
        contactForceFct.reset(new ResultFunctorIntegrate<Double>(contactForceDensityFunc, feFct, contactForce));
    resultFunctors_[SMOOTH_CONTACT_FORCE] = contactForceFct;
    availResults_.insert(contactForce);
  
    
    // === SMOOTH DEFORMATION ENERGY DENSITY ===
    shared_ptr<ResultInfo> defEnergyDens(new ResultInfo);
    defEnergyDens->resultType = SMOOTH_DEFORM_ENERGY_DENS;
    defEnergyDens->dofNames = "";
    defEnergyDens->unit = MapSolTypeToUnit(SMOOTH_DEFORM_ENERGY_DENS);
    defEnergyDens->entryType = ResultInfo::SCALAR;
    defEnergyDens->definedOn = ResultInfo::ELEMENT;
    shared_ptr<CoefFunctionFormBased> dedFunc;
    if( isComplex_ ) {
      dedFunc.reset(new CoefFunctionBdBKernel<Complex>(feFct, 0.5));
    } else {
      dedFunc.reset(new CoefFunctionBdBKernel<Double>(feFct, 0.5));
    }
    DefineFieldResult( dedFunc, defEnergyDens );
    stiffFormCoefs_.insert(dedFunc);
  

    // === SMOOTH DEFORMATION ENERGY===
    shared_ptr<ResultInfo> defEnergy(new ResultInfo);
    defEnergy->resultType = SMOOTH_DEFORM_ENERGY;
    defEnergy->dofNames = "";
    defEnergy->unit = MapSolTypeToUnit(SMOOTH_DEFORM_ENERGY);
    defEnergy->entryType = ResultInfo::SCALAR;
    defEnergy->definedOn = ResultInfo::REGION;
    availResults_.insert(defEnergy);
    shared_ptr<ResultFunctor> energyFunc;
    if( isComplex_ ) {
      energyFunc.reset(new EnergyResultFunctor<Complex>(feFct, defEnergy, 0.5));
    } else {
      energyFunc.reset(new EnergyResultFunctor<Double>(feFct, defEnergy, 0.5));
    }
    resultFunctors_[SMOOTH_DEFORM_ENERGY] = energyFunc;

    // DefineFieldResult( energyFunc, defEnergy );

    stiffFormFunctors_.insert(energyFunc);


    // === SMOOTH JACOBIAN ===
    shared_ptr<ResultInfo> jac(new ResultInfo);
    jac->resultType = SMOOTH_JACOBIAN;
    jac->dofNames = "";
    jac->unit =  "";
    jac->entryType = ResultInfo::SCALAR;
    jac->definedOn = ResultInfo::ELEMENT;
    shared_ptr<CoefFunctionGeom> coefFctJac;
    coefFctJac.reset(new CoefFunctionGeom("Jacobian"));

    DefineFieldResult( coefFctJac, jac );
  }
  void SmoothPDE::ReadContact(StdVector<std::string>& surfList1,
                                  StdVector<std::string>& surfList2,
                                  StdVector<std::string>& volumeList,
                                  StdVector<std::string>& contactLawList,
                                  StdVector<bool>& useSurfaceMidpointsList) {

    // Check if the node is defined
    PtrParamNode bcNode = myParam_->Get("contactList", ParamNode::PASS);
    if (bcNode) {
      StdVector<shared_ptr<EntityList> > ent;
      StdVector<PtrCoefFct > coef;
      std::string elemName = "contact";
      ParamNodeList elems = myParam_->Get("contactList")->GetList(elemName);
      
      ent.Resize(elems.GetSize());
      coef.Resize(elems.GetSize());
      
      // Define xml input
      std::string surfName1;
      std::string surfName2;
      std::string volumeName;
      std::string contactLaw;
      bool useSurfaceMidpoints;

      surfList1.Clear();
      surfList2.Clear();
      volumeList.Clear();
      contactLawList.Clear();
      useSurfaceMidpointsList.Clear();
      
      //std::vector<std::pair<std::string,std::string>> surfList;
      for (UInt i = 0; i < elems.GetSize(); ++i) {
        PtrParamNode xml = elems[i];
        try {
          // read xml input
          surfName1 = xml->Get("Surface1")->As<std::string>();
          surfName2 = xml->Get("Surface2")->As<std::string>();
          volumeName = xml->Get("Volume")->As<std::string>();
          contactLaw = xml->Get("contactLaw")->As<std::string>();
          useSurfaceMidpoints = xml->Get("useSurfaceMidpoints")->As<bool>();
          
          // check if we only have surface elements in the given region
          EntityList::ListType listType = EntityList::ELEM_LIST;
          if( ptGrid_->GetEntityDim(surfName1) == ptGrid_->GetDim()-1 && ptGrid_->GetEntityDim(surfName2) == ptGrid_->GetDim()-1 ) {
          listType = EntityList::SURF_ELEM_LIST;
          } else {
            throw Exception("contactForce can only be evaluated on surface regions");
          }

          // finish checks
          switch( ptGrid_->GetEntityType(surfName1) ) {
            case EntityList::NAMED_NODES:
            ent[i] = ptGrid_->GetEntityList( EntityList::NODE_LIST, surfName1);
            break;
          case EntityList::REGION:
          case EntityList::NAMED_ELEMS:
            ent[i] = ptGrid_->GetEntityList( listType, surfName1 );
            break;
          case EntityList::NO_TYPE:
            EXCEPTION("No entities with name '" << surfName1 << "' known");
            break;
          }
          // check just defined type of entitylist
          if (ent[i]->GetType() == EntityList::NODE_LIST) {
            EXCEPTION("contactForce must be defined on elements")
          }
          switch( ptGrid_->GetEntityType(surfName2) ) {
            case EntityList::NAMED_NODES:
            ent[i] = ptGrid_->GetEntityList( EntityList::NODE_LIST, surfName2);
            break;
          case EntityList::REGION:
          case EntityList::NAMED_ELEMS:
            ent[i] = ptGrid_->GetEntityList( listType, surfName2 );
            break;
          case EntityList::NO_TYPE:
            EXCEPTION("No entities with name '" << surfName2 << "' known");
            break;
          }
          // check just defined type of entitylist
          if (ent[i]->GetType() == EntityList::NODE_LIST) {
            EXCEPTION("contactForce must be defined on elements")
          }

          // add everything to the lists
          surfList1.push_back(surfName1);
          surfList2.push_back(surfName2);
          volumeList.push_back(volumeName);
          contactLawList.push_back(contactLaw);
          useSurfaceMidpointsList.push_back(useSurfaceMidpoints);


        } catch (Exception &e) {
          RETHROW_EXCEPTION(e, pdename_ << ": Could not read definition for '" << elemName
                                      << "' on entities '" << surfName1 << "' and' " << surfName2 <<"'");
        }
      }
    }
  }

  
  std::map<SolutionType, shared_ptr<FeSpace> > SmoothPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode)
  {
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    
    if( formulation == "default" || formulation == "H1" )
    {
      PtrParamNode potSpaceNode = infoNode->Get("smoothDisplacement");
      crSpaces[SMOOTH_DISPLACEMENT] = FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[SMOOTH_DISPLACEMENT]->Init(solStrat_);
    }
    else
      EXCEPTION( "The formulation " << formulation << "of the smooth PDE is not known!" );
    return crSpaces;
  }
} // end namespace CoupledField
