// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "StdPDE.hh"

#include "MatVec/vector.hh"
#include "MatVec/generatematvec.hh"
#include "Driver/stdSolveStep.hh"
#include "Driver/transientdriver.hh"

#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"

// headers for Paramhandling
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/CFSOLASParams.hh"
#include "Driver/assemble.hh"

#include "OLAS/algsys/algebraicSys.hh"

namespace CoupledField {

  StdPDE::StdPDE(Grid *aptgrid, PtrParamNode paramNode ) :
    BasePDE(paramNode),
    ptgrid_(aptgrid),
    subType_(),
    numCouplingBcs_(0),
    nonLin_(false),
    nonLinMaterial_(false),
    isHysteresis_(false),
    totalFormulation_(false),
    isIterCoupled_(false),
    updateCouplingBCs_(false),
    ptCoupling_(NULL),
    iterCoupledCounter_(0),
    TS_alg_(NULL),
    effectiveMass_(false),
    diagMass_(false),
    firstTimeStepStatic_(true),
    isInstationary_(false),
    needsAlgsys_(true),
    isAlwaysStatic_(false),
    dim_(ptgrid_->GetDim()), 
    isaxi_(param->Get("domain")->Get("geometryType")->As<std::string>() == "axi"),
    isComplex_(false),    
    needSolPrev_(false),
    fracDamping_(false),
    fracMemory_(0),
    isBiotSavart_(false),
    inType_(NOTUSED),
    isIncrFormulation_(false),
    updatedLagrangeForm_(false),
    ComputeRHSforHarm_(false),
    assemble_(NULL),
    solveStep_(NULL),
    algsys_(NULL),
    isSetInitialCondition_(false),
    InitialCondition_(0.0)
  {
  }
  
  StdPDE::~StdPDE() 
  {
    
  }

  bool StdPDE::HasPeriodicBC()
  {
    for(UInt i = 0; i < constraints_.GetSize(); i++)
      if(constraints_[i]->periodic) return true;

    return false;
  }

  // ======================================================
  // GRID SECTION (Meshing, ...) 
  // ======================================================

  void StdPDE::CreateMatrices_Solver() {
    
    
    // ==============================
    // create and initialize matrices 
    // ==============================

    // create algebraic system and intialize matrices
    algsys_->CreateLinSys();
    algsys_->InitMatrix();
    
    // Check for analysistype
    if ( analysistype_ != EIGENFREQUENCY ) {
      
      // create solver and preconditioner
      algsys_->CreateSolver();
      algsys_->CreatePrecond();
      
    } else {
      // create eigenvalue solver
      algsys_->CreateEigenSolver();
    }
        
    // now reset AlgebraicSystem 
    algsys_->InitRHS();
  	algsys_->InitSol();
  }


  BaseSolveStep * StdPDE::GetSolveStep() {
    
    return solveStep_;
  }
  

  // ======================================================
  // ALGSYS SECTION (SOLVER, ...) 
  // ======================================================

  PtrParamNode StdPDE::FindLinearSystem(const std::string& sysName) {

    PtrParamNode pn, linSysNode;
    
    pn = param->GetByVal("sequenceStep", "index", sequenceStep_, ParamNode::PASS);
    linSysNode = pn->Get("linearSystems", ParamNode::INSERT);
    pn = linSysNode->GetByVal("system", "name", sysName, ParamNode::INSERT);
    
    // If no system with the specified name could be found in XML file
    // we just generate a new ParamNode.
    //WARN("Check, if <linearSystems> node is created properly");
    //    if(!pn) {
//      
//      
//      linSysNode->GetChildren().Push_back(PtrParamNode(new ParamNode());
//      pn = linSysNode->GetChildren().Last();
//      pn->SetName("system");
//      pn->GetChildren().Push_back(new ParamNode());
//      PtrParamNode nameNode = pn->GetChildren().Last();
//      nameNode->SetName("name");
//      nameNode->SetValue(sysName);
//    }
    
    return pn;
  }

  Double StdPDE::GetFracDampMatrixCoeff(RegionIdType regionId) {
    

    Double coeff;

    // pre factor of fractional derivative (same for all algorithms)
    TransientDriver * driver = NULL; 
    driver = dynamic_cast<TransientDriver*>(domain->GetSingleDriver() );
    
    if( driver == NULL) {
      EXCEPTION( "Fractional damping only possible for transient simulation!" );
    }
    Double dt = driver->GetDeltaT();

    Double y;
    materials_[regionId]->GetScalar(y,FRACTIONAL_EXPONENT,Global::REAL);

    coeff = std::exp(-(y-1.0) * std::log(dt));

    // needed for formulation with only MASS and STIFFNESS matrix
    // pre factor of Newmark time stepping scheme
    Double beta = TS_alg_->GetNewmarkBeta();
    coeff *= 1.0 / (beta*dt*dt);

    /*
    (*debug) << std::endl << "Parameters in GetFracDampMatrixCoeff are:"
             << std::endl << "regionId: " << regionId
             << std::endl << "y    = " << y
             << std::endl << "dt   = " << dt
             << std::endl << "beta = " << beta
             << std::endl << "coeff= " << coeff << std::endl << std::endl;
    */ 
    return coeff;
  }

  DampingType StdPDE::GetDamping(RegionIdType reg_id) const
  {
    std::map<RegionIdType,DampingType>::const_iterator it = dampingList_.find(reg_id);

    return it != dampingList_.end() ? it->second : NONE;
  }


  shared_ptr<ResultInfo> StdPDE::GetResultInfo( SolutionType solType ) {
    
    shared_ptr<ResultInfo> res;
    ResultSet::const_iterator it = availResults_.begin();
    for( ; it != availResults_.end(); ++it )  {
      if( (*it)->resultType == solType) {
        res = *it;
        break;
      }
    }

    // check, if resultinfo was found
    if( !res.get() ) {
      EXCEPTION( "A result with solutionType '" << SolutionTypeEnum.ToString(solType)
                 << "' was not found for " << pdename_ );
    }

    return res;
  }

 void StdPDE::ReadDisplacementAndUpdateGrid( UInt step)
 {
   /* only update grid if langrange type has been set */
   if ( !updatedLagrangeForm_ )
   {
     return;
   }
   /* do not set new grid in the first step */
   if ( step != 0 )
   {
     for ( UInt nreg = 0; nreg < regions4GridDisplacements_.GetSize(); nreg++ )
     {
       ResultHandler* resultHandler = domain->GetResultHandler();
       shared_ptr<BaseResult> gridDisplacement = resultHandler->GetResult( fileName4GridDisplacements_,
           1,
           step,
           MECH_DISPLACEMENT,        
           regions4GridDisplacements_[nreg] );

       Result<Double> *result =
         dynamic_cast<Result<Double>*>(&(*gridDisplacement));
       if (result == NULL)
       {
         EXCEPTION("Cannot read result 'Grid-Displacements' from input id '"
             <<  fileName4GridDisplacements_ << "'");
       }
       Vector<Double>& resVec = result->GetVector();
       shared_ptr<EntityList> nodesList = gridDisplacement->GetEntityList();
       StdVector<UInt> nodes;

       EntityIterator it;

       it = nodesList->GetIterator();
       for( it.Begin(); !it.IsEnd(); it++ )
       {
         nodes.Push_back(it.GetNode());
       }

       ptgrid_->SetNodeOffset(nodes, resVec);
     }
   }
 }


  // ******************
  //   ReadOlasParams
  // ******************
  void StdPDE::ReadOlasParams( std::string sysName ) {


    // Log to .las file
    (*cla) <<  " --- CFS: Setting parameters for linear system '"
           << sysName << " ---" << std::endl;

    // Set parameters for OLAS
    std::string amExpert = "no";
    param->GetValue( "override", amExpert, ParamNode::PASS );

    PtrParamNode linSysNode;
    PtrParamNode temp = param->GetByVal("sequenceStep", "index", sequenceStep_);
    temp = temp->Get("linearSystems", ParamNode::INSERT);
    linSysNode = temp ->GetByVal("system", "name", sysName, ParamNode::INSERT );
    CFSOLASParams::SetParams( sysName, linSysNode, 
                              analysistype_, assemble_,
                              (amExpert == "yes") );

  }

  //============================================================================================
  //FeFunction Methods
  //============================================================================================

  shared_ptr<BaseFeFunction> StdPDE::GetFeFunction( SolutionType solType ) {
    
    //TODO> We need to find a more failsafe way to store the entity names associated with a 
    //FeFunction
    shared_ptr<BaseFeFunction> feFct;
    SolutionType mySolType;
    if( feFunctions_.find(solType) == feFunctions_.end()){
      //ok so it could be that we are looking for a postProc Result
      if(postProcResults_.find(solType) != postProcResults_.end()){
        mySolType = postProcResults_[solType];
      }else{
          EXCEPTION( "A FeFunction descriptor with solutionType '" << SolutionTypeEnum.ToString(solType)
                      << "' was not found for " << pdename_ );
      }
    }else{
      mySolType = solType;
    }
    if(feFunctions_.find(mySolType) != feFunctions_.end()){
      feFct = feFunctions_[mySolType];
    }else{
      EXCEPTION("StdPDE::GetFeFunction: Could not find the corresponding FeFunction for the given entity name\n \
                         Did you specify all Regions, Surfregions and NamedNodes in the xml?");
    }
    return feFct;

  }

} // end of namespace
