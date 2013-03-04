// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <math.h>
#include <stdlib.h>
#include <complex>
#include <fstream>
#include <set>
#include <utility>

#include "DataInOut/ParamHandling/CFSOLASParams.hh"
// headers for Paramhandling
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/resultHandler.hh"
#include "Domain/domain.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "Driver/singleDriver.hh"
#include "Driver/stdSolveStep.hh"
#include "Driver/transientdriver.hh"
#include "General/Enum.hh"
#include "MatVec/vector.hh"
#include "Materials/baseMaterial.hh"
#include "OLAS/algsys/basesystem.hh"
#include "PDE/eqnMap.hh"
#include "StdPDE.hh"
#include "Utils/basenodestoresol.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/result.hh"
#include "Utils/tools.hh"

namespace CoupledField {


class BaseSolveStep;
class SingleVector;
struct Elem;

  StdPDE::StdPDE(Grid *aptgrid, PtrParamNode paramNode ) :
    BasePDE(paramNode),
    ptgrid_(aptgrid),
    numPDENodes_(0),
    numElems_(0),
    subType_(),
    numCouplingBcs_(0),
    nonLin_(false),
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
    isaxi_(domain->IsAxisymmetric()),
    isComplex_(false),    
    needSolPrev_(false),
    sol_(NULL),  
    solVec_(NULL),
    rhsVec_(NULL),
    solPrev_(NULL),
    solVecPrev_(NULL),
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

  const Vector<Double>& StdPDE::getDeriv(DERIVType derivType) const {
  
    if ( TS_alg_ == NULL ) {
      EXCEPTION( pdename_ << ":getDeriv: No derivative defined for this PDE" );
    }
    return TS_alg_->GetDeriv(derivType);
  }
  
  const Vector<Double>& StdPDE::getOld(TIMEStepType timeStepType) const {

    if ( TS_alg_ == NULL ) {
      EXCEPTION( pdename_ << ":getOld: No time stepping defined for this PDE");
    }
    return TS_alg_->GetOld(timeStepType);
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

  // calculates L2-norm of RHS regarding dirichlet entries due to penalty 
  // formulation by setting them 0
  Double StdPDE::RhsL2Norm(Vector<Double>& actRHS)
  {

    Integer eqnNr;
  
    // Eliminate inhom. dirichlet node from RHS (due to penalty formulation)
    for ( UInt i = 0; i < idBcs_.GetSize(); i++ ) {

      // Get grip of current idBC
      InhomDirichletBc const & actBc = *idBcs_[i];

      // Get entity iterator
      EntityIterator it = actBc.entities->GetIterator();
      
      for ( it.Begin(); !it.IsEnd(); it++ ) {
        eqnNr = eqnMap_->GetEqn( *actBc.result, it, actBc.dof );
        if ( eqnNr != 0 ) {
          actRHS[(eqnNr-1)] = 0.0;
        }
      }
    } 
    

    return actRHS.NormL2();
  }


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


  // real valued method (for TRANSIENT and STATIC)
  void StdPDE::GetSolVecOfElement( Vector<Double>& elemSol,
                                   const EntityIterator& it,
                                   shared_ptr<ResultInfo> res ) {


    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *res, it );


    elemSol.Resize( eqns.GetSize() );
    elemSol.Init(0);
    NodeStoreSol<Double> * solhelp = 
      dynamic_cast<NodeStoreSol<Double>*>(sol_);
    Vector<Double> & sol = solhelp->GetAlgSysVector();
    
    for( UInt i = 0; i < eqns.GetSize(); i++ ) {
      if ( eqns[i] != 0 ) {
        elemSol[i] = sol[abs(eqns[i])-1];
      } else {
        elemSol[i] = 0.0;
      }
     }
  }


  // complex valued method (for HARMONIC)
  void StdPDE::GetSolVecOfElement( Vector<Complex>& elemSol,
                                   const EntityIterator& it,
                                   shared_ptr<ResultInfo> res ) {


    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *res, it );


    elemSol.Resize( eqns.GetSize() );
    elemSol.Init( Complex(0.0, 0.0) );
    NodeStoreSol<Complex> * solhelp = 
      dynamic_cast<NodeStoreSol<Complex>*>(sol_);
    Vector<Complex> & sol = solhelp->GetAlgSysVector();
    
    for( UInt i = 0; i < eqns.GetSize(); i++ ) {
      if ( eqns[i] != 0 ) {
        elemSol[i] = sol[abs(eqns[i])-1];
      } else {
        elemSol[i] = Complex(0.0, 0.0);
      }
     }
  }


  // complex valued method (for HARMONIC)
  void StdPDE::GetRHSVecOfElement( Vector<Double>& elemRHS,
                                   const EntityIterator& it,
                                   shared_ptr<ResultInfo> res ) {


    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *res, it );


    elemRHS.Resize( eqns.GetSize() );
    elemRHS.Init( 0.0 );
    Vector<Double> rhs;
    algsys_->GetRHSVal(rhs);


    for( UInt i = 0; i < eqns.GetSize(); i++ ) {
      if ( eqns[i] != 0 ) {
        elemRHS[i] = rhs[abs(eqns[i])-1];
      } else {
        elemRHS[i] = 0.0;
      }
     }
  }


  // complex valued method (for HARMONIC)
  void StdPDE::GetRHSVecOfElement( Vector<Complex>& elemRHS,
                                   const EntityIterator& it,
                                   shared_ptr<ResultInfo> res ) {


    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *res, it );


    elemRHS.Resize( eqns.GetSize() );
    elemRHS.Init( Complex(0.0, 0.0) );
    Vector<Complex> rhs;
    algsys_->GetRHSVal(rhs);


    for( UInt i = 0; i < eqns.GetSize(); i++ ) {
      if ( eqns[i] != 0 ) {
        elemRHS[i] = rhs[abs(eqns[i])-1];
      } else {
        elemRHS[i] = Complex(0.0, 0.0);
      }
     }
  }

  
  // real valued method (for TRANSIENT)
  void StdPDE::GetDerivSolVecOfElement(Vector<Double>& sol,
                                       const EntityIterator& it,
                                       shared_ptr<ResultInfo> res) {

    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *res, it );
    sol.Resize( eqns.GetSize() );
    sol.Init( 0.0 ); 
    
    if (  analysistype_ == TRANSIENT) {
      const Vector<Double> & sol_der1 = getDeriv(FIRST_DERIV);
        
      for( UInt i = 0; i < eqns.GetSize(); i++ ) {
        if ( eqns[i] != 0 ) {
          sol[i] = sol_der1[abs(eqns[i])-1];
        } else {
          sol[i] = 0.0;
        }
      }
    }
  }



  // complex valued method (for HARMONIC)
  void StdPDE::GetDerivSolVecOfElement(Vector<Complex>& sol,
                                       const EntityIterator& it,
                                       shared_ptr<ResultInfo> res) {


    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *res, it );
    
    sol.Resize( eqns.GetSize() );
    sol.Init( 0.0 );
    
    // we obtain from assemble: frequency =  2*PI*actFreq
    Double omega = solveStep_->GetActFreq() * 2 * PI;
    Complex jomega = Complex(0.0,omega);
    
    if ( analysistype_ == HARMONIC ) {
      NodeStoreSol<Complex> * solhelp = 
          dynamic_cast<NodeStoreSol<Complex>*>(sol_);
      const Vector<Complex> & solAtNode = solhelp->GetAlgSysVector();

      for( UInt i = 0; i < eqns.GetSize(); i++ ) {
        if ( eqns[i] != 0 ) {
          sol[i] = jomega * solAtNode[abs(eqns[i])-1];
        } else {
          sol[i] = Complex(0.0, 0.0);
        }
      }
    }
  }



  // real valued method (for TRANSIENT)
  void StdPDE::GetDeriv2SolVecOfElement( Vector<Double>& sol,
                                         const EntityIterator& it,
                                         shared_ptr<ResultInfo> res) {


    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *res, it );
    
    sol.Resize( eqns.GetSize() );
    sol.Init( 0.0 );
    
    if ( analysistype_ == TRANSIENT ) {
      const Vector<Double> & sol_der2 = getDeriv(SECOND_DERIV);
      for( UInt i = 0; i < eqns.GetSize(); i++ ) {
        if ( eqns[i] != 0 ) {
          sol[i] = sol_der2[abs(eqns[i])-1];
        } else {
          sol[i] = 0.0;
        }
      }
    }
  }


  // real valued method (for HARMONIC)
  void StdPDE::GetDeriv2SolVecOfElement( Vector<Complex>& sol,
                                         const EntityIterator& it,
                                         shared_ptr<ResultInfo> res) {


    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *res, it );
    
    sol.Resize( eqns.GetSize() );
    sol.Init( 0.0 );
    
    // we obtain from assemble: frequency 
    Double omega = solveStep_->GetActFreq() * 2 * PI;
    
    if ( analysistype_ == HARMONIC ) {
      NodeStoreSol<Complex> * solhelp = 
          dynamic_cast<NodeStoreSol<Complex>*>(sol_);
      Vector<Complex> & solAtNode = solhelp->GetAlgSysVector();

      for( UInt i = 0; i < eqns.GetSize(); i++ ) {
        if ( eqns[i] != 0 ) {
          sol[i] = - omega * omega *solAtNode[abs(eqns[i])-1];
        } else {
          sol[i] = Complex(0.0, 0.0);
        }
      }
    }
  }
  
  void StdPDE::GetAnyDerivSolVecOfElement(Vector<Double>& sol, const EntityIterator& it, shared_ptr<ResultInfo> res, DERIVType derivative){
    switch(derivative){
    case NO_DERIVTYPE: GetSolVecOfElement(sol, it, res); break;
    case FIRST_DERIV: GetDerivSolVecOfElement(sol, it, res); break;
    case SECOND_DERIV: GetDeriv2SolVecOfElement(sol, it, res); break;
    }
  }

  void StdPDE::GetAnyDerivSolVecOfElement(Vector<Complex>& sol, const EntityIterator& it, shared_ptr<ResultInfo> res, DERIVType derivative){
    switch(derivative){
    case NO_DERIVTYPE: GetSolVecOfElement(sol, it, res); break;
    case FIRST_DERIV: GetDerivSolVecOfElement(sol, it, res); break;
    case SECOND_DERIV: GetDeriv2SolVecOfElement(sol, it, res); break;
    }
  }

  // real valued method (for TRANSIENT ): returns previous solution
  //                                      for element
  void StdPDE::GetPrevSolVecOfElement( Vector<Double>& elemSol,
                                       const EntityIterator& it,
                                       shared_ptr<ResultInfo> res ) {

    if ( solPrev_ == NULL ) 
      EXCEPTION("Previous Solution not defined");

    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *res, it );


    elemSol.Resize( eqns.GetSize() );
    elemSol.Init(0);
    NodeStoreSol<Double> * solhelp = 
      dynamic_cast<NodeStoreSol<Double>*>(solPrev_);
    Vector<Double> & sol = solhelp->GetAlgSysVector();
    
    for( UInt i = 0; i < eqns.GetSize(); i++ ) {
      if ( eqns[i] != 0 ) {
        elemSol[i] = sol[abs(eqns[i])-1];
      } else {
        elemSol[i] = 0.0;
      }
     }
  }

  //stores an algsys_ vector into an StdVector
  void StdPDE::StoreAlgsysToVec(Vector<Double>& vec, Double * pt) {


    //const UInt numElems = numPDENodes_ * dofspernode_;
    UInt numElems = eqnMap_->GetNumEqns();

    vec.Resize(numElems);
    for (UInt i=0; i<numElems; i++) {
      vec[i] = pt[i];
    }
  }

  SingleVector *  StdPDE::GetPrevSolutionVector() {
    return solVecPrev_;
  }

  shared_ptr<ResultInfo> StdPDE::GetResultInfo( SolutionType solType ) {
    
    shared_ptr<ResultInfo> res;
    for( UInt i = 0; i < results_.GetSize(); i++ ) {
      if( results_[i]->resultType == solType) {
        res = results_[i];
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
     std::map<RegionIdType, GridDisplData>::iterator gridDisplIter \
       = gridDisplData_.begin();
     for (; gridDisplIter != gridDisplData_.end(); ++gridDisplIter)
     {
       RegionIdType regId = gridDisplIter->first;
       GridDisplData gridDispData = gridDisplIter->second;
       ResultHandler* resultHandler = domain->GetResultHandler();
       shared_ptr<BaseResult> gridDisplacement = resultHandler->GetResult( gridDispData.fileName4GridDisplacements_,
           1,
           step,
           gridDispData.solType,        
           ptgrid_->GetRegion().ToString( regId ));

       Result<Double> *result =
         dynamic_cast<Result<Double>*>(&(*gridDisplacement));
       if (result == NULL)
       {
         EXCEPTION("Cannot read result 'Grid-Displacements' from input id '"
             <<  gridDispData.fileName4GridDisplacements_ << "'");
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

} // end of namespace
