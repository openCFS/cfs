#include "StdPDE.hh"

#include "Utils/vector.hh"
#include "Driver/stdSolveStep.hh"
#include "Driver/transientdriver.hh"

#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"

// headers for Paramhandling
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/CFSOLASParams.hh"
#include "Driver/assemble.hh"

namespace CoupledField {


  StdPDE::StdPDE(Grid *aptgrid, WriteResults *aOutFile, TimeFunc *aTimeFunc)
    : BasePDE(),
      numPDENodes_(0),
      numElems_(0),
      nonLin_(false),
      nonLinMaterial_(false),
      isHysteresis_(false),
      incStopCrit_(1e-2), 
      residualStopCrit_(1e-3),
      TS_alg_(NULL),
      effectiveMass_(false),
      diagMass_(false),
      firstTimeStepStatic_(true),
      isAlwaysStatic_(false),
      isaxi_(false),
      isComplex_(false),
      fracDamping_(false),
      sol_(NULL),     
      isIncrFormulation_(false),
      ComputeRHSforHarm_(false),
      solveStep_(NULL) {

    ENTER_FCN( "StdPDE::StdPDE");

    // =====================================================================
    // set file pointers
    // =====================================================================
    outFile_    = aOutFile;
    ptTimeFunc_ = aTimeFunc;
    ptgrid_     = aptgrid;
    assemble_   = NULL;
    algsys_     = NULL;
    
    // =====================================================================
    // set analysis parameters
    // =====================================================================
    couplingBCsCounter_ = 0;
    isIterCoupled_ = false;
    updateCouplingBCs_ = false;
    dim_ = ptgrid_->GetDim();
    iterCoupledCounter_ = 0;
    effectiveMass_ = false;

    // =====================================================================
    // various parameters
    // =====================================================================
    needsAlgsys_ = true;


  }
  
  StdPDE::~StdPDE() {
    ENTER_FCN( "StdPDE::~StdPDE" );
  }



  const Vector<Double>& StdPDE::getS1() const {
  
    ENTER_FCN( "StdPDE::getS1" );
  
    if ( TS_alg_ != NULL ) {
      return TS_alg_->GetDeriv1();
    }
    else {
      (*error) << pdename_ << ":getS1: No timestepping defined for this PDE";
      Error( __FILE__, __LINE__ );

      // Only a dummy line for compiler
      return TS_alg_->GetDeriv1();      
    }
  }
  
  
  const Vector<Double>& StdPDE::getS2() const {
    
    ENTER_FCN( "StdPDE::getS2" );
    
    if ( TS_alg_ != NULL ) {
      return TS_alg_->GetDeriv2();
    }
    else {
      (*error) << pdename_ << ":getS2: No timestepping defined for this PDE";
      Error( __FILE__, __LINE__ );

      // Only a dummy line for compiler
      return TS_alg_->GetDeriv2();
    }
  }

  // ======================================================
  // GRID SECTION (Meshing, ...) 
  // ======================================================

  // calculates L2-norm of RHS regarding dirichlet entries due to penalty 
  // formulation by setting them 0
  Double StdPDE::RhsL2Norm(Vector<Double>& actRHS)
  {
    ENTER_FCN( "StdPDE::RhsL2Norm" );

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
    
    ENTER_FCN( "StdPDE::CreateMatrices_Solver" );
    
    // ==============================
    // create and initialize matrices 
    // ==============================

    // create algebraic system and intialize matrices
    SETPROFILE("Before CreatLinSys()");
    algsys_->CreateLinSys();
    SETPROFILE("After CreatLinSys()");
    algsys_->InitMatrix();
    
    // Check for analysistype
    if ( analysistype_ != EIGENFREQUENCY ) {
      
      // create solver and preconditioner
      SETPROFILE("Before CreateSolver()");
      algsys_->CreateSolver();
      SETPROFILE("Before CreatePrecond()");
      algsys_->CreatePrecond();
      
    } else {
      // create eigenvalue solver
      algsys_->CreateEigenSolver();
    }
        
    // now reset AlgebraicSystem 
    algsys_->InitRHS();
    algsys_->InitSol();
    SETPROFILE("-- Finished CreateMatrices_Solver--");
  }


  BaseSolveStep * StdPDE::GetSolveStep() {
    
    ENTER_FCN( "StdPDE::GetSolveStep" );
    return solveStep_;
  }
  

  // ======================================================
  // ALGSYS SECTION (SOLVER, ...) 
  // ======================================================

  Double StdPDE::GetFracDampMatrixCoeff(RegionIdType regionId) {
    
    ENTER_FCN( "StdPDE::GetFracDampMatrixCoeff" );

    Double coeff;

    // pre factor of fractional derivative (same for all algorithms)
    TransientDriver * driver = NULL; 
    driver = dynamic_cast<TransientDriver*>(domain->GetSingleDriver() );
    
    if( driver == NULL) {
      (*error) << "Fractional damping only possible for transient simulation!";
      Error( __FILE__, __LINE__ );
    }
    Double dt = driver->GetTimeStep();

    Double y;
    materials_[regionId]->GetScalar(y,FRACTIONAL_EXPONENT,REAL);

    coeff = std::exp(-(y-1.0) * std::log(dt));

    // needed for formulation with only MASS and STIFFNESS matrix
    // pre factor of Newmark time stepping scheme
    Double beta = TS_alg_->GetNewmarkBeta();
    coeff *= 1.0 / (beta*dt*dt);

#ifdef DEBUG
    (*debug) << std::endl << "Parameters in GetFracDampMatrixCoeff are:"
             << std::endl << "regionId: " << regionId
             << std::endl << "y    = " << y
             << std::endl << "dt   = " << dt
             << std::endl << "beta = " << beta
             << std::endl << "coeff= " << coeff << std::endl << std::endl;
#endif

    return coeff;
  }
  

  // real valued method (for TRANSIENT and STATIC)
  void StdPDE::GetSolVecOfElement( Vector<Double>& elemSol,
                                   const EntityIterator& it ) {

    ENTER_FCN( "StdPDE::GetSolVecOfElement" );

    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *results_[0], it );


    elemSol.Resize( eqns.GetSize() );
    elemSol.Init(0);
    NodeStoreSol<Double> * solhelp = 
      dynamic_cast<NodeStoreSol<Double>*>(sol_);
    Vector<Double> sol = solhelp->GetAlgSysVector();
    
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
                                   const EntityIterator& it ) {

    ENTER_FCN( "StdPDE::GetSolVecOfElement" );
    Warning( "Implement", __FILE__, __LINE__);

    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *results_[0], it );


    elemSol.Resize( eqns.GetSize() );
    elemSol.Init( Complex(0.0, 0.0) );
    NodeStoreSol<Complex> * solhelp = 
      dynamic_cast<NodeStoreSol<Complex>*>(sol_);
    Vector<Complex> sol = solhelp->GetAlgSysVector();
    
    for( UInt i = 0; i < eqns.GetSize(); i++ ) {
      if ( eqns[i] != 0 ) {
        elemSol[i] = sol[abs(eqns[i])-1];
      } else {
        elemSol[i] = Complex(0.0, 0.0);
      }
     }
  }

  
  // real valued method (for TRANSIENT)
  void StdPDE::GetDerivSolVecOfElement(Vector<Double>& sol,
                                       const EntityIterator& it) {

    ENTER_FCN( "StdPDE::GetDerivSolVecOfElement" );
    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *results_[0], it );
    sol.Resize( eqns.GetSize() );
    sol.Init( 0.0 ); 
    
    if (  analysistype_ == TRANSIENT) {
      const Vector<Double> & sol_der1 = getS1();
        
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
                                        const EntityIterator& it) {

    ENTER_FCN( "StdPDE::GetDerivSolVecOfElement" );

    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *results_[0], it );
    
    sol.Resize( eqns.GetSize() );
    sol.Init( 0.0 );
    
    // we obtain from assemble: frequency =  2*PI*actFreq
    Double omega = solveStep_->GetActFreq() * 2 * PI;
    Complex jomega = Complex(0.0,omega);
    
    if ( analysistype_ == HARMONIC ) {
      NodeStoreSol<Complex> * solhelp = 
 	dynamic_cast<NodeStoreSol<Complex>*>(sol_);
      Vector<Complex> solAtNode = solhelp->GetAlgSysVector();

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
                                         const EntityIterator& it) {

    ENTER_FCN( "StdPDE::GetDeriv2SolVecOfElement" );

    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *results_[0], it );
    
    sol.Resize( eqns.GetSize() );
    sol.Init( 0.0 );
    
    if ( analysistype_ == TRANSIENT ) {
      const Vector<Double> & sol_der2 = getS2();
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
                                         const EntityIterator& it) {

    ENTER_FCN( "StdPDE::GetDeriv2SolVecOfElement" );

    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *results_[0], it );
    
    sol.Resize( eqns.GetSize() );
    sol.Init( 0.0 );
    
    // we obtain from assemble: frequency 
    Double omega = solveStep_->GetActFreq() * 2 * PI;
    
    if ( analysistype_ == HARMONIC ) {
      NodeStoreSol<Complex> * solhelp = 
 	dynamic_cast<NodeStoreSol<Complex>*>(sol_);
      Vector<Complex> solAtNode = solhelp->GetAlgSysVector();

      for( UInt i = 0; i < eqns.GetSize(); i++ ) {
        if ( eqns[i] != 0 ) {
          sol[i] = - omega * omega *solAtNode[abs(eqns[i])-1];
        } else {
          sol[i] = Complex(0.0, 0.0);
        }
      }
    }
  }


  void StdPDE::sortStresses(Vector<Double>& unsorted, Vector<Double>& sorted){

    ENTER_FCN( "StdPDE::SortStresses" );

    //soring according to capa (unv) notation
    //our notation is Voigt: Txx Tyy Tzz Tyz Txz Txy
    sorted.Resize(6);

    if (subType_ == "3d") {
      //Txx Txy Tyy Txz Tyz Tzz
      sorted[0] = unsorted[0];
      sorted[1] = unsorted[5];
      sorted[2] = unsorted[1];
      sorted[3] = unsorted[4];
      sorted[4] = unsorted[3];
      sorted[5] = unsorted[2];
    }
    else if (subType_ == "axi") {

      //Tphiphi 0 Trr 0 Trz Tzz
      sorted[0] = unsorted[3];
      sorted[1] = 0.0;
      sorted[2] = unsorted[0];
      sorted[3] = 0.0;
      sorted[4] = unsorted[2];
      sorted[5] = unsorted[1];
    }
    else {
      //0 0 Tyy 0 Tyz Tzz
      sorted[0] = 0.0;
      sorted[1] = 0.0;
      sorted[2] = unsorted[0];
      sorted[3] = 0.0;
      sorted[4] = unsorted[2];
      sorted[5] = unsorted[1];
    }
  }


  void StdPDE::sortStresses( Vector<Complex> &unsorted,
                             Vector<Complex> &sorted ) {
    ENTER_FCN( "StdPDE::SortStresses" );

    //soring according to capa (unv) notation
    //our notation is Voigt: Txx Tyy Tzz Tyz Txz Txy

    sorted.Resize(6);
    if (subType_ == "3d") {

      //Txx Txy Tyy Txz Tyz Tzz
      sorted[0] = unsorted[0];
      sorted[1] = unsorted[5];
      sorted[2] = unsorted[1];
      sorted[3] = unsorted[4];
      sorted[4] = unsorted[3];
      sorted[5] = unsorted[2];
    }

    else if (subType_ == "axi") {

      //Tphiphi 0 Trr 0 Trz Tzz
      sorted[0] = unsorted[3];
      sorted[1] = 0.0;
      sorted[2] = unsorted[0];
      sorted[3] = 0.0;
      sorted[4] = unsorted[2];
      sorted[5] = unsorted[1];
    }

    else {

      //0 0 Tyy 0 Tyz Tzz
      sorted[0] = 0.0;
      sorted[1] = 0.0;
      sorted[2] = unsorted[0];
      sorted[3] = 0.0;
      sorted[4] = unsorted[2];
      sorted[5] = unsorted[1];
    }
  }


  //stores an algsys_ vector into an StdVector
  void StdPDE::StoreAlgsysToVec(Vector<Double>& vec, Double * pt) {

    ENTER_FCN( "StdPDE::StoreAlgsysToVec" );

    //const UInt numElems = numPDENodes_ * dofspernode_;
    UInt numElems = eqnMap_->GetNumEqns();

    vec.Resize(numElems);
    for (UInt i=0; i<numElems; i++) {
      vec[i] = pt[i];
    }
  }

  CFSVector *  StdPDE::GetSolutionVector() {
    ENTER_FCN( "StdPDE::GetSolutionVector" );
    return solVec_;
  }

  // ======================================================
  // POSTPROCESSING  
  // ======================================================
  void StdPDE::ComputeVolDefSurf(StdVector<RegionIdType> &surfRegions, 
                                 StdVector<std::string> &strDir) {

    ENTER_FCN( "StdPDE::ComputeVolDefSurf" );

    // function uses index starting at zero
    Vector<Double>  subDomVolReal;
    Vector<Complex> subDomVolComplex;
    if ( analysistype_ == HARMONIC ) {
      subDomVolComplex.Resize(surfRegions.GetSize());
      subDomVolComplex.Init();
    }
    else {
      subDomVolReal.Resize(surfRegions.GetSize());
      subDomVolReal.Init();
    }

    UInt dof, dir;

    // Loop over all specified surfaces
    for ( UInt actSF = 0; actSF < surfRegions.GetSize(); actSF++ ) {

      // Determine co-ordinate direction for surface computation
      if ( strDir[actSF] == "ux" ) {
        dir = 1;
      }
      else if ( strDir[actSF] == "uy" ) {
        dir = 2;
      }
      else if ( strDir[actSF] == "uz" ) {
        dir = 3;
      }
      else {
        (*error) << "ComputeVolDefSurf: dof = '" << dir << "' is not "
                 << "allowed! Must be one of 'ux', 'uy' or 'uz'";
        Error( __FILE__, __LINE__ );
      }

      //we start from zero!
      dof = dir - 1;

      NodeStoreSol<Complex> * solHarmonic;
      NodeStoreSol<Double> * solTransient;
      if ( analysistype_ == HARMONIC ) {
        solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);
        subDomVolComplex[actSF] = 0;   
      }
      else {
        solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);;
        subDomVolReal[actSF] = 0;   
      }

      StdVector<SurfElem*> elemssd;
      ptgrid_->GetSurfElems(elemssd,surfRegions[actSF]);

      for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
        BaseFE * ptSurfEl = elemssd[actEl]->ptElem;
        StdVector<UInt> connecth = elemssd[actEl]->connect;
        
        Matrix<Double> ptSurfCoord;
        ptgrid_->GetElemNodesCoord( ptSurfCoord, connecth, true );

        //get the deformed solution
        if (analysistype_ == HARMONIC ) {
          Vector<Complex> disp(ptSurfEl->GetNumNodes());
          for (UInt lnode=0; lnode < ptSurfEl->GetNumNodes(); lnode++) {
            solHarmonic->Get(connecth[lnode]-1, dof, disp[lnode]);
          }
          // extract to volume element
          subDomVolComplex[actSF] += 
            ComputeVolElem(ptSurfEl,ptSurfCoord,disp); 
        }
        else {
          Vector<Double> disp(ptSurfEl->GetNumNodes());
          for (UInt lnode=0; lnode < ptSurfEl->GetNumNodes(); lnode++) {
            solTransient->Get(connecth[lnode]-1,dof,disp[lnode]);
          }
          // extract to volume element
          subDomVolReal[actSF] += ComputeVolElem(ptSurfEl,ptSurfCoord,disp); 
        }
      }
    }

    std::string resulttype = "DeformedSurfVolume";
    std::string unit = "(m^3)";
    std::string analysis;
    Double analysisVal;

    // convert region Ids to names
    StdVector<std::string> regionNames;
    ptgrid_->RegionIdToName( regionNames, surfRegions );
    
    if ( analysistype_ == HARMONIC ) {
      subDomVolReal.Resize(surfRegions.GetSize());
      subDomVolReal.Init();
      for (UInt actSF = 0; actSF < surfRegions.GetSize(); actSF++) {
        subDomVolReal[actSF] = abs( subDomVolComplex[actSF] );
      }
 
      analysis    = "Frequency:";
      analysisVal = solveStep_->GetActFreq();
      Info->WriteResult( pdename_, resulttype, regionNames, subDomVolReal, 
                         unit, analysis, analysisVal );
    }
    else {
      analysis    = "Time:";
      analysisVal = solveStep_->GetActTime();
      Info->WriteResult( pdename_, resulttype, regionNames, subDomVolReal, 
                         unit, analysis, analysisVal );
    }
  }

  
  Double StdPDE::ComputeVolElem(BaseFE * surfEl, Matrix<Double>& surfCoord, 
                                Vector<Double> disp) {

    ENTER_FCN( "StdPDE::ComputeVolElem" );

    Double elemVol;
    Double averageDis;
    UInt nrSurfNodes = surfEl->GetNumNodes();


    //compute average displacedment
    averageDis = 0;
    for (UInt i=0; i<nrSurfNodes; i++) {
      averageDis += disp[i]; 
    }
    averageDis /= (Double)nrSurfNodes;

    //compute the deformed volume    
    elemVol = averageDis * surfEl->CalcVolume(surfCoord,isaxi_);
      
    return elemVol;
  }


  Complex StdPDE::ComputeVolElem(BaseFE * surfEl, Matrix<Double>& surfCoord, 
                                 Vector<Complex> disp) {

    ENTER_FCN( "StdPDE::ComputeVolElem" );

    Complex elemVol;
    Complex averageDis;
    UInt nrSurfNodes = surfEl->GetNumNodes();


    //compute average displacedment
    averageDis = 0;
    for (UInt i=0; i<nrSurfNodes; i++) {
      averageDis += disp[i]; 
    }
    averageDis /= (Double)nrSurfNodes;

    //compute the deformed volume    
    elemVol = averageDis * surfEl->CalcVolume(surfCoord,isaxi_);
      
    return elemVol;
  }


  // ******************
  //   ReadOlasParams
  // ******************
  void StdPDE::ReadOlasParams( std::string sysName ) {

    ENTER_FCN( "StdPDE::ReaOlasParams" );

    // Log to .las file
    (*cla) <<  " --- CFS: Setting parameters for linear system '"
           << sysName << " ---" << std::endl;

    // Set parameters for OLAS
    std::string amExpert;
    params->Get( "override", amExpert, "expert" );
    CFSOLASParams::SetParams( sysName, params, olasParams_, analysistype_,
                              (amExpert == "yes") );

  }

} // end of namespace
