#include "StdPDE.hh"

#include "Utils/vector.hh"
#include "Driver/stdSolveStep.hh"

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
      nonLin_(FALSE),
      incStopCrit_(1e-2), 
      residualStopCrit_(1e-3),
      materialData_(NULL),
      TS_alg_(NULL),
      effectiveMass_(FALSE),
      firstTimeStepStatic_(TRUE),
      isAlwaysStatic_(FALSE),
      isaxi_(FALSE),
      isComplex_(FALSE),
      sol_(NULL),     
      isIncrFormulation_(FALSE),
      ComputeRHSforHarm_(FALSE),
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
    eqnData_    = NULL;
    
    // =====================================================================
    // set analysis parameters
    // =====================================================================
    couplingBCsCounter_ = 0;
    numDirichletBCs_ = 0;
    isIterCoupled_ = FALSE;
    updateCouplingBCs_ = FALSE;
    dim_ = ptgrid_->GetDim();
    geoUpdate_ = FALSE;
    iterCoupledCounter_ = 0;
    effectiveMass_ = FALSE;

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

  void StdPDE::setBCs_id_phase_(UInt i, Double & phase) {

    ENTER_FCN("setBSs_PDE::setBCs_id_phase");
    
    if (bcs_id_phase_.GetSize() <= i)
      Error("no such index in Vector bcs_id_phase_",__FILE__,__LINE__);
    bcs_id_phase_[i]=phase;
  }

  // ======================================================
  // GRID SECTION (Meshing, ...) 
  // ======================================================

  void StdPDE::GetElemCoords( const StdVector<UInt> connect, 
                              Matrix<Double> &coordMat ) {

    ENTER_FCN( "StdPDE::GetElemCoords" );

    UInt pdeNode;
    ptgrid_->GetElemNodesCoord(coordMat, connect);
  
    if (deltCoords_.GetSizeRow() != 0 && geoUpdate_ == TRUE) {
      for (UInt i=0; i<coordMat.GetSizeRow(); i++) {
        for (UInt j=0; j<coordMat.GetSizeCol(); j++) {
          pdeNode = eqnData_->Mesh2PDENode(connect[j]);
          coordMat(i,j) += deltCoords_(i, pdeNode - 1);
        }
      }
    }
  }

  // calculates L2-norm of RHS regarding dirichlet entries due to penalty 
  // formulation by setting them 0
  Double StdPDE::RhsL2Norm(Vector<Double>& actRHS)
  {
    ENTER_FCN( "StdSolveStep::RhsL2Norm" );

    UInt dof, eqnDof;
    Integer eqnNr;
  
    StdVector<UInt> nodes;
  
    UInt dofsPerEQN = eqnData_->GetNumDofsPerEQN();

    // Eliminate homogeneous dirichlet node from RHS (due to penalty formulation).
    // In the case of a block system, there might be still some homogeneous dirichlet
    // entries left
    for ( UInt i = 0; i < bcs_hd_.GetSize(); i++ ) {
      dof = 1;
      if ( dofspernode_ > 1 ) {
        dof = domain->GetCoordSystem()->GetVecComponent( homDirichDof_[i] );
      }

      ptgrid_->GetNodesByName( nodes, bcs_hd_[i] );

      for ( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
        eqnData_->Node2EQN(nodes[iNode],dof,eqnNr,eqnDof);
        if ( eqnNr != 0 ) {
          actRHS[(eqnNr-1)*dofsPerEQN + eqnDof-1] = 0.0;
        }
      }
    }

    // Eliminate inhom. dirichlet node from RHS (due to penalty formulation)
    for ( UInt i = 0; i < bcs_id_.GetSize(); i++ ) {
      dof = 1;
      if ( dofspernode_ > 1 ) {
        dof = domain->GetCoordSystem()->GetVecComponent( inhomDirichDof_[i] );
      }

      ptgrid_->GetNodesByName( nodes, bcs_id_[i] );

      for ( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
        eqnData_->Node2EQN( nodes[iNode], dof, eqnNr, eqnDof );
        if ( eqnNr != 0 ) {
          actRHS[(eqnNr-1)*dofsPerEQN + eqnDof-1] = 0.0;
        }
      }
    } 

    return actRHS.NormL2();
  }

  void StdPDE::WriteRestart(const UInt nstep) {

    ENTER_FCN( "StdPDE::WriteRestart" );

    if (pdename_=="mechanic" || pdename_=="acoustic")
      {
        std::string simName = commandLine->GetSimName();
        std::string restartFileName_ = simName+"_"+pdename_+".restart";
        std::ofstream writeTo(restartFileName_.c_str(), std::ios_base::out|std::ios_base::trunc);
        GetMemento(memento_);
        writeTo << nstep << std::endl;
        writeTo << memento_;
      }

    else 
      {
        Error( "Restart functionality only for mechanic and acoustic tested", __FILE__, __LINE__ );
      }

  }

  void StdPDE::ReadRestart(UInt &startStep) {

    ENTER_FCN( "StdPDE::ReadRestart" );

    if (pdename_=="mechanic" || pdename_=="acoustic")
      {
        std::string simName = commandLine->GetSimName();
        std::string transFromTo = "standard";
        Double actFrequency=0.0;
        char buffer[64];

        std::string restartFileName_ = simName+"_"+pdename_+".restart";
        std::ifstream readRestart(restartFileName_.c_str(), std::ios_base::in );
        if (!readRestart) {
          Error( "Can not find the restart file from which i should read ?!?!", 
                 __FILE__, __LINE__ );
        }
        readRestart.seekg(0, std::ios_base::beg);
        readRestart.getline( buffer, 64, '\n' );
        std::sscanf(buffer,"%u", &startStep );
    
        UInt rhsSize = eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN();

        memento_.SetSize(rhsSize);
        memento_.analysisType_ = analysistype_;
    
        // Readind the restart file via an overwritten input operator >>
        readRestart >> memento_;
        readRestart.close();
        SetMemento( memento_, transFromTo, actFrequency );
      }

    else 
      {
        Error( "Restart functionality only for mechanic and acoustic tested", __FILE__, __LINE__ );
      }

  }

  void StdPDE::GetMemento(PDEMemento & memento) {

    ENTER_FCN( "StdPDE::GetMemento" );

    // first get memento of coupling object
    if (isIterCoupled_) {
      ptCoupling_->GetMemento(memento.couplingMemento_);
    }

    // then write own data to PDEMemento
    memento.analysisType_ = analysistype_;

    if (memento.sol_)
      delete memento.sol_;

    if ( analysistype_ == STATIC || analysistype_ == TRANSIENT ) {
      // --- Real values --
      // Set solution
      memento.sol_ = new Vector<Double>;
      dynamic_cast<Vector<Double>&>(*(memento.sol_)) =
        dynamic_cast<NodeStoreSol<Double>&>(*(sol_)).GetAlgSysVector();

      if (analysistype_ == TRANSIENT) {
        // Set first derivative
        memento.solDeriv1_ = getS1();   
          
        // Set second derivative
        memento.solDeriv2_ = getS2();
      }
    }
    else {

      // --- Complex values --      
      // Set solution
      memento.sol_ = new Vector<Complex>;
      dynamic_cast<Vector<Complex>&>(*(memento.sol_)) =
        dynamic_cast<NodeStoreSol<Complex>&>(*(sol_)).GetAlgSysVector();
    } 

    // now memento is initialized
    memento.isSet_ = TRUE;
  }


  void StdPDE::SetMemento( PDEMemento &memento, std::string transFromTo,
                           Double frequency) {

    ENTER_FCN( "StdPDE::SetMemento" );
  
    // if there is no information in the memento just leave
    if ( memento.isSet_ == FALSE ) {
      return;
    }
  
    if ( analysistype_ == STATIC || analysistype_ == TRANSIENT ) {
      
      // convert solution to transient StoreSolution type
      Vector<Double> & solVec = 
        (dynamic_cast<NodeStoreSol<Double> &>(*sol_)).GetAlgSysVector();
      

      if ( transFromTo == "complexToReal" ) {

        // --- transform complex values to real one --
        Vector<Complex>& mementoVec = 
          dynamic_cast<Vector<Complex>&>(*(memento.sol_));

        for ( UInt i=0; i < solVec.GetSize(); i++ ) {
          solVec[i] = mementoVec[i].real();
        }

        if (analysistype_ == TRANSIENT) {
          // Set first and second derivative
          memento.solDeriv1_.Resize(mementoVec.GetSize());
          memento.solDeriv2_.Resize(mementoVec.GetSize());

          Complex val;
          for ( UInt i=0; i<mementoVec.GetSize(); i++ ) {
            
            val   = mementoVec[i];
            memento.solDeriv1_[i] = - 2*PI*frequency * val.imag();
            memento.solDeriv2_[i] = 
              - 4 * PI * PI * frequency * frequency * val.real();
          }

          TS_alg_->SetDeriv1(memento.solDeriv1_);
          TS_alg_->SetDeriv2(memento.solDeriv2_);

        }
      }
      else {

        Vector<Double>& mementoVec = 
          dynamic_cast<Vector<Double>&>(*(memento.sol_));

        for (UInt i = 0; i< solVec.GetSize(); i++ ) {
          solVec[i] = mementoVec[i];
        }

        // if previous step was transient and the current step is also
        // then give the time derivative to the timestepping algorithm
        if (analysistype_ == TRANSIENT
            && memento.analysisType_ == TRANSIENT) {
          TS_alg_->SetDeriv1(memento.solDeriv1_);
          TS_alg_->SetDeriv2(memento.solDeriv2_);
        }
      }
    }

    else {

      // --- Complex values --      
      // Set solution

       // convert solution to transient StoreSolution type
      Vector<Complex> & solVec = 
        (dynamic_cast<NodeStoreSol<Complex> &>(*sol_)).GetAlgSysVector();
      
      Vector<Complex>& mementoVec = 
          dynamic_cast<Vector<Complex>&>(*(memento.sol_));

      for (UInt i = 0; i< solVec.GetSize(); i++ ) {
        solVec[i] = mementoVec[i];
      }
    }

    // If PDE is coupled, set the according coupling memento
    if (isIterCoupled_) {
      ptCoupling_->SetMemento(memento.couplingMemento_);
    }
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
      // crate eigenvalue solver
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
  
  UInt StdPDE::GetTimeStepCounter() {
    ENTER_FCN( "StdPDE::GetTimeStepCounter" );
    return solveStep_->GetActStep();
  }



  //=========================================================================
  // Implemented TransformSol4Slice
  //=========================================================================
  /*
    Has the task to shift the solution vector and the vectors for the first- and
    second derivatives. Furtheron, the entries of these vectors which are not
    calculated for the new slice are initialized with zeros:
  */
  //void StdPDE::TransformSol4Slice(UInt &  nodeShift, UInt & shiftFactor, 
  //				  const UInt flag) {
  void StdPDE::TransformSol4Slice(UInt & shiftFactor, UInt & nodeShift,
		UInt & elemgrid, Double &  meshsize, const UInt flag){

    ENTER_FCN( "StdPDE::TransformSol4Slice" );

    ptgrid_->TransformGridStruct(shiftFactor, nodeShift,
					elemgrid, meshsize, flag);
   
    if (flag) 
      return;

    //perform the transformation of solution
    NodeStoreSol<Double>  & solHelp = dynamic_cast<NodeStoreSol<Double>&>(*sol_);
  
    Vector<Double> & sol = solHelp.GetAlgSysVector();
    Vector<Double> solD1 = getS1();
    Vector<Double> solD2 = getS2();
    
    //get equation numbers
    UInt dof = 1;
    Integer eqnNrFrom, eqnNrTo; 
    UInt eqnDof, nodeFrom;

//     for (UInt node=1; node <= 10*shiftFactor; node++) {
//       eqnData_->Node2EQN(node,dof,eqnNrTo, eqnDof);
//       sol[eqnNrTo-1]   = 0;
//       solD1[eqnNrTo-1] = 0;
//       solD2[eqnNrTo-1] = 0;
//     }

    for (UInt node=1; node <=numPDENodes_-nodeShift; node++) {
      eqnData_->Node2EQN(node,dof,eqnNrTo, eqnDof);
      nodeFrom = node + nodeShift;
      eqnData_->Node2EQN(nodeFrom,dof,eqnNrFrom, eqnDof);
      
      sol[eqnNrTo-1]   = sol[eqnNrFrom-1];
      solD1[eqnNrTo-1] = solD1[eqnNrFrom-1];
      solD2[eqnNrTo-1] = solD2[eqnNrFrom-1];
    }

    // set remaining nodes to zero
    for (UInt node=numPDENodes_-nodeShift+1; node <=numPDENodes_; node++) {
      eqnData_->Node2EQN(node,dof,eqnNrTo, eqnDof);
      sol[eqnNrTo-1]   = 0;
      solD1[eqnNrTo-1] = 0;
      solD2[eqnNrTo-1] = 0;
    }
    
    TS_alg_->SetDeriv1(solD1);
    TS_alg_->SetDeriv2(solD2);
    
    
    //    WriteResultsInFile(0,0);
    
  }


  /*
   * save important nodes for a post analysis
   * ATENTION: The first time SaveNodes is called in
   * transient4slicedriver timeStep has the value meshsize
   */
  void StdPDE::SaveNodes(const UInt shiftFactor, const Double timeStep,
			 const UInt numShift, const Integer nodeShift, 
			 const UInt maxnumelemz_) {
  
    ENTER_FCN( "StdPDE::TransformSol4Slice" );
    
    NodeStoreSol<Double>  & solHelp = dynamic_cast<NodeStoreSol<Double>&>(*sol_);
    Vector<Double> & sol = solHelp.GetAlgSysVector();

    Integer dof = 1;
    Integer filesave, eqn;
    UInt local, global, eqnDof;

    if(nodeShift == -1){
      //First time the function SaveNodes is called
      std::string S="mkdir -p saveNodes";
      system(S.c_str());

      //Write the actual meshsize in mesh.dat
      std::fstream x("saveNodes/mesh.dat", std::ios::out|std::ios::app);
      x << timeStep << '\t' << shiftFactor << '\n';
      return;
    }
  
    // Calculate the global Node number
    for(UInt i = 0; i<=maxnumelemz_; i++){
      //access the global node number
      local  = i*shiftFactor+1;
      global = local + nodeShift*numShift;
      filesave = i + (nodeShift*numShift)/shiftFactor;
  
      eqnData_->Node2EQN(local,dof,eqn, eqnDof);
      //Store the solution to the right position of the solution Vector.
      //testing if shifts have been made
  
      // First of all built the whole filename
      std::string namePrefix="saveNodes/focuse";
      std::string namePostfix = ".save";
      std::string totalName;

      // Create complete filename
      std::string quantString;

      //Enum2String((double) global, quantString);
      //should be a "good" conversion from Integer to string

      std::ostringstream o;

      o<<filesave;
      totalName = namePrefix + o.str();
      totalName += "node";
      totalName += namePostfix;
  
      //Öffnen der Datei
      std::fstream x(totalName.c_str(), std::ios::out|std::ios::app);

      x << timeStep << '\t'  <<  sol[eqn-1] << '\n'; 
    }
  }

  // ======================================================
  // ALGSYS SECTION (SOLVER, ...) 
  // ======================================================

  Double StdPDE::GetFracDampMatrixCoeff(UInt actSD) {
    
    ENTER_FCN( "StdPDE::GetFracDampMatrixCoeff" );

    Double coeff;

    // pre factor of fractional derivative (same for all algorithms)
    Double y  = materialData_[actSD].GetDampingBeta();
    Double dt = TS_alg_->GetTimeStep();
    coeff = std::exp(-(y-1.0) * std::log(dt));

    // needed for formulation with only MASS and STIFFNESS matrix
    // pre factor of Newmark time stepping scheme
    Double beta = TS_alg_->GetNewmarkBeta();
    coeff *= 1.0 / (beta*dt*dt);

//     std::cerr<< std::endl << "Parameters in GetFracDampMatrixCoeff are:"
//              << std::endl << "actSD: " << actSD
//              <<	std::endl << "y    = " << y
//              <<	std::endl << "dt   = " << dt
//              <<	std::endl << "beta = " << beta
//              << std::endl << "coeff= " << coeff << std::endl << std::endl;

    return coeff;
  }
  

  // real valued method (for TRANSIENT and STATIC)
  void StdPDE::GetSolVecOfElement( Vector<Double>& elemSol,
                                   StdVector<UInt>& connecth ) {

    ENTER_FCN( "StdPDE::GetSolVecOfElement" );

    // displacement of element nodes
    elemSol.Resize(dofspernode_ * connecth.GetSize());
    elemSol.Init(0);
    Integer eqnNr; 
    UInt eqnDof;
    UInt dofsPerEQN = eqnData_->GetNumDofsPerEQN();

    NodeStoreSol<Double> * solhelp = 
      dynamic_cast<NodeStoreSol<Double>*>(sol_);
    Vector<Double> sol = solhelp->GetAlgSysVector();
  
    for(UInt actNode=0; actNode<connecth.GetSize(); actNode++) {
      for(UInt actDof=0; actDof < dofspernode_; actDof++) {
        eqnData_->Node2EQN(connecth[actNode],actDof+1,eqnNr,eqnDof);
        if (eqnNr!= 0) {
          elemSol[actDof + actNode*dofspernode_] =
            sol[eqnDof-1 + dofsPerEQN*(abs(eqnNr-1))];
        }
        else {
          elemSol[actDof + actNode*dofspernode_] = 0.0;
        }
      }
    }
  }


  // complex valued method (for HARMONIC)
  void StdPDE::GetSolVecOfElement( Vector<Complex>& elemSol,
                                   StdVector<UInt>& connecth ) {

    ENTER_FCN( "StdPDE::GetSolVecOfElement" );

    // displacement of element nodes
    elemSol.Resize(dofspernode_ * connecth.GetSize());
    elemSol.Init(0);
    Integer eqnNr; 
    UInt eqnDof;
    UInt dofsPerEQN = eqnData_->GetNumDofsPerEQN();

    NodeStoreSol<Complex> * solhelp = 
      dynamic_cast<NodeStoreSol<Complex>*>(sol_);
    Vector<Complex> sol = solhelp->GetAlgSysVector();
  
    for(UInt actNode=0; actNode<connecth.GetSize(); actNode++) {
      for(UInt actDof=0; actDof < dofspernode_; actDof++) {
        eqnData_->Node2EQN(connecth[actNode],actDof+1,eqnNr,eqnDof);
        if (eqnNr!= 0) {
          elemSol[actDof + actNode*dofspernode_] =
            sol[eqnDof-1 + dofsPerEQN*(abs(eqnNr-1))];
        }
        else {
          elemSol[actDof + actNode*dofspernode_] = 0.0;
        }
      }
    }
  }

  
  // real valued method (for TRANSIENT)
  void StdPDE::GetDerivSolVecOfElement(Vector<Double>& sol,
                                       StdVector<UInt>& connecth) {

    ENTER_FCN( "StdPDE::GetDerivSolVecOfElement" );

    // displacement of element nodes
    sol.Resize(dofspernode_ * connecth.GetSize());
    sol.Init(0);
    Integer eqnNr;
    UInt  eqnDof;
    UInt dofsPerEQN = eqnData_->GetNumDofsPerEQN();
  
    if ( analysistype_ == TRANSIENT ) {
      const Vector<Double> & sol_der1 = getS1();
    
      for( UInt actNode = 0; actNode < connecth.GetSize(); actNode++ ) {
        for( UInt actDof = 0; actDof < dofspernode_; actDof++ ) {
          eqnData_->Node2EQN(connecth[actNode],actDof+1,eqnNr,eqnDof);
          if ( eqnNr != 0 ) {
            sol[actDof + actNode*dofspernode_] =
              sol_der1[eqnDof-1 + dofsPerEQN*(abs(eqnNr-1))];
          }
          else {
            sol[actDof + actNode*dofspernode_] = 0.0;
          }
        }
      }
    }
  }

  // complex valued method (for HARMONIC)
  void StdPDE::GetDerivSolVecOfElement(Vector<Complex>& sol,
                                       StdVector<UInt>& connecth) {

    ENTER_FCN( "StdPDE::GetDerivSolVecOfElement" );

    // displacement of element nodes
    sol.Resize(dofspernode_ * connecth.GetSize());
    sol.Init(0);
    Integer eqnNr;
    UInt  eqnDof;
    UInt dofsPerEQN = eqnData_->GetNumDofsPerEQN();

    // we obtain from assemble: frequency =  2*PI*actFreq
    Double omega = assemble_->GetActFrequency();
    Complex jomega = Complex(0.0,omega);

    if ( analysistype_ == HARMONIC ) {
      NodeStoreSol<Complex> * solhelp = 
	dynamic_cast<NodeStoreSol<Complex>*>(sol_);
      Vector<Complex> solAtNode = solhelp->GetAlgSysVector();

      for( UInt actNode = 0; actNode < connecth.GetSize(); actNode++ ) {
        for( UInt actDof = 0; actDof < dofspernode_; actDof++ ) {
          eqnData_->Node2EQN(connecth[actNode],actDof+1,eqnNr,eqnDof);
          if ( eqnNr != 0 ) {
            sol[actDof + actNode*dofspernode_] =
              solAtNode[eqnDof-1 + dofsPerEQN*(abs(eqnNr-1))] * jomega;
          }
          else {
            sol[actDof + actNode*dofspernode_] = 0.0;
          }
        }
      }
    }
  }



  // real valued method (for TRANSIENT)
  void StdPDE::GetDeriv2SolVecOfElement( Vector<Double>& sol,
                                         StdVector<UInt>& connecth ) {

    ENTER_FCN( "StdPDE::GetDeriv2SolVecOfElement" );

    // displacement of element nodes
    sol.Resize(dofspernode_ * connecth.GetSize());
    sol.Init(0);
    Integer eqnNr; 
    UInt  eqnDof;
    UInt dofsPerEQN = eqnData_->GetNumDofsPerEQN();

    if (analysistype_ == TRANSIENT) {
      const Vector<Double> & sol_der2 = getS2();
    
      for( UInt actNode = 0; actNode < connecth.GetSize(); actNode++ ) {
        for( UInt actDof = 0; actDof < dofspernode_; actDof++ ) {
          eqnData_->Node2EQN(connecth[actNode],actDof+1,eqnNr,eqnDof);
          if (eqnNr!= 0) {
            sol[actDof + actNode*dofspernode_] =
              sol_der2[eqnDof-1 + dofsPerEQN*(abs(eqnNr-1))];
          }
          else {
            sol[actDof + actNode*dofspernode_] = 0.0;
          }
        }
      }
    }
  }


  // real valued method (for HARMONIC)
  void StdPDE::GetDeriv2SolVecOfElement( Vector<Complex>& sol,
                                         StdVector<UInt>& connecth ) {

    ENTER_FCN( "StdPDE::GetDeriv2SolVecOfElement" );

    // displacement of element nodes
    sol.Resize(dofspernode_ * connecth.GetSize());
    sol.Init(0);
    Integer eqnNr; 
    UInt  eqnDof;
    UInt dofsPerEQN = eqnData_->GetNumDofsPerEQN();

    // we obtain from assemble: frequency =  2*PI*actFreq
    Double omega = assemble_->GetActFrequency();

    if ( analysistype_ == HARMONIC ) {
      NodeStoreSol<Complex> * solhelp = 
	dynamic_cast<NodeStoreSol<Complex>*>(sol_);
      Vector<Complex> solAtNode = solhelp->GetAlgSysVector();

      for( UInt actNode = 0; actNode < connecth.GetSize(); actNode++ ) {
        for( UInt actDof = 0; actDof < dofspernode_; actDof++ ) {
          eqnData_->Node2EQN(connecth[actNode],actDof+1,eqnNr,eqnDof);
          if (eqnNr!= 0) {
            sol[actDof + actNode*dofspernode_] = - omega * omega *
	      solAtNode[eqnDof-1 + dofsPerEQN*(abs(eqnNr-1))];
          }
          else {
            sol[actDof + actNode*dofspernode_] = 0.0;
          }
        }
      }
    }
  }


  void StdPDE::GetDerivSolOfElement( Matrix<Double>& sol,
                                     StdVector<UInt>& connect_PDE ) {

    ENTER_FCN( "StdPDE::GetDerivSolOfElement" );

    // displacement of element nodes
    sol.Resize(dofspernode_, connect_PDE.GetSize());

    const Vector<Double>& sol_der1 = getS1();  

    for( UInt actNode = 0; actNode < connect_PDE.GetSize(); actNode++ ) {
      for( UInt actDof = 0; actDof < dofspernode_; actDof++) {
        sol[actDof][actNode] =
          sol_der1[actDof + dofspernode_*(connect_PDE[actNode]-1)];
      }
    }
  }

  void StdPDE::sortStresses(Vector<Double>& unsorted, Vector<Double>& sorted){

    ENTER_FCN( "StdPDE::SortStresses" );

    //soring according to capa (unv) notation
    //our notation is Voigt: Txx Tyy Tzz Tyz Txz Txy

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
    UInt numElems = eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN();

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
    if ( analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC ) {
      subDomVolComplex.Resize(surfRegions.GetSize());
    }
    else {
      subDomVolReal.Resize(surfRegions.GetSize());
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
      if ( analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC ) {
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
        GetElemCoords(connecth, ptSurfCoord);

        //get the deformed solution
        if (analysistype_ == HARMONIC || analysistype_==MULTIHARMONIC) {
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
    
    if ( analysistype_ == HARMONIC || analysistype_== MULTIHARMONIC ) {
      subDomVolReal.Resize(surfRegions.GetSize());
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
