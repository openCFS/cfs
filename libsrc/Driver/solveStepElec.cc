#include <fstream>
#include <iostream>
#include <string>

#include "solveStepElec.hh"
#include "Forms/forms_header.hh"
#include "Utils/preisach.hh"
#include "assemble.hh"

#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"

namespace CoupledField {

  SolveStepElec::SolveStepElec(StdPDE& apde) : StdSolveStep(apde)
  {

    ENTER_FCN( "SolveStepElec::SolveStepElec" );
    doInit_ = TRUE;
  }


  SolveStepElec::~SolveStepElec() {
    ENTER_FCN( "SolveStepElec::~SolveStepElec" );
  }
 

  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  void SolveStepElec:: PreStepStatic(const UInt kstep, const Double asteptime,
                                     const Boolean reset)
  {
    ENTER_FCN( "SolveStepElec::PreStepStatic" );

    if (isIterCoupled_)     
      algsys_->InitSol();
  
    if (geoUpdate_) {
      algsys_->InitRHS();
      algsys_->InitSol();
      algsys_->InitMatrix();
      assemble_->SetReassemble();   
    }


    //for hysteresis modeling
    Boolean hystModel = FALSE;
    //  Boolean hystModel = TRUE;
    UInt numElems = PDE_.getPDE_numElems();

    if (hystModel) {
      if (kstep==1) {
        Eprevious_.Resize(numElems);
        Dprevious_.Resize(numElems);
        epsDiff_.Resize(subdoms_.GetSize(),numElems);
        Eprevious_.Init(0);
        Dprevious_.Init(0);
        epsDiff_.Init(8.854e-12);
      }

      if (doInit_) {
        //set the Preisach values; should be obtained from the xml_file
        Esat_     = 2.0e6;
        Psat_     = 0.04;
        Ec_       = 0.9e6;
        Dir_      = 2;
        isVirgin_ = TRUE;   
        hyst_ = new Preisach(numElems, Esat_, Psat_, Ec_, isVirgin_);
      
        doInit_ = FALSE;
      }
      else {
        DoUpdateHyst();
      }
    }
  }



  // time is used for a series of static calculations
  // don't get confused with REAL transient simulations!
  void SolveStepElec::SolveStepStatic(const UInt kstep, const Double asteptime,
                                      const Boolean reset) {

    ENTER_FCN( "SolveStepElec::SolveStepStatic" );
  
    lasttimecalc_ = asteptime;
    laststepcalc_ = kstep;
  
    Boolean nonLin = FALSE;
    //  Boolean nonLin = TRUE;
    if (nonLin) {
      StepStaticNonLinEpsDiff(kstep,asteptime,reset);
    }
    else {
      StepStaticLin(kstep,asteptime,reset);
    }

  }


  void SolveStepElec::StepStaticNonLinEpsDiff(const UInt kstep, const Double asteptime,
                                              const Boolean reset) {

    ENTER_FCN( "SolveStepElec::StepStaticNonLin" );

    laststepcalc_ = kstep;
    lasttimecalc_ = asteptime;

    Boolean performOneMoreStep;
    UInt iterationCounter=0;
  
    Vector<Double> newSol(eqnData_->GetNumEQNs());
    Vector<Double> oldSol(eqnData_->GetNumEQNs());
    Vector<Double> solPrev(eqnData_->GetNumEQNs());
    Vector<Double> incrSol(eqnData_->GetNumEQNs());

    Double* actRHS;
    Double* solPtr;

    Vector<Double> coeff;

    oldSol.Init(0);
    newSol.Init(0);
  
    //save solution opf previous time step  
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
    solhelp->GetAlgSysVector(solPrev);

    // Update extrema-list just for first time step
    if (laststepcalc_ == 1) {
      DoUpdateHyst();
    }

    //clear RHS
    algsys_->InitRHS();

    //set BCs
    SetBCs(lasttimecalc_);

    //compute constant part of RHS (div D(t) )
    // ComputeConstPartRHS();

    // stores this as linear part of RHS
    algsys_->GetRHSVal( actRHS );
    StoreAlgsysToVec(RhsLinVal_, actRHS );

    do {
      iterationCounter++;
      // for every time step write out number of iteration loops to standard out
      if (iterationCounter == 1)
        std::cout << std::endl << "Time step:   "  << kstep 
                  << "  ,Iterations: " << iterationCounter << std::endl;
      else 
        std::cout << "Iter:  " << iterationCounter << std::endl;
    
#ifdef DEBUG
      *debug << std::endl
             << "====================================================== "
             << std::endl
             <<   "Nonlinear Elecs: Perform internal loop no. "
             << iterationCounter << std::endl;      
#endif
        
      // set solution of previous iteration
      if (iterationCounter == 1 ) {
        oldSol = solPrev;
      }
      else {
        oldSol = newSol;
      }
    
      // Set linear part of RHS
      //    algsys_->InitRHS(RhsLinVal_.GetPointer());

      //clear RHS
      algsys_->InitRHS();

      //compute differential permittivity
      if ( iterationCounter != 1 || laststepcalc_ ==1 ) {
        ComputeDiffEpsilon();
      }

      //set up the new system-matrix
      algsys_->InitMatrix(SYSTEM);
      assemble_->SetReassemble();   

      assemble_->SetMaterialArray( &epsDiff_ ); 
      assemble_->AssembleMatrices();

      algsys_->ConstructEffectiveMatrix(matrix_factor_);

      // set changing part of RHS
      algsys_->UpdateRHS(SYSTEM,solPrev.GetPointer());

      //compute residual for incremental solution
      //    coeff = -oldSol;
      //algsys_->UpdateRHS(SYSTEM,coeff.GetPointer());

      //   algsys_->Print(SYSTEM);

      // build in the Dirichlet vales in system mmatrix and rhs
      algsys_->BuildInDirichlet();

      algsys_->Print(SYSTEM);

      //get RHS
      Vector<Double> RHS;
      algsys_->GetRHSVal( actRHS );
      StoreAlgsysToVec(RHS, actRHS );       
      Double residualNorm = RhsL2Norm( RHS );

      algsys_->SetupSolver();
      algsys_->SetupPrecond();
    
      algsys_->Solve();
      algsys_->GetSolutionVal( solPtr );
      StoreAlgsysToVec(incrSol, solPtr );

      //Double alpha = 1;
      //    newSol = oldSol + incrSol * alpha;
      newSol = incrSol;

      //    std::cout << "New Solution:\n" << newSol << std::endl;
      //    std::cout << "Old Solution:\n" << oldSol << std::endl;
#ifdef DEBUG
      *debug << std::endl
             << "New Solution:" << std::endl << newSol << std::endl;
#endif

      //put new solution to sol_
      sol_->SetAlgSysVector(newSol);  
    
      // compute L2-Norm of error between last incremental solution and
      // actual incremental solution
      Double solIncrL2Norm=0;
      for (UInt i=0; i<newSol.GetSize(); i++)
        solIncrL2Norm += (newSol[i]-oldSol[i])*(newSol[i]-oldSol[i]);
    
      solIncrL2Norm = sqrt(solIncrL2Norm);
      Double actSolL2Norm = newSol.NormL2();
    
      Double incrementalErr;
      if (actSolL2Norm > 1)
        incrementalErr = solIncrL2Norm / actSolL2Norm;
      else
        incrementalErr = solIncrL2Norm;
    
      // output of norms and data
      nonLinLogging_ = TRUE;
      if ( nonLinLogging_ == TRUE ) {
        Info->WriteNonLinIter(pdename_, iterationCounter, residualNorm,
                              incrementalErr);
      }

      //    std::cout << "ResNorm=" << residualNorm << "  incrNorm=" 
      //        << incrementalErr << std::endl;
    
      // boolean variable, holds condition if another iteration step
      //   is necessary
      incStopCrit_ = 1e-2;
      nonLinMaxIter_ = 200;
      performOneMoreStep =  ( (incrementalErr > incStopCrit_) || (residualNorm > 1e-5))
        && (residualNorm > 1e-10) ;
    
    } while ( performOneMoreStep && iterationCounter < nonLinMaxIter_ );  

    if ( iterationCounter >=  nonLinMaxIter_) {
      Error(" Number of nonlinear iterations too larger");
    }
  
  }


  void SolveStepElec::PostStepStatic(const UInt kstep, const Double asteptime)
  {
    ENTER_FCN( "SolveStepElec::PostStepStatic" );

    if (isIterCoupled_)
      (*iterCoupledCounter_)++;


#ifdef ADAPTGRID
    if (flags->CalcErrorMap_)
      {
        Double         totalErr;
        ElemStoreSol<Double>  Sol_Mesh;
        Vector<Double> solVec;
      
        ptError_=new SpaceErrorEstimator();

        ptError_->Init(this);
      
        sol_->TransformNodeSolution(Sol_Mesh,sol_,eqnData_,ptgrid_);

        sol_->GetCompleteVector(solVec);

        //  solVec.Resize(sol_.size());
        //       int i;
        //       for (i=0; i<sol_.size(); i++)
        //        solVec[i]=sol_[0][i];

        ptError_->CalcErrorMap(solVec,subdoms_,ptgrid_,errorMap_,totalErr);
      
        std::cout << " total error of calculation:: " << totalErr << std::endl;
        *data << errorMap_ << std::endl;
      }
#endif
  }




  void SolveStepElec::AddPolarizationToRHS() {
  
    ENTER_FCN( "SolveStepElec::AddPolarizationToRHS" );

    //we assume, that the actual solution is stored in sol_!
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    GradientFieldOp<Double> * FieldOp = new GradientFieldOp<Double>(ptgrid_, &PDE_, eqnData_,
                                                                    *solhelp, ELEC_POTENTIAL, 
                                                                    isaxi_);

    Vector<Double> LCoord, Efield;
    Double Ecomp, Pval, PvalReduced;
    UInt comp = Dir_ - 1;
    UInt pdeElem=1;

    Matrix<Double>     ptCoord;
    Vector<Double>     sol, solderiv1, solderiv2, rhs;
    BaseFE             *ptElem;
    StdVector<UInt> connect;
    StdVector<Integer>  connect_PDE;
  
    BaseForm * rhsInt = new  PiezoPolarizationInt(Dir_, PDE_.getPDE_dofspernode(), isaxi_);

    for (UInt actSD=0; actSD<subdoms_.GetSize(); actSD++) {
      StdVector<Elem*> elemssd;
      ptgrid_->GetVolElems(elemssd,subdoms_[actSD]);
    
      for (UInt iel=0; iel < elemssd.GetSize(); iel++) {

        //compute the electric field intensity
        elemssd[iel]->ptElem->GetCoordMidPoint(LCoord);
        FieldOp->CalcElemGradField( Efield, elemssd[iel], LCoord, 1);

        //get correct component of electric field for scalar Preisach model
        Ecomp = Efield[comp]; 

        //compute polarization
        Pval = hyst_->computeValue(Ecomp,pdeElem);
     
        //      Pval = Ecomp * 8.85500e-11;
        PvalReduced = (8.85500e-12 - 5.84000E-09)*Ecomp + Pval;

        //      std::cerr << Ecomp << "  " << Pval << "  " << PvalReduced << std::endl;

        ptElem  = elemssd[iel]->ptElem;
        connect = elemssd[iel]->connect;
        GetElemCoords(connect, ptCoord);
      
        rhsInt->SetElemPtr(ptElem);
        rhsInt->SetFactor(PvalReduced);
        rhsInt->CalcElemVector(ptCoord, rhs);
      
        //get equation numbers 
        eqnData_->Node2EQN(connect, connect_PDE);

        //assemble
        algsys_->SetElementRHS(&rhs[0],  pdeId1_, connect_PDE.GetPointer(), connect_PDE.GetSize());

        pdeElem++;
      }  
    }

  }


  void SolveStepElec::DoUpdateHyst() {
  
    ENTER_FCN( "SolveStepElec::DoUpdateHyst" );

    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    GradientFieldOp<Double> * FieldOp = new GradientFieldOp<Double>(ptgrid_, &PDE_, eqnData_,
                                                                    *solhelp, ELEC_POTENTIAL, 
                                                                    isaxi_);

    Vector<Double> LCoord, Efield;
    StdVector<Elem*> elemssd;
    UInt pdeElem=1;
    Double Ecomp, Pval;
    Integer comp = Dir_-1;

    // loop over all subdomains
    for (UInt isd=0; isd<subdoms_.GetSize(); isd++) {
      // get vector of Elem of subdomain with color: subdoms[isd]
      ptgrid_->GetVolElems(elemssd,subdoms_[isd]);
      
      // loop over elements of subdomain
      for (UInt iel=0; iel< elemssd.GetSize(); iel++)        {
        elemssd[iel]->ptElem->GetCoordMidPoint(LCoord);

        //compute electric field
        FieldOp->CalcElemGradField( Efield, elemssd[iel], LCoord, 1);

        //get correct component of electric field for scalar Preisach model
        //and invoke the update MinMaxList method
        Ecomp = Efield[comp]; 
        hyst_->updateMinMaxList(Ecomp, pdeElem);

        Pval = hyst_->computeValue(Ecomp,pdeElem);

        //      std::cerr << Ecomp << "    " << Pval << std::endl;

        Eprevious_[pdeElem-1] = Ecomp;
        Dprevious_[pdeElem-1] = Pval; //8.85419E-12 * Ecomp + Pval;

        pdeElem++;
      }
    }

    delete FieldOp;

  }

  void SolveStepElec::ComputeConstPartRHS() {
  
    ENTER_FCN( "SolveStepElec::ComputeConstPartRHS" );

    //we assume, that the actual solution is stored in sol_!
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    GradientFieldOp<Double> * FieldOp = new GradientFieldOp<Double>(ptgrid_, &PDE_, eqnData_,
                                                                    *solhelp, ELEC_POTENTIAL, 
                                                                    isaxi_);

    Vector<Double> LCoord, Efield;
    Double Ecomp, Pval, Dval;
    Integer comp = Dir_ - 1;
    UInt pdeElem=1;

    Matrix<Double>     ptCoord;
    Vector<Double>     Dfield, rhs;
    BaseFE             *ptElem;
    StdVector<UInt> connect;
    StdVector<Integer>   connect_PDE;
  
    ElecPolarizationInt *rhsInt = new   ElecPolarizationInt(isaxi_);

    for (UInt actSD=0; actSD<subdoms_.GetSize(); actSD++) {
      StdVector<Elem*> elemssd;
      ptgrid_->GetVolElems(elemssd,subdoms_[actSD]);
    
      for (UInt iel=0; iel < elemssd.GetSize(); iel++) {

        //compute the electric field intensity
        elemssd[iel]->ptElem->GetCoordMidPoint(LCoord);
        FieldOp->CalcElemGradField( Efield, elemssd[iel], LCoord, 1);

        //get correct component of electric field for scalar Preisach model
        Ecomp = Efield[comp]; //.NormL2(); //[comp]; 

        //compute polarization
        Pval = hyst_->computeValue(Ecomp,pdeElem);

        //      Pval = Ecomp*2e-7;

        //compute dielectric displacement
        Dval = Pval; //8.85419E-12 * Ecomp + Pval;

        if (Ecomp > 0 ) {
          Dfield  = Efield * ( Dval/Ecomp);
        }
        else {
          Dfield = Efield;
        }
        ptElem  = elemssd[iel]->ptElem;
        connect = elemssd[iel]->connect;
        GetElemCoords(connect, ptCoord);
      
        rhsInt->SetElemPtr(ptElem);
        rhsInt->SetSrcVec(Dfield);
        rhsInt->CalcElemVector(ptCoord, rhs);
      
        //get equation numbers 
        eqnData_->Node2EQN(connect, connect_PDE);

        //assemble
        algsys_->SetElementRHS(&rhs[0], pdeId1_, connect_PDE.GetPointer(), connect_PDE.GetSize());

        pdeElem++;
      }  
    }

  }

  void SolveStepElec::ComputeDiffEpsilon() {
  
    ENTER_FCN( "SolveStepElec::ComputeDiffEpsilon" );

    //we assume, that the actual solution is stored in sol_!
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    GradientFieldOp<Double> * FieldOp = new GradientFieldOp<Double>(ptgrid_, &PDE_, eqnData_,
                                                                    *solhelp, ELEC_POTENTIAL, 
                                                                    isaxi_);

    Vector<Double> LCoord, Efield;
    Double Ecomp, Pval, Dval, dE, dD, eps;
    Integer comp = Dir_ - 1;
    UInt pdeElem=1;

    for (UInt actSD=0; actSD<subdoms_.GetSize(); actSD++) {
      StdVector<Elem*> elemssd;
      ptgrid_->GetVolElems(elemssd,subdoms_[actSD]);
    
      for (UInt iel=0; iel < elemssd.GetSize(); iel++) {

        //compute the electric field intensity
        elemssd[iel]->ptElem->GetCoordMidPoint(LCoord);
        FieldOp->CalcElemGradField( Efield, elemssd[iel], LCoord, 1);

        //get correct component of electric field for scalar Preisach model
        Ecomp = Efield[comp]; //.NormL2(); //[comp]; 

        //compute polarization
        Pval = hyst_->computeValue(Ecomp,pdeElem);

        //      Pval = Ecomp*2e-7;

        //compute dielectric displacement
        Dval = Pval; //8.85419E-12 * Ecomp + Pval;

        //compute differential epsilon
        dE = Ecomp - Eprevious_[pdeElem-1];
        dD = Dval - Dprevious_[pdeElem-1];
        if ( (abs(dD) < 1e-12) || (abs(dE) < 1e-10) ) {
          eps = materialData_[actSD].GetPermittivity(2,2);
          if (eps < 8.854e-12) {
            eps = 8.854e-12;
          }
          epsDiff_[actSD][iel] = eps;
        }
        else {
          epsDiff_[actSD][iel] = dD / dE;
        }

        pdeElem++;
      }  
    }

    //  std::cout << "EpsDiff: " << epsDiff_ << std::endl;
  }


  void SolveStepElec::StepStaticNonLin(const UInt kstep, const Double asteptime,
                                       const Boolean reset) {

    ENTER_FCN( "SolveStepElec::StepStaticNonLin" );

    //   laststepcalc_ = kstep;
    //   lasttimecalc_ = asteptime;

    //   UInt job;
    //   Boolean performOneMoreStep;
    //   UInt iterationCounter=0;
  
    //   Vector<Double> newSol(numPDENodes_);
    //   Vector<Double> oldSol(numPDENodes_);
  
    //   // just update dirichlet values
    //   job = 3;
  
    //   // if first time step, setup system matrix
    //   if (laststepcalc_ == 1) {
    //     assemble_->AssembleMatrices();
    //     algsys_->ConstructEffectiveMatrix(matrix_factor_);
    
    //     //set job to 1: build in dirichlet BCs and compute preconditioner
    //     job = 1;
    //   }  
  
    //   NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    //   // set BCs, if effective mass matrix formulation, values of BCs depend on 
    //   //  predictors, so predictors have to be computed beforehand
    //   SetBCs(lasttimecalc_);

    //   // set old solution  
    //   newSol = solhelp->GetAlgSysVector();

    //   // stores this as linear part of RHS
    //   // StoreAlgsysToVec(RhsLinVal_, algsys_->GetRHSVal() );

    //   // Update extrema-list just for first time step
    //   if (laststepcalc_ == 1) {
    //     DoUpdateHyst();
    //   }

    //   do {
    //     iterationCounter++;
    //     // for every time step write out number of iteration loops to standard out
    //     if (iterationCounter == 1)
    //       std::cout << std::endl << "Time step:   "  << kstep 
    //              << "  ,Iterations: " << iterationCounter << std::endl;
    //     else 
    //       std::cout      << "  " << iterationCounter << std::endl;
    
    // #ifdef DEBUG
    //     *debug << std::endl
    //         << "====================================================== "
    //         << std::endl
    //         <<   "Nonlinear Elecs: Perform internal loop no. "
    //         << iterationCounter << std::endl;      
    // #endif
        
    //     // set solution of previous iteration
    //     oldSol = newSol;
    
    //     if ( iterationCounter>1 ) 
    //       job = 3;
    
    //     // Set linear part of RHS
    //     //    algsys_->InitRHS(RhsLinVal_.GetPointer());
    //     algsys_->InitRHS();
    
    //     //std::cout << "  ,Iterations: " << iterationCounter << std::endl;
    //     // put nonlinear part to RHS
    //     AddPolarizationToRHS();

    //     // calculation of residual error (takes care for Dirichlet BCs========
    //     Vector<Double> actRHS;
    //     StoreAlgsysToVec(actRHS, algsys_->GetRHSVal() );       
  
    //     // std::cerr << "RHS\n" << actRHS << std::endl;

    //     algsys_->BuildInDirichlet();

    //     StoreAlgsysToVec(actRHS, algsys_->GetRHSVal() ); 
    //     //  std::cerr << "RHS\n" << actRHS << std::endl;

    //     if ( job == 1 ) {
    //       algsys_->SetupSolver(job);
    //       algsys_->SetupPrecond(job);
    //     }
    
    //     algsys_->Solve();
    
    //     // store new solution in newSol
    //     StoreAlgsysToVec(newSol, algsys_->GetSolutionVal() );

    //     //   Double relax_val = 0.1;
    //     //    newSol = oldSol*(1-relax_val) +  newSol*relax_val;
     
    //     //   std::cout << "NewSol:\n " << newSol << std::endl;

    // #ifdef DEBUG
    //     *debug << std::endl
    //         << "New Solution:" << std::endl << newSol << std::endl;
    // #endif

    //     //put new solution to sol_
    //     sol_->SetAlgSysVector(newSol);  
    
    //     // compute L2-Norm of error between last incremental solution and
    //     // actual incremental solution
    //     Double solIncrL2Norm=0;
    //     for (UInt i=0; i<newSol.GetSize(); i++)
    //       solIncrL2Norm += (newSol[i]-oldSol[i])*(newSol[i]-oldSol[i]);
    
    //     solIncrL2Norm = sqrt(solIncrL2Norm);
    //     Double actSolL2Norm = newSol.NormL2();
    
    //     Double incrementalErr;
    //     if (actSolL2Norm > 1)
    //       incrementalErr = solIncrL2Norm / actSolL2Norm;
    //     else
    //       incrementalErr = solIncrL2Norm;
    
    //     // output of norms and data
    //     nonLinLogging_ = TRUE;
    //     if ( nonLinLogging_ == TRUE ) {
    //       Info->WriteNonLinIter(pdename_, iterationCounter, incrementalErr,
    //                          incrementalErr);
    //     }
    
    //     // boolean variable, holds condition if another iteration step
    //     //   is necessary
    //     incStopCrit_ = 1e-4;
    //     nonLinMaxIter_ = 200;
    //     performOneMoreStep =    (incrementalErr > incStopCrit_);
    
    //   } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  
  
  
  }

} // end of namespace

