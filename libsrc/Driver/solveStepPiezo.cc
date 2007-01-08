#include <fstream>
#include <iostream>
#include <string>

#include "solveStepPiezo.hh"
#include "Utils/preisach.hh"
#include "Utils/jiles.hh"
#include "assemble.hh"
#include "Forms/gradfieldop.hh"
#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"
#include "Domain/domain.hh"

namespace CoupledField {

  SolveStepPiezo::SolveStepPiezo(StdPDE& apde) : StdSolveStep(apde)
  {

    ENTER_FCN( "SolveStepPiezo::SolveStepPiezo" );
    doInit_ = true;
  }


  SolveStepPiezo::~SolveStepPiezo() {
    ENTER_FCN( "SolveStepPiezo::~SolveStepPiezo" );
  }
 

  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  void SolveStepPiezo:: PreStepTrans()
  {
    ENTER_FCN( "SolveStepPiezo::PreStepStatic" );


    // due to coupling-pdes, the RHS has to be initialized BEFORE 
    // the coupling forces are assembled to the RHS
    algsys_->InitRHS();

    UInt numElems = PDE_.getPDE_numElems();

    if (isHyst_) {
      if (actStep_==1) {
        Eprevious_.Resize(numElems);
        Dprevious_.Resize(numElems);
        epsDiff_.Resize(subdoms_.GetSize(),numElems);
        Eprevious_.Init(0);
        Dprevious_.Init(0);
        epsDiff_.Init(8.854e-12);
      }

      if (doInit_) {
        //set the Preisach values; should be obtained from the xml_file
        StdVector<Elem*> elemssd;
        bool isVirgin;
        hyst_.Resize(subdoms_.GetSize());
        hyst_.Init(NULL);

        for (UInt iSD=0; iSD<subdoms_.GetSize(); iSD++) {
          Double Esat, Psat;
          Double Ec = 0; //currently not used
          UInt dir, numSDElems;

          ptgrid_->GetVolElems(elemssd,subdoms_[iSD]);
          numSDElems = elemssd.GetSize(); 

          std::string hystType;
          materialData_[iSD]->GetScalar(hystType, HYST_MODEL);
          if ( hystType == "preisach" ) {
            materialData_[iSD]->GetScalar(Esat,E_SATURATION,REAL);
            materialData_[iSD]->GetScalar(Psat,P_SATURATION,REAL);
            materialData_[iSD]->GetScalar((Integer&)dir,P_DIRECTION,INTEGER);
            isVirgin = true; 

	    //            hyst_[iSD] = new Preisach(numSDElems, Esat, Psat, isVirgin);
          }
          else if (hystType == "jiles") {
            Double a, alpha, k, c;
	    materialData_[iSD]->GetScalar(a,A_JILES,REAL);
	    materialData_[iSD]->GetScalar(alpha,ALPHA_JILES,REAL);
	    materialData_[iSD]->GetScalar(k,K_JILES,REAL);
	    materialData_[iSD]->GetScalar(c,C_JILES,REAL);
	    materialData_[iSD]->GetScalar(Psat,P_SATURATION,REAL);

            hyst_[iSD] = new Jiles(numSDElems, Psat, a, alpha, k, c);
            hyst_[iSD]->SetTimeStepVal(TS_alg_->GetTimeStep());
          }
        }
        doInit_ = false;
      }
      else {
        //   std::cout << "Do Update in PreStep" << std::endl;
        DoUpdateHyst();
      }
    }
  }



  // time is used for a series of static calculations
  // don't get confused with REAL transient simulations!
  void SolveStepPiezo::SolveStepTrans() {

    ENTER_FCN( "SolveStepPiezo::SolveStepTrans" );

    std::cout<<"Solve Step Trans"<<std::endl;  
    if (isHyst_) {
      StepTransNonLinEpsDiff();
    }
    else if (nonLin_){
      std::cout<<"We do a transient nonlinear calculation"<<std::endl;
      StepTransNonLin();
    }
    else if (nonLinMaterial_){
      std::cout<<"We do a transient nonlinear calculation"<<std::endl;
      StepTransNonLinMaterial();
    }

    else {
      std::cout<<"We do a transient linear calculation"<<std::endl;
      StepTransLin();
    }

  }

  void SolveStepPiezo::StepTransMaterialNonLin(){
    
    ENTER_FCN( "SolveStepPiezo::StepTransMaterialNonLin" );

      std::cout<<"We do a step of transient nonlinear calculation"<<std::endl;


  }

  void SolveStepPiezo::StepTransNonLinEpsDiff() {

    ENTER_FCN( "SolveStepPiezo::StepTransNonLinEpsDiff" );

    bool performOneMoreStep;
    UInt iterationCounter=0;
  
    UInt numEqns = eqnMap_->GetNumEqns();
    Vector<Double> newSol(numEqns);
    Vector<Double> oldSol(numEqns);
    Vector<Double> solPrev(numEqns);
    Vector<Double> incrSol(numEqns);

    Double* actRHS;
    Double* solPtr;

    Vector<Double> coeff;

    oldSol.Init(0);
    newSol.Init(0);
  
    //save solution of previous time step  
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
    solhelp->GetAlgSysVector(solPrev);

    // Update extrema-list just for first time step
    if (actStep_ == 1) {
      DoUpdateHyst();
    }

    //clear RHS
    algsys_->InitRHS();

    //set BCs
    PDE_.SetBCs(actTime_);

    // stores this as linear part of RHS
    algsys_->GetRHSVal( actRHS );
    StoreAlgsysToVec(RhsLinVal_, actRHS );

    do {
      iterationCounter++;
      // for every time step write out number of iteration loops to standard out
      if (iterationCounter == 1)
        std::cout << std::endl << "Time step:   "  << actStep_ 
                  << "  ,Iterations: " << iterationCounter << std::endl;
      else 
        std::cout << "Iter:  " << iterationCounter << std::endl;
    
#ifdef DEBUG
      *debug << std::endl
             << "====================================================== "
             << std::endl
             <<   "Nonlinear Piezos: Perform internal loop no. "
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
      if ( iterationCounter != 1 || actStep_ ==1 ) {
        ComputeDiffEpsilon();
      }

      //set up the new matrices
      algsys_->InitMatrix();
      //assemble_->SetReassemble();   

      //assemble_->SetMaterialArray( &epsDiff_ ); 
      Warning( "Assemble::SetMaterialArray not implemented!",
               __FILE__, __LINE__ );
      assemble_->AssembleMatrices();

      algsys_->ConstructEffectiveMatrix(matrix_factor_);

      // set changing part of RHS
      algsys_->UpdateRHS(SYSTEM,solPrev.GetPointer());

      //compute residual for incremental solution
      //    coeff = -oldSol;
      //algsys_->UpdateRHS(SYSTEM,coeff.GetPointer());

      // build in the Dirichlet vales in system mmatrix and rhs
      algsys_->BuildInDirichlet();

      //   algsys_->Print(SYSTEM);

      //get RHS
      Vector<Double> RHS;
      algsys_->GetRHSVal( actRHS );
      StoreAlgsysToVec(RHS, actRHS );       
      Double residualNorm = PDE_.GetRhsL2Norm( RHS );

      algsys_->SetupSolver();
      algsys_->SetupPrecond();
    
      algsys_->Solve();
      algsys_->GetSolutionVal( solPtr );
      StoreAlgsysToVec(incrSol, solPtr );

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
      if ( nonLinLogging_ == true ) {
        Info->WriteNonLinIter(pdename_, iterationCounter, residualNorm,
                              incrementalErr);
      }

      // boolean variable, holds condition if another iteration step
      //   is necessary
      performOneMoreStep =  ( (incrementalErr > incStopCrit_) || (residualNorm > residualStopCrit_) )
        && (residualNorm > 1e-14 || iterationCounter == 1 ) ;
                              
    } while ( performOneMoreStep && iterationCounter < nonLinMaxIter_ );  

    if ( iterationCounter >=  nonLinMaxIter_) {
      Error(" Number of nonlinear iterations too larger", __FILE__, __LINE__ );
    }
  
    //  std::cout << "Do Update after Sol" << std::endl;
    DoUpdateHyst();
  }


  void SolveStepPiezo::DoUpdateHyst() {
  
    ENTER_FCN( "SolveStepPiezo::DoUpdateHyst" );

    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    // Get correct result-type
    GradientFieldOp<Double> * FieldOp = new GradientFieldOp<Double>(ptgrid_, &PDE_,
                                                                    eqnMap_,
                                                                    *solhelp, ELEC_POTENTIAL, 
                                                                    results_[0], isaxi_);

    Vector<Double> LCoord, Efield;
    StdVector<Elem*> elemssd;
    UInt pdeElem=1;
    Double Ecomp, Pval;
    UInt comp;

    // loop over all subdomains
    for (UInt isd=0; isd<subdoms_.GetSize(); isd++) {

      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion( subdoms_[isd] );
      EntityIterator it = actSDList.GetIterator();
      
      //get direction of polarization
      materialData_[isd]->GetScalar((Integer&)comp,P_DIRECTION,INTEGER);
      comp -= 1;

      // loop over elements of subdomain
      for ( it.Begin(); !it.IsEnd(); it++ ) {
        it.GetElem()->ptElem->GetCoordMidPoint(LCoord);

        //compute electric field
        FieldOp->CalcElemGradField( Efield, it, LCoord, 1);

        //get correct component of electric field for scalar Preisach model
        //and invoke the update MinMaxList method
        Ecomp = Efield[comp]; 
        //      pdeElem = 0;
        hyst_[isd]->updateMinMaxList(Ecomp, pdeElem);

        Pval = hyst_[isd]->computeValue(Ecomp,pdeElem);

        //      std::cerr << Ecomp << "    " << Pval << std::endl;

        Eprevious_[pdeElem-1] = Ecomp;
        Dprevious_[pdeElem-1] = Pval; //8.85419E-12 * Ecomp + Pval;

        pdeElem++;
      }
    }

    delete FieldOp;

  }


  void SolveStepPiezo::ComputeDiffEpsilon() {
  
    ENTER_FCN( "SolveStepPiezo::ComputeDiffEpsilon" );

    //we assume, that the actual solution is stored in sol_!
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    
    GradientFieldOp<Double> * FieldOp = 
      new GradientFieldOp<Double>(ptgrid_, &PDE_, eqnMap_,
                                  *solhelp, ELEC_POTENTIAL, 
                                  results_[0],isaxi_);

    Vector<Double> LCoord, Efield;
    Double Ecomp, Pval, Dval, dE, dD, eps;
    UInt comp;
    UInt pdeElem=1;

    for (UInt actSD=0; actSD<subdoms_.GetSize(); actSD++) {
      ElemList actSDList(domain->GetGrid() );
      actSDList.SetRegion( subdoms_[actSD] );
      EntityIterator it = actSDList.GetIterator();
      
      //get direction of polarization
      materialData_[actSD]->GetScalar((Integer&)comp,P_DIRECTION,INTEGER);
      comp -= 1;
      UInt iel = 0;
      for ( it.Begin(); !it.IsEnd(); it++, iel++ ) {
        
        //compute the electric field intensity
        it.GetElem()->ptElem->GetCoordMidPoint(LCoord);
        FieldOp->CalcElemGradField( Efield, it, LCoord, 1);

        //get correct component of electric field for scalar Preisach model
        Ecomp = Efield[comp]; //.NormL2(); //[comp]; 

        //compute polarization
        //      pdeElem = 0;
        Pval = hyst_[actSD]->computeValue(Ecomp,pdeElem);

        //      Pval = Ecomp*2e-7;

        //compute dielectric displacement
        Dval = Pval; //8.85419E-12 * Ecomp + Pval;

        //compute differential epsilon
        dE = Ecomp - Eprevious_[pdeElem-1];
        dD = Dval - Dprevious_[pdeElem-1];
        if ( (abs(dD) < 1e-12) || (abs(dE) < 1e-10) ) {
          materialData_[actSD]->GetScalar(eps,ELEC_PERMITTIVITY,REAL);
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



} // end of namespace

