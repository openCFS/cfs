// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_use_mpcci.hh>

#include "iterSolveStep.hh"

#include "MatVec/basematrix.hh"
#include "PDE/StdPDE.hh"
#include "CoupledPDE/itercoupledpde.hh"
#include "CoupledPDE/pdecoupling.hh"

namespace CoupledField
{
  //! Derived class for step-wise solving of iterative coupled StdPDEs

  IterSolveStep::IterSolveStep(IterCoupledPDE &apde) 
    : BaseSolveStep(),
      rPDE_(apde),
      rCouplings_(apde.Couplings_)
  {
    
    startStep_ = 1;

  }

  IterSolveStep::~IterSolveStep()
  {
  }


  //----------------------- STATIC--------------------------------------------

  void IterSolveStep::SolveStepStatic(PtrParamNode analysis_id, const bool reAssembleMatrices)
  {
  
    SingleVector *val, *oldVal;
    UInt iter = 0;
    UInt counter = 0;
    bool normsReached = false;


    while (iter < rPDE_.maxiter_ &&  (! normsReached)) {

      if (rPDE_.nonLinLogging_) {
        
        Info->PrintF(rPDE_.pdename_,"\n"); 
        Info->PrintF(rPDE_.pdename_, 
                     " COUPLED ITERATION %i =================================\n",
                     iter+1);
      }
      
      counter = 0;
      normsReached = true;
      
      for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
        
        if (rPDE_.nonLinLogging_)
          Info->PrintF(rPDE_.pdename_, " Processing PDE %s\n", 
                       (rPDE_.PDEs_[i]->GetName()).c_str());
        
        rPDE_.PDEs_[i]->GetSolveStep()->SetActTime(actTime_);
        rPDE_.PDEs_[i]->GetSolveStep()->SetActStep(actStep_);
        rPDE_.PDEs_[i]->GetSolveStep()->PreStepStatic();
        rPDE_.PDEs_[i]->CalcInputCoupling();
        rPDE_.PDEs_[i]->GetSolveStep()->SolveStepStatic(analysis_id);
        rPDE_.PDEs_[i]->GetSolveStep()->PostStepStatic();
        rPDE_.PDEs_[i]->CalcOutputCoupling();
        
        // Calculate Norms
        for (UInt k=0; k<rCouplings_[i]->GetNumOutputCouplings(); k++) {
          
          rCouplings_[i]->GetOutputValues(k, val);
          rCouplings_[i]->GetOutputOldValues(k, oldVal);
          
          rPDE_.norms_[counter] = 
            CalcNorm(rCouplings_[i]->GetOutputNormType(k), *val, *oldVal);
          
          if (rPDE_.nonLinLogging_) {
            
            Info->PrintF(rPDE_.pdename_, " %s : Norm of %s = %g\n", 
                         (rCouplings_[i]->GetPDE()->GetName()).c_str(),
                         (SolutionTypeEnum.ToString(rCouplings_[i]->GetOutputQuantity(k))).c_str(),
                         rPDE_.norms_[counter]);
          }
          
          if (rPDE_.norms_[counter] > rCouplings_[i]->GetOutputEpsilon(k) && 
              rCouplings_[i]->GetOutputNormType(k) != NO_NORM)
            normsReached = false;
          
          //copy values of new solution to old one
          dynamic_cast<Vector<Double>&>(*oldVal) = 
            dynamic_cast<Vector<Double>&>(*val);
          
        }
        counter++;            
      } // end of for-loop
      
      iter++;
      
      if(rPDE_.nonLinLogging_)
        Info->PrintF(rPDE_.pdename_, "\n");
      
    } // end of while-loop
    
  }


  //----------------------- TRANSIENT-----------------------------------------
  void IterSolveStep::SolveStepTrans(PtrParamNode analysis_id)
  {

    UInt iter = 0;

    bool normsReached = false;
    
    // In the beginning of each time step
    // the coupling data has to be reseted
    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
      
      rPDE_.PDEs_[i]->ResetCoupling();
      rPDE_.PDEs_[i]->GetSolveStep()->SetStartStep(startStep_);
    }

    //in case of FSI the predictor of the mechanic time integration is called first
//    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
//      rPDE_.PDEs_[i]->GetSolveStep()->PredictorStep();
//    }

   
    while (iter < rPDE_.maxiter_ &&  (! normsReached)) {
      
      if (rPDE_.nonLinLogging_) {
        
        Info->PrintF(rPDE_.pdename_,"\n"); 
        Info->PrintF(rPDE_.pdename_,
                     " COUPLED ITERATION %i =================================\n", 
                     iter+1);
      }

      UInt counter = 0;
      normsReached = true;

#ifdef MpCCI
      // whenever iter == 0 the old time step has converged
      // in CalcInputCoupling of mpcciPDE CFS++ tells with the 
      // flag converged_==true that a new time step begins
      //For the first iteration of the first time step this is ignored 
      //because of the flag flagFirstTime in CalcInputCoupling of mpcciPDE
      if (iter == 0)
        rPDE_.PDEs_[0]->converged_ = true; 
      else
        rPDE_.PDEs_[0]->converged_ = false; 
#endif

      for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
        
        
        if (rPDE_.nonLinLogging_)
          Info->PrintF(rPDE_.pdename_, " Processing PDE %s\n", 
                       (rPDE_.PDEs_[i]->GetName()).c_str());

       
          
#ifdef MpCCI
        // check for a HALTCFS File
        // if there exist a file with name HALTCFS in the executing directory
        // than CFS++ will create a HALT file such that FASTEST also stops
        std::ifstream readHALTCFS("HALTCFS", std::ios_base::in );
        if (readHALTCFS && actStep_== numTimeStep_) {
          readHALTCFS.close();
          std::ofstream stopFASTEST3D("../HALT", std::ios_base::out );
        }
#endif

        rPDE_.PDEs_[i]->GetSolveStep()->SetActTime(actTime_);
        rPDE_.PDEs_[i]->GetSolveStep()->SetActStep(actStep_);
        rPDE_.PDEs_[i]->GetSolveStep()->PreStepTrans();
        rPDE_.PDEs_[i]->CalcInputCoupling();
        rPDE_.PDEs_[i]->GetSolveStep()->SolveStepTrans(analysis_id);

        rPDE_.PDEs_[i]->CalcOutputCoupling();
              
        // Calculate Norms
        for (UInt k=0; k<rCouplings_[i]->GetNumOutputCouplings(); k++) {
            
          SingleVector *val, *oldVal;
          rCouplings_[i]->GetOutputValues(k, val);
          rCouplings_[i]->GetOutputOldValues(k, oldVal);
          rPDE_.norms_[counter] = 
            CalcNorm(rCouplings_[i]->GetOutputNormType(k), *val, *oldVal);
            
          if (rPDE_.nonLinLogging_) {
              
            Info->PrintF(rPDE_.pdename_, " %s : Norm of %s = %g\n", 
                         (rCouplings_[i]->GetPDE()->GetName()).c_str(),
                         (SolutionTypeEnum.ToString(rCouplings_[i]->GetOutputQuantity(k))).c_str(), 
                         rPDE_.norms_[counter]);

            Info->PrintF(rPDE_.pdename_, " actStep_ = %d\n", actStep_);
            Info->PrintF(rPDE_.pdename_, " numTimeStep_ = %d\n", numTimeStep_);
          }
          if (rPDE_.norms_[counter] > rCouplings_[i]->GetOutputEpsilon(k)) 
            normsReached = false;
                  
          dynamic_cast<Vector<Double>&>(*oldVal) = 
            dynamic_cast<Vector<Double>&>(*val);
          //*oldVal = *val;
          counter++;              
        }
      } // end of for-loop

      iter++;
      if (rPDE_.nonLinLogging_)
        Info->PrintF(rPDE_.pdename_, "\n"); 

    } // end of while-loop

    // now we are converged and can compute any postprocessing-quantities
    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++)
      rPDE_.PDEs_[i]->GetSolveStep()->PostStepTrans();


#ifdef MpCCI
    if (actStep_== numTimeStep_) {
      rPDE_.PDEs_[0]->converged_ = true;
      rPDE_.PDEs_[0]->CalcInputCoupling();
    }
#endif
    
  } 

  //----------------------- HARMONIC---------------------------------------
  void IterSolveStep::SolveStepHarmonic(PtrParamNode analysis_id)
  {
    EXCEPTION("Harmonic iterative coupling is not yet implemented"); 
  }


  void IterSolveStep::SetActTime( const Double actTime )
  {

    actTime_ = actTime;

    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
      actAnalysisType_ = rPDE_.PDEs_[i]->GetAnalysisType();

      if ( actAnalysisType_ == BasePDE::TRANSIENT )
        rPDE_.PDEs_[i]->GetSolveStep()->SetActTime(actTime);
    }
  }

  void IterSolveStep::SetActFreq( const Double actFreq )
  {

    actFreq_ = actFreq;
    
    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
      actAnalysisType_ = rPDE_.PDEs_[i]->GetAnalysisType();

      if ( actAnalysisType_ == BasePDE::HARMONIC )
        rPDE_.PDEs_[i]->GetSolveStep()->SetActFreq(actFreq);
    }
  }

  void IterSolveStep::SetActStep( const UInt actStep )
  {
    
    actStep_ = actStep;

    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
      actAnalysisType_ = rPDE_.PDEs_[i]->GetAnalysisType();

      if ( actAnalysisType_ == BasePDE::TRANSIENT )
        rPDE_.PDEs_[i]->GetSolveStep()->SetActStep(actStep);
      // Dirty Hack!!
      else if ( actAnalysisType_ == BasePDE::HARMONIC )
        rPDE_.PDEs_[i]->GetSolveStep()->SetActStep(1);
    }
  }


  void IterSolveStep::SetNumTimeSteps( UInt numTimeStep)
  {

    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
      actAnalysisType_ = rPDE_.PDEs_[i]->GetAnalysisType();

      if ( actAnalysisType_ == BasePDE::TRANSIENT ) {
        rPDE_.PDEs_[i]->GetSolveStep()->SetNumTimeSteps(numTimeStep);
        numTimeStep_=numTimeStep;
      }
    }
  }

  void IterSolveStep::SetStartStep( const UInt startStep )
  {
    
    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
      actAnalysisType_ = rPDE_.PDEs_[i]->GetAnalysisType();

      if ( actAnalysisType_ == BasePDE::TRANSIENT )
        rPDE_.PDEs_[i]->GetSolveStep()->SetStartStep(startStep);
    }
  }


  Double IterSolveStep::CalcNorm(NormType normtype, SingleVector & val,
                                 SingleVector & oldval)
  {

    // ATTENTION: Currently only working with Double-values
    // will be changed as soon as dynamic type information
    // is available

    SingleVector * delta = NULL;
 
		// initialize or receive compiler warning
    Double norm = 0.0, valNorm2 = 0.0;
  
    // Distinguish complex and real case
    if ( val.GetEntryType() == BaseMatrix::COMPLEX ) {

      delta = new Vector<Complex>;
      delta->Resize( val.GetSize() );
      delta->Init( );
    } else {
      delta = new Vector<Double>;
      delta->Resize( val.GetSize() );
      delta->Init();
    }
   
    
    // Calculate difference
    delta->Add(1.0, val, -1.0, oldval );

//    Vector<Double> & val_vec = dynamic_cast<Vector<Double>& >(val);
//    Vector<Double> & oldval_vec = dynamic_cast<Vector<Double>& >(oldval);
//    delta = val_vec - oldval_vec;

    switch (normtype)
      {
      case NO_NORM:
        return 0;
        break;
      
      case L2ABS:
        norm = delta->NormL2();
        break;

      case L2REL:
        valNorm2 =  val.NormL2();
        if (valNorm2 > 0)
          norm = delta->NormL2() / valNorm2;
        else
          norm = delta->NormL2();

        break;
      }

    // Delete delta-vector
    delete delta;

    return norm;
  }
} // end of namespace
