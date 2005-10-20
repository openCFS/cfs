#include "iterSolveStep.hh"

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
    ENTER_FCN( "IterSolveStep::IterSolveStep" );
    
    startStep_ = 1;

  }

  IterSolveStep::~IterSolveStep() {
    
    ENTER_FCN( "IterSolveStep::~IterSolveStep" );

  }


  //----------------------- STATIC---------------------------------------

  void IterSolveStep::SolveStepStatic( const Boolean updatesysmat ) {

    ENTER_FCN ( "IterSolveStep::SolveStepStatic" );
  
    CFSVector *val, *oldVal;
    UInt iter = 0;
    UInt counter = 0;
    Boolean normsReached = FALSE;
    std::string quantityConv;


    while (iter < rPDE_.maxiter_ &&  (! normsReached))
      {
        if (rPDE_.nonLinLogging_)
          {
            Info->PrintF(rPDE_.pdename_,"\n"); 
            Info->PrintF(rPDE_.pdename_, " COUPLED ITERATION %i =================================\n", 
                         iter+1);
          }
        
        counter = 0;
        normsReached = TRUE;
      
        for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++)
          {
            if (rPDE_.nonLinLogging_)
              Info->PrintF(rPDE_.pdename_, " Processing PDE %s\n", 
                           (rPDE_.PDEs_[i]->GetName()).c_str());

            // Only solve current PDE, if the corresponding
            // flag in 'solvePDE_' is set to TRUE
            if (rPDE_.solvePDE_[i] == TRUE) {
              rPDE_.PDEs_[i]->GetSolveStep()->SetActTime(actTime_);
              rPDE_.PDEs_[i]->GetSolveStep()->SetActStep(actStep_);
              rPDE_.PDEs_[i]->GetSolveStep()->PreStepStatic(updatesysmat);
              rPDE_.PDEs_[i]->CalcInputCoupling();
              rPDE_.PDEs_[i]->GetSolveStep()->SolveStepStatic(updatesysmat);
              rPDE_.PDEs_[i]->GetSolveStep()->PostStepStatic();
              rPDE_.PDEs_[i]->CalcOutputCoupling();

              // Calculate Norms
              for (UInt k=0; k<rCouplings_[i]->GetNumOutputCouplings(); k++)
                {
                  rCouplings_[i]->GetOutputValues(k, val);
                  rCouplings_[i]->GetOutputOldValues(k, oldVal);
                  
                  rPDE_.norms_[counter] = CalcNorm(rCouplings_[i]->GetOutputNormType(k), *val, *oldVal);

                  if (rPDE_.nonLinLogging_)
                    {
                      Enum2String(rCouplings_[i]->GetOutputQuantity(k), quantityConv);
                      Info->PrintF(rPDE_.pdename_, " %s : Norm of %s = %g\n", 
                                   (rCouplings_[i]->GetPDE()->GetName()).c_str(),
                                    quantityConv.c_str(), rPDE_.norms_[counter]);
                    }
                  
                  if (rPDE_.norms_[counter] > rCouplings_[i]->GetOutputEpsilon(k) && 
                      rCouplings_[i]->GetOutputNormType(k) != NO_NORM)
                    normsReached = FALSE;
                
                  //copy values of new solution to old one
                  *oldVal = *val;

                } // end of if 
              counter++;            
            }
          }

        iter++;
        
        if(rPDE_.nonLinLogging_)
          Info->PrintF(rPDE_.pdename_, "\n"); 
      }

    // now we are converged and can compute any postprocessing-quantities
    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++)
      rPDE_.PDEs_[i]->PostProcess();

  }


  //----------------------- TRANSIENT---------------------------------------
  void IterSolveStep::SolveStepTrans( const Boolean updatesysmat )
  {
    ENTER_FCN( "IterSolveStep::SolveStepTrans" );

    UInt iter = 0;

    Boolean normsReached = FALSE;
    std::string quantityConv;
    
    // In the beginning of each time step
    // the coupling data has to be reseted
    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++)
      {
        rPDE_.PDEs_[i]->ResetCoupling();
        rPDE_.PDEs_[i]->GetSolveStep()->SetStartStep(startStep_);
      }

   
    while (iter < rPDE_.maxiter_ &&  (! normsReached))
      {
        if (rPDE_.nonLinLogging_)
          {
            Info->PrintF(rPDE_.pdename_,"\n"); 
            Info->PrintF(rPDE_.pdename_, " COUPLED ITERATION %i =================================\n", 
                         iter+1);
          }

        UInt counter = 0;
        normsReached = TRUE;

#ifdef MpCCI
	// whenever iter == 0 the old time step has converged
	// in CalcInputCoupling of mpcciPDE CFS++ tells with the 
	// flag converged_==true that a new time step begins
	//For the first iteration of the first time step this is ignored 
	//because of the flag flagFirstTime in CalcInputCoupling of mpcciPDE
	if (iter == 0)
	  rPDE_.PDEs_[0]->converged_ = TRUE; 
	else
	  rPDE_.PDEs_[0]->converged_ = FALSE; 
#endif

      
        for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++)
          {


            if (rPDE_.nonLinLogging_)
              Info->PrintF(rPDE_.pdename_, " Processing PDE %s\n", 
                           (rPDE_.PDEs_[i]->GetName()).c_str());

            // Only solve current PDE, if the corresponding
            // flag in 'solvePDE_' is set to TRUE
            if (rPDE_.solvePDE_[i] == TRUE) {
              rPDE_.PDEs_[i]->GetSolveStep()->SetActTime(actTime_);
              rPDE_.PDEs_[i]->GetSolveStep()->SetActStep(actStep_);
              rPDE_.PDEs_[i]->GetSolveStep()->PreStepTrans( updatesysmat );
              rPDE_.PDEs_[i]->CalcInputCoupling();
              rPDE_.PDEs_[i]->GetSolveStep()->SolveStepTrans( updatesysmat );
              rPDE_.PDEs_[i]->GetSolveStep()->PostStepTrans();
              rPDE_.PDEs_[i]->CalcOutputCoupling();
              
              // Calculate Norms
              for (UInt k=0; k<rCouplings_[i]->GetNumOutputCouplings(); k++)
                {
                  CFSVector *val, *oldVal;
                  rCouplings_[i]->GetOutputValues(k, val);
                  rCouplings_[i]->GetOutputOldValues(k, oldVal);
                  rPDE_.norms_[counter] = CalcNorm(rCouplings_[i]->GetOutputNormType(k), *val, *oldVal);

                  if (rPDE_.nonLinLogging_)
                    {
                      Enum2String(rCouplings_[i]->GetOutputQuantity(k), quantityConv);
                      Info->PrintF(rPDE_.pdename_, " %s : Norm of %s = %g\n", 
                                   (rCouplings_[i]->GetPDE()->GetName()).c_str(),
                                   quantityConv.c_str(), rPDE_.norms_[counter]);
                    }
                  if (rPDE_.norms_[counter] > rCouplings_[i]->GetOutputEpsilon(k)) 
                    normsReached = FALSE;
                  
                  *oldVal = *val;
                } // end if
              counter++;              
            }
          }

        iter++;
        if (rPDE_.nonLinLogging_)
          Info->PrintF(rPDE_.pdename_, "\n"); 
      } // end while

#ifdef MpCCI
    if (actStep_==numTimeStep_)
      {
	rPDE_.PDEs_[0]->converged_ = TRUE;
	rPDE_.PDEs_[0]->CalcInputCoupling();
      }
#endif

  } 
  
  //----------------------- HARMONIC---------------------------------------
  void IterSolveStep::SolveStepHarmonic( const Boolean reset ) {

    ENTER_FCN( "IterSolveStep::SolveStepHarmonic" );

    Error( "Harmonic iterative coupling is not yet implemented", 
           __FILE__, __LINE__);
  }


  void IterSolveStep::SetActTime( const Double actTime ) {
    ENTER_FCN( "IterSolveStep::SetActTime() ");
    actTime_ = actTime;

    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++)
      rPDE_.PDEs_[i]->GetSolveStep()->SetActTime(actTime);
  }

  void IterSolveStep::SetActFreq( const Double actFreq ) {
    ENTER_FCN( "IterSolveStep::SetActFreq() ");

    actFreq_ = actFreq;
    
    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++)
      rPDE_.PDEs_[i]->GetSolveStep()->SetActFreq(actFreq);
  }

  void IterSolveStep::SetActStep( const UInt actStep ) {
    ENTER_FCN( "IterSolveStep::SetActStep() ");
    
    actStep_ = actStep;

    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
      rPDE_.PDEs_[i]->GetSolveStep()->SetActStep(actStep);
    }
  }


  void IterSolveStep::SetTimeStep( Double dt) {
    ENTER_FCN( "->GetSolveStep() ");

    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++)
      rPDE_.PDEs_[i]->GetSolveStep()->SetTimeStep(dt);
  }

  Double IterSolveStep::CalcNorm(NormType normtype, CFSVector & val, CFSVector & oldval)
  {
    ENTER_FCN( "IterSolveStep::CalcNorm" );

    // ATTENTION: Currently only working with Double-values
    // will be changed as soon as dynamic type information
    // is available

    Vector<Double> delta;
 
    Double norm, valNorm2;
  

    Vector<Double> & val_vec =\
      dynamic_cast<Vector<Double>& >(val);

    Vector<Double> & oldval_vec =\
      dynamic_cast<Vector<Double>& >(oldval);
  
    delta = val_vec - oldval_vec;

    switch (normtype)
      {
      case NO_NORM:
        return 0;
        break;
      
      case L2ABS:
        norm = delta.NormL2();
        break;

      case L2REL:
        valNorm2 =  val_vec.NormL2();
        if (valNorm2 > 0)
          norm = delta.NormL2() / valNorm2;
        else
          norm = delta.NormL2();

        break;
      }

    return norm;
  }
} // end of namespace
