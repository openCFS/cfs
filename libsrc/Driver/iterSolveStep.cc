#include "iterSolveStep.hh"

#include "PDE/StdPDE.hh"

namespace CoupledField
{
  //! Derived class for step-wise solving of iterative coupled StdPDEs

  IterSolveStep::IterSolveStep(IterCoupledPDE &apde) 
    : BaseSolveStep(),
      rPDE_(apde),
      rCouplings_(apde.Couplings_),
      actlevel_(apde.actlevel_){

    ENTER_FCN( "IterSolveStep::IterSolveStep" );
    

  }

  IterSolveStep::~IterSolveStep() {
    
    ENTER_FCN( "IterSolveStep::~IterSolveStep" );

  }


  //----------------------- STATIC---------------------------------------

  void IterSolveStep::SolveStepStatic(const Integer kstep, const Double aTime,
				       const Integer level,
				       const Boolean updatesysmat ) {

    ENTER_FCN ( "IterSolveStep::SolveStepStatic" );
  
    CFSVector *val, *oldVal;
    Integer iter = 0;
    Integer counter = 0;
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
      
	for (Integer i=0; i<rPDE_.PDEs_.GetSize(); i++)
	  {
	    if (rPDE_.nonLinLogging_)
	    Info->PrintF(rPDE_.pdename_, " Processing PDE %s\n", 
			 (rPDE_.PDEs_[i]->GetName()).c_str());

	    // Only solve current PDE, if the corresponding
	    // flag in 'solvePDE_' is set to TRUE
	    if (rPDE_.solvePDE_[i] == TRUE) {
	      rPDE_.PDEs_[i]->GetSolveStep()->PreStepStatic(kstep,aTime,actlevel_,updatesysmat);
	      rPDE_.PDEs_[i]->CalcInputCoupling();
	      rPDE_.PDEs_[i]->GetSolveStep()->SolveStepStatic(kstep,aTime,actlevel_,updatesysmat);
	      rPDE_.PDEs_[i]->GetSolveStep()->PostStepStatic(kstep,aTime,actlevel_);
	      rPDE_.PDEs_[i]->CalcOutputCoupling();
	      
	      // Calculate Norms
	      for (Integer k=0; k<rCouplings_[i]->GetNumOutputCouplings(); k++)
		{
		  rCouplings_[i]->GetOutputValues(k, val);
		  rCouplings_[i]->GetOutputOldValues(k, oldVal);
		  rPDE_.norms_[counter] = CalcNorm(rCouplings_[i]->GetOutputNormType(k), *val, *oldVal);

		  if (rPDE_.nonLinLogging_)
		    {
		      Enum2String(rCouplings_[i]->GetOutputQuantity(k), quantityConv);
		      Info->PrintF(rPDE_.pdename_, " %s : Norm of %s = %g\n", 
				   (rCouplings_[i]->GetPDEName()).c_str(),
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
    for (Integer i=0; i<rPDE_.PDEs_.GetSize(); i++)
      rPDE_.PDEs_[i]->PostProcess(actlevel_);

  }


  //----------------------- TRANSIENT---------------------------------------
  void IterSolveStep::SolveStepTrans(const Integer kstep, const Double asteptime, const Integer level, 
				      const Boolean updatesysmat)
  {
    ENTER_FCN( "IterSolveStep::SolveStepTrans" );

    Double  steptime  = asteptime;

    Integer iter = 0;
    Boolean normsReached = FALSE;
    std::string quantityConv;

    // In the beginning of each time step
    // the coupling data has to be reseted
    for (Integer i=0; i<rPDE_.PDEs_.GetSize(); i++)
      rPDE_.PDEs_[i]->ResetCoupling();

    while (iter < rPDE_.maxiter_ &&  (! normsReached))
      {
	if (rPDE_.nonLinLogging_)
	  {
	    Info->PrintF(rPDE_.pdename_,"\n"); 
	    Info->PrintF(rPDE_.pdename_, " COUPLED ITERATION %i =================================\n", 
			 iter+1);
	  }

	Integer counter = 0;
	normsReached = TRUE;
      
	for (Integer i=0; i<rPDE_.PDEs_.GetSize(); i++)
	  {
	    if (rPDE_.nonLinLogging_)
	      Info->PrintF(rPDE_.pdename_, " Processing PDE %s\n", 
			   (rPDE_.PDEs_[i]->GetName()).c_str());

	    // Only solve current PDE, if the corresponding
	    // flag in 'solvePDE_' is set to TRUE
	    if (rPDE_.solvePDE_[i] == TRUE) {
	      
	      rPDE_.PDEs_[i]->GetSolveStep()->PreStepTrans(kstep, steptime, level, updatesysmat);
	      rPDE_.PDEs_[i]->CalcInputCoupling();
	      rPDE_.PDEs_[i]->GetSolveStep()->SolveStepTrans(kstep, steptime, level, updatesysmat);
	      rPDE_.PDEs_[i]->GetSolveStep()->PostStepTrans(kstep, steptime, level);
	      rPDE_.PDEs_[i]->CalcOutputCoupling();
	      
	      // Calculate Norms
	      for (Integer k=0; k<rCouplings_[i]->GetNumOutputCouplings(); k++)
		{
		  CFSVector *val, *oldVal;
		  rCouplings_[i]->GetOutputValues(k, val);
		  rCouplings_[i]->GetOutputOldValues(k, oldVal);
		  rPDE_.norms_[counter] = CalcNorm(rCouplings_[i]->GetOutputNormType(k), *val, *oldVal);

		  if (rPDE_.nonLinLogging_)
		    {
		      Enum2String(rCouplings_[i]->GetOutputQuantity(k), quantityConv);
		      Info->PrintF(rPDE_.pdename_, " %s : Norm of %s = %g\n", 
				   (rCouplings_[i]->GetPDEName()).c_str(),
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
      }
  } 
  
  //----------------------- HARMONIC---------------------------------------
  void IterSolveStep::SolveStepHarmonic(const Integer freqStep, const Double frequency, 
					Integer level, const Boolean reset) {

    ENTER_FCN( "IterSolveStep::SolveStepHarmonic" );

    Error( "Harmonic iterative coupling is not yet implemented", 
	   __FILE__, __LINE__);
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
