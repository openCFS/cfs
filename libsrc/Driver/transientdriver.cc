#include <fstream>
#include <iostream>
#include <string>

#include "transientdriver.hh"
#include "actimeerror.hh"

namespace CoupledField
{

TransientDriver :: TransientDriver(Domain<Point2D> * adomain)
:BaseDriver<Point2D>(adomain)
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::TransientDriver" << std::endl;
#endif

  // read step data from input file
  ptdomain_->GetInFile()->ReadStepData(numstep_,isavebegin_,isaveend_,isaveincr_,firstdt_);

}

TransientDriver :: ~TransientDriver()
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::~TransientDriver" << std::endl;
#endif

}


void TransientDriver :: SetupMatricesPDE(const Integer pdenumber,const Integer type)
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::SetUpMatricesPDE" << std::endl;
#endif
 
  ptdomain_->GetPDE(pdenumber)->SetupMatrices(type);

}

void TransientDriver :: SolveProblem()
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::SolveProblem" << std::endl;
#endif

  Integer level=0;
  Integer pdenumber  = 0;
  Integer nsys = 0;
  Double steptime=firstdt_;
  Integer stepsave=isavebegin_-1;

  Double dt=firstdt_;
  Boolean updatesysmat=FALSE;

  ptdomain_->GetPDE(pdenumber)->CalcParameters(dt);
  ptdomain_->GetPDE(pdenumber)->SetMatrixFactors();

  Integer nstep;
  for (nstep = 0; nstep<numstep_; nstep++)
    {
      ptdomain_->GetPDE(pdenumber)->SolveStepTrans(ptdomain_->GetBCs(), nstep, steptime, level, updatesysmat);

   // writing results in output-file
    if (nstep == stepsave && (nstep < isaveend_))
      { 
        ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
        stepsave+=isaveincr_;
      }

   steptime+=dt;
   }
}

void TransientDriver :: SolveProblemAdapt()
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::SolveProblemAdapt" << std::endl;
#endif

  Integer level=0;
  Integer pdenumber  = 0;
  Integer nsys = 0;
  Double steptime=firstdt_;

  // calculation of end-time
  Double endtime=numstep_*firstdt_;

  std::cout << "EndTime" << endtime << std::endl;

//  TimeErrorEstimator ** ptTimeError=new TimeErrorEstimation * [pdenumber];
  TimeErrorEstimator * ptTimeError;
 
  ptTimeError=ptdomain_->GetPDE(pdenumber)->CreatePtTimeError();

  Double dt=firstdt_;
  Boolean resetsysmat=FALSE;

  ptdomain_->GetPDE(pdenumber)->CalcParameters(dt);
  ptdomain_->GetPDE(pdenumber)->SetMatrixFactors();

  if (InfoPrint)
   (*infofile) << "# step   " << " timestep " << std::endl << 0 << "  " << dt << std::endl;
  
Integer nstep=0;
do
  {
      ptdomain_->GetPDE(pdenumber)->SolveStepTrans(ptdomain_->GetBCs(), nstep, steptime, level, resetsysmat);

   // writing results in output-file
      ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();

      nstep++;
      steptime+=dt; 

  // test error
   if (ptTimeError->TestError(dt))
      {
         Double prev_dt=dt;
         ptTimeError->ChangeStep(dt);

         ptdomain_->GetPDE(pdenumber)->CalcParameters(dt);
         ptdomain_->GetPDE(pdenumber)->SetMatrixFactors();

         resetsysmat=TRUE;

       // print info in file.info
         if (InfoPrint)
         { 
            (*infofile) << " steps" << prev_dt << " " << dt << std::endl;

        if (prev_dt < dt)
          (*infofile) <<  nstep+1 << "  " << steptime << "     change    coarse" << std::endl;
            else
           (*infofile) <<  nstep+1 << "  " << steptime << "     change    refine" << std::endl;
         }

      }
   else
     {
       // print info in file.info
         if (InfoPrint)
          (*infofile) << nstep << "  " << steptime << std::endl;
     }
   } while (steptime <= endtime);

  std::cout << " number of steps " << nstep << std::endl;
}

} // end of namespace
