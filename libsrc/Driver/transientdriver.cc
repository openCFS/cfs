#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "transientdriver.hh"
#include "actimeerror.hh"
#include "vector.hh"

namespace CoupledField
{

TransientDriver :: TransientDriver(Domain * adomain)
:BaseDriver(adomain)
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::TransientDriver" << std::endl;
#endif

  // get time steps information from conf-file

  conf->get("numsteps",numstep_);
  conf->get("firstdt", firstdt_);
  conf->get("stepsavebeg",isavebegin_);
  conf->get("stepsaveend",isaveend_);
  conf->get("stepsaveincr",isaveincr_);
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

   Vector<Double> derAdt=ptdomain_->GetPDE(pdenumber)->getS2()-ptdomain_->GetPDE(pdenumber)->getS2old();

   Double size=ptdomain_->GetPDE(pdenumber)->getSize();
   Double error1=sqrt(derAdt*derAdt/size)*dt*dt*0.0833;

// l2 norm solution
   Vector<Double> sol=ptdomain_->GetPDE(pdenumber)->getS();

   std::cout << " Error " << error1 << " " << sol.norm_2() << " ratio " << error1/sol.norm_2() << std::endl; 

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

//  TimeErrorEstimator ** ptTimeError=new TimeErrorEstimation * [pdenumber];
  TimeErrorEstimator * ptTimeError;
 
  ptTimeError=ptdomain_->GetPDE(pdenumber)->CreatePtTimeError();

  Double dt=firstdt_;
  Boolean resetsysmat=FALSE;

  ptdomain_->GetPDE(pdenumber)->CalcParameters(dt);
  ptdomain_->GetPDE(pdenumber)->SetMatrixFactors();

/*  if (InfoPrint)
    (*infofile) << "# step      dt      timestep " << std::endl << 0 << "  " << dt << std::endl;
*/

 Integer nstep=0;

 Integer startrepeat;
 conf->get("startrepeat",startrepeat,"Acoustic");  

 Integer i;
 Double steptime_prev=steptime;

 for (i=0; i<startrepeat; i++)
 {
 ptdomain_->GetPDE(pdenumber)->SolveStepTrans(ptdomain_->GetBCs(), nstep, steptime, level, resetsysmat);

 // writing results in output-file
 ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
 if (InfoPrint)
  (*infofile) << nstep << " " << dt << " " << steptime << std::endl;

 ptTimeError->CalcThirdDer();
 ptTimeError->CalcError(dt);

 nstep++;
 steptime_prev=steptime;
 steptime+=dt;
 }

  ptdomain_->GetPDE(pdenumber)->SolveStepTrans(ptdomain_->GetBCs(), nstep, steptime, level, resetsysmat);

do
  {
  // test error
   if (ptTimeError->TestError(dt))
      {
        Double prev_dt=dt;
        ptTimeError->ChangeStep(dt);

//        if (prev_dt < dt)
//       {
       ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
   if (InfoPrint)
    (*infofile) << nstep << " " << dt << " " << steptime << " coarse" << std::endl;
       nstep++;
       steptime+=dt;
       std::cout << " COARSE " << std::endl;
//       }
//       else
//       {
//        steptime=steptime-prev_dt;
//       ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
//   if (InfoPrint)
//    (*infofile) << nstep << " " << dt << " " << steptime << " refine" << std::endl;
//       steptime+=dt;
//       nstep++;
//       std::cout << " REFINE " << std::endl;
//       }

       
        std::cout << " dt " << dt << " step " << nstep << std::endl;

        ptdomain_->GetPDE(pdenumber)->CalcParameters(dt);
        ptdomain_->GetPDE(pdenumber)->SetMatrixFactors();

        resetsysmat=TRUE;

       // print info in file.info
/*         if (InfoPrint)
         { 
        if (prev_dt < dt)
          (*infofile) <<  nstep+1 << "       " << dt << "     " << steptime << "     change    coarse" << std::endl;
            else
           (*infofile) <<  nstep+1<< "       " << dt << "     " << steptime << "     change    refine" << std::endl;
         }
*/
      }
   else
     {
      ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
   if (InfoPrint)
    (*infofile) << nstep << " " << dt << " " << steptime << std::endl;
      nstep++;
      steptime+=dt;

/*       // print info in file.info
         if (InfoPrint)
          (*infofile) << nstep << "         " << dt << "      " << steptime << std::endl;
*/
   
  }

  ptdomain_->GetPDE(pdenumber)->SolveStepTrans(ptdomain_->GetBCs(), nstep, steptime, level, resetsysmat);

  std::cout << steptime << " steptime " << dt << " dt " << std::endl;
  } 
while ( steptime <= endtime);

  ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
   if (InfoPrint)
    (*infofile) << nstep << " " << steptime << " " << dt << std::endl;

  std::cout << " number of steps " << nstep << std::endl;
}

} // end of namespace
