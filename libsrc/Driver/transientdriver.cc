#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "transientdriver.hh"
#include "actimeerror.hh"
#include "vector.hh"
#include "spaceerror.hh"

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

/*
void TransientDriver :: SolveProblemAdapt()
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::SolveProblemAdapt" << std::endl;
#endif

  Integer i,level=0;
  Integer pdenumber  = 0;
  Integer nsys = 0;
  Double steptime=firstdt_;
  Double steptime_prev=0;
  Integer nstep=0;

  // calculation of end-time
  Double endtime=numstep_*firstdt_;

//  TimeErrorEstimator ** ptTimeError=new TimeErrorEstimation * [pdenumber];
  TimeErrorEstimator * ptTimeError;
 
  ptTimeError=ptdomain_->GetPDE(pdenumber)->CreatePtTimeError();

  Double dt=firstdt_;
  Boolean resetsysmat=FALSE;

  ptdomain_->GetPDE(pdenumber)->CalcParameters(dt);
  ptdomain_->GetPDE(pdenumber)->SetMatrixFactors();

 Integer startrepeat;
 conf->get("startrepeat",startrepeat,"Acoustic");  

 for (i=0; i<startrepeat; i++)
 {
 ptdomain_->GetPDE(pdenumber)->SolveStepTrans(ptdomain_->GetBCs(), nstep, steptime, level, resetsysmat);

 ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();

 ptTimeError->CalcThirdDer();
 ptTimeError->CalcError(dt);

 nstep++;
 steptime_prev=steptime;
 steptime+=dt;
 }

  ptdomain_->GetPDE(pdenumber)->SolveStepTrans(ptdomain_->GetBCs(), nstep, steptime, level, resetsysmat);

  ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();

do
  {
   if (ptTimeError->TestError(dt))
      {
        Double prev_dt=dt;
        ptTimeError->ChangeStep(dt);

       nstep++;
       steptime_prev=steptime;
       steptime+=dt;

       ptdomain_->GetPDE(pdenumber)->CalcParameters(dt);
       ptdomain_->GetPDE(pdenumber)->SetMatrixFactors();

       resetsysmat=TRUE;
      }
   else
     {
      nstep++;
      steptime_prev=steptime;
      steptime+=dt;
     resetsysmat=FALSE;
     }

  ptdomain_->GetPDE(pdenumber)->SolveStepTrans(ptdomain_->GetBCs(), nstep, steptime, level, resetsysmat);
  ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();

  std::cout << steptime << " steptime " << dt << " dt " << std::endl;
  } 
while ( steptime <= endtime);

  std::cout << " number of steps " << nstep << std::endl;
}
*/

void TransientDriver :: SolveProblemAdapt()
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::SolveProblemAdapt" << std::endl;
#endif

  Integer i,level=0;
  Integer pdenumber  = 0;
  Integer nsys = 0;
  Double steptime=firstdt_;
  Double steptime_prev=0;
  Integer nstep=0;

  // calculation of end-time
  Double endtime=numstep_*firstdt_;

//  TimeErrorEstimator ** ptTimeError=new TimeErrorEstimation * [pdenumber];
  TimeErrorEstimator * ptTimeError;

  ptTimeError=ptdomain_->GetPDE(pdenumber)->CreatePtTimeError();

  Double dt=firstdt_;
  Boolean resetsysmat=FALSE;

  ptdomain_->GetPDE(pdenumber)->CalcParameters(dt);
  ptdomain_->GetPDE(pdenumber)->SetMatrixFactors();

 Integer startrepeat;
 conf->get("startrepeat",startrepeat,"Acoustic");

 for (i=0; i<startrepeat; i++)
 {
 ptdomain_->GetPDE(pdenumber)->SolveStepTrans(ptdomain_->GetBCs(), nstep, steptime, level, resetsysmat);

 ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
 if (InfoPrint)
  (*infofile) << steptime << " " << dt << std::endl;

 // ptTimeError->CalcThirdDer();
 ptTimeError->CalcError(dt);

 nstep++;
 steptime_prev=steptime;
 steptime+=dt;
 }

  ptdomain_->GetPDE(pdenumber)->SolveStepTrans(ptdomain_->GetBCs(), nstep, steptime, level, resetsysmat);

do
  {
   if (ptTimeError->TestError(dt))
      {
	std::cout << " Test error is truw " << std::endl;
        Double prev_dt=dt;
        ptTimeError->ChangeStep(dt);

       if (prev_dt < dt)
       {
       ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();

 if (InfoPrint)
  (*infofile) << steptime << " " << dt << std::endl;

       nstep++;
       steptime_prev=steptime;
       steptime+=dt;
       std::cout << " COARSE " << std::endl;
      }
      else
      {
	std::cout << " we do refinement " << steptime << " dt " << dt << " " <<std::endl;
       steptime=steptime_prev;
       steptime+=dt;
       std::cout << " REFINE " << std::endl;
      }
        ptdomain_->GetPDE(pdenumber)->CalcParameters(dt);
        ptdomain_->GetPDE(pdenumber)->SetMatrixFactors();

        resetsysmat=TRUE;
      }
   else
     {
      ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
   if (InfoPrint)
    (*infofile) << nstep << " " << dt << " " << steptime << std::endl;
      nstep++;
      steptime_prev=steptime;
      steptime+=dt;
     resetsysmat=FALSE;

  }

  ptdomain_->GetPDE(pdenumber)->SolveStepTrans(ptdomain_->GetBCs(), nstep, steptime, level, resetsysmat);

  std::cout << steptime << " steptime " << dt << " dt " << std::endl;
  }
while ( steptime <= endtime);

  ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
 if (InfoPrint)
  (*infofile) << steptime << " " << dt << std::endl;

  std::cout << " number of steps " << nstep << std::endl;
}

void TransientDriver :: SolveProblemAdaptSpace()
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::SolveProblemAdaptSpace" << std::endl;
#endif

  Integer level=0;
  Integer pdenumber = 0;
  Integer nsys = 0;
  Double steptime=firstdt_;
  Integer stepsave=isavebegin_-1;

  Integer maxnumrepeat, numrepeat=0;
  conf->get("maxnumrepeat",maxnumrepeat,"SpaceAdaptivity");

  Double dt=firstdt_;
  Boolean updatesysmat=FALSE;

  SpaceErrorEstimator * ptSpaceError;
  ptSpaceError=ptdomain_->GetPDE(pdenumber)->CreatePtSpaceError();

  ptdomain_->GetPDE(pdenumber)->CalcParameters(dt);
  ptdomain_->GetPDE(pdenumber)->SetMatrixFactors();

  Integer nstep;
  for (nstep = 0; nstep<numstep_; nstep++)
    {
      ptdomain_->GetPDE(pdenumber)->SolveStepTrans(ptdomain_->GetBCs(), nstep, steptime, level, updatesysmat);

      numrepeat=0;
      while (ptSpaceError->TestError() && numrepeat != maxnumrepeat)
      {
       ptSpaceError->RefineMesh();
       ptdomain_->Update();

  /*
  Char * name="testref";
  WriteResults * ptOut=new WriteResultsUnverg(name);
  ptOut->Init(ptdomain_->GetGrid());
  ptOut->WriteGrid(1);
  if (ptOut) delete ptOut;
  exit(1);
  */

      ptdomain_->GetPDE(pdenumber)->RestoreSol();
  
      ptdomain_->GetPDE(pdenumber)->SolveStepTransNewMesh(ptdomain_->GetBCs(), nstep, steptime, level);

      numrepeat++;     
       }
     
    if (nstep == stepsave && (nstep < isaveend_))
      {
        ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
        stepsave+=isaveincr_;
      }
   steptime+=dt;
   }
}
 
} // end of namespace
