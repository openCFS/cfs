#include <fstream>
#include <iostream>

#include <general_head.hh>
#include <utils_head.hh>
#include <datainout_head.hh>
#include <elements_head.hh>
#include <forms_head.hh>
#include <linalg_head.hh>
#include <domain_head.hh>
#include <pde_head.hh>

#include "driver.hh"

namespace CoupledField
{

Driver::Driver(FileType * const aptFileType, Integer anummesh, Material * aptMaterial)
{
#ifdef TRACE
  (*trace) << "entering Driver::Driver" << std::endl;
#endif
  ptFileType=aptFileType;

  ptFileType->ReadNumStepsAndTimeSteps(numsteps, dt0);
//  ptFileType->ReadOutputOptions(SaveDer1, SaveDer2);
  SaveDer1=FALSE; SaveDer2=FALSE;

  ptgrid=new GridInterfaceCFS<Point2D>(ptFileType);
  ptMaterial=aptMaterial;
}

void Driver::SolveNewmarkMethod(OutResultUnverg * ptUnverg)
{
#ifdef TRACE
  (*trace) << "entering Driver :: SolveNewmarkMethod" << std::endl;
#endif
  ptUnverg->Create(ptgrid,0);

//  Double endtime=1.0;   ////////////////////////////////////////////
  Double t=0;
  
  Double epsilon=1e-25; 
  ptAcPDE=new AcousticPDE(epsilon,dt0,ptgrid,ptMaterial, ptFileType);

  Integer i;
  for (i=0; i<numsteps; i++) 
    {
      t+=dt0;

      std::cout << "begin of step number " << i << std::endl;

      ptAcPDE->SolveNewmarkMethodStatic(t);
  
      std::cout << i << " " << ptAcPDE->getS().size() << std::endl;
      ptUnverg->Dataset55(" fluid potential", ptAcPDE->getS(), i+1, t);
      if (SaveDer1) 
	ptUnverg->Dataset55(" fluid potential, 1st deriv., ", ptAcPDE->getS1(), i+1, t); 
      if (SaveDer2)
	ptUnverg->Dataset55(" fluid potential, 2nd deriv., ", ptAcPDE->getS2(), i+1, t);
    }
}

Driver :: ~Driver()
{
#ifdef TRACE
  (*trace) << "entering Driver::~Driver" << std::endl;
#endif
  if (ptgrid) delete ptgrid;
  if (ptAcPDE) delete ptAcPDE;
}

}
