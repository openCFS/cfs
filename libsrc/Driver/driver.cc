#include <fstream>
#include <iostream>
#include <string>

#include "interface_gridlib.hh"
#include "acousticPDE.hh"
#include "driver.hh"

namespace CoupledField
{

template<class Dim>
Driver<Dim>::Driver(FileType * const aptFileType, Integer anummesh, Material * aptMaterial)
{
#ifdef TRACE
  (*trace) << "entering Driver::Driver" << std::endl;
#endif
  ptFileType=aptFileType;

  ptFileType->ReadNumStepsAndTimeSteps(numsteps, dt0);
//  ptFileType->ReadOutputOptions(SaveDer1, SaveDer2);
  SaveDer1=FALSE; SaveDer2=FALSE;

  ptgrid=new InterfaceGridlib<Dim>(ptFileType);
//  ptgrid=new GridInterfaceCFS<Dim>(ptFileType);
  ptgrid->Read();
  std::cout << " we do subdivosion " << std::endl;

  ptMaterial=aptMaterial;
}

template<class Dim>
void Driver<Dim>::SolveNewmarkMethod(OutResultUnverg<Dim> * ptUnverg)
{
#ifdef TRACE
  (*trace) << "entering Driver :: SolveNewmarkMethod" << std::endl;
#endif
 
/// Save the grid before a uniform refinement in a separate unverg-file 
//  OutResultUnverg * ptUnvergPreGrid=new OutResultUnverg("grid_pre"); 
//  ptUnvergPreGrid->Create(ptgrid,0);  
//  if (ptUnvergPreGrid) delete ptUnvergPreGrid;

  ptUnverg->Create(ptgrid,0);

//  Double endtime=1.0;   ////////////////////////////////////////////
  Double t=0;
  
  Double epsilon=1e-25; 
  ptAcPDE=new AcousticPDE(epsilon,dt0,ptgrid,0,ptMaterial, ptFileType);

  Integer i;
  for (i=0; i<numsteps; i++) 
    {
      t+=dt0;

      ptAcPDE->SolveNewmarkMethodStatic(t);
  
      ptUnverg->Dataset55(" fluid potential", ptAcPDE->getS(), i+1, t);
      if (SaveDer1) 
	ptUnverg->Dataset55(" fluid potential, 1st deriv., ", ptAcPDE->getS1(), i+1, t); 
      if (SaveDer2)
	ptUnverg->Dataset55(" fluid potential, 2nd deriv., ", ptAcPDE->getS2(), i+1, t);
    }
}

template<class Dim>
Driver<Dim> :: ~Driver()
{
#ifdef TRACE
  (*trace) << "entering Driver::~Driver" << std::endl;
#endif
  if (ptgrid) delete ptgrid;
  if (ptAcPDE) delete ptAcPDE;
}

}
