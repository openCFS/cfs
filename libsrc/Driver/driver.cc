#include <fstream>
#include <iostream>
#include <string>

//#include "interface_gridlib.hh"
#include "acousticPDE.hh"
#include "driver.hh"

#include "outGMV.hh"

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
//  numsteps=30; dt0=2.000000000000E-07;

  ptFileType->ReadOutputOptions(SaveDer1, SaveDer2);
  SaveDer1=FALSE; SaveDer2=FALSE;

//  ptgrid=new InterfaceGridlib<Dim>(ptFileType);
  ptgrid=new GridInterfaceCFS<Dim>(ptFileType);
  ptgrid->Read();

  ptMaterial=aptMaterial;
}

template<class Dim>
void Driver<Dim>::SolveNewmarkMethod(OutResultUnverg<Dim> * ptUnverg)
{
#ifdef TRACE
  (*trace) << "entering Driver :: SolveNewmarkMethod" << std::endl;
#endif
 
/// Save the grid before a uniform refinement in a separate unverg-file 
/*
  OutResultUnverg<Dim> * ptUnvergPreGrid=new OutResultUnverg<Dim>("grid_pre"); 
  ptUnvergPreGrid->Create(ptgrid,0);  
  if (ptUnvergPreGrid) delete ptUnvergPreGrid;
  ptgrid->SubdivideUniform(0);
*/
   ptUnverg->Create(ptgrid,0);

   OutGMV<Dim> * ptGMV=new OutGMV<Dim>("test",ptgrid);
   ptGMV->Write(0);
   delete ptGMV;

//  Double endtime=1.0;   ////////////////////////////////////////////
  Double t=0;
  
  Double epsilon=1e-25; 
  ptAcPDE=new AcousticPDE<Dim>(epsilon,dt0,ptgrid,0,ptMaterial, ptFileType);

  Integer i;
  for (i=0; i<numsteps; i++) 
    {
      t+=dt0;

      ptAcPDE->SolveNewmarkMethodStatic(t);

      PrintResultsUnverg(ptUnverg, ptAcPDE,i,t);  
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

template<class Dim>
void Driver<Dim> :: PrintResultsUnverg(OutResultUnverg<Dim> * ptUnverg, PDE * ptPDE, const Integer step, const Double t)
{

      ptUnverg->Dataset55(" fluid potential", ptPDE->getS(), step+1, t);

      if (SaveDer1)
      ptUnverg->Dataset55(" fluid potential, 1st deriv., ", ptPDE->getS1(),
step+1, t);

      if (SaveDer2)
      ptUnverg->Dataset55(" fluid potential, 2nd deriv., ", ptPDE->getS2(),
step+1, t);

}

} // end of namespace
