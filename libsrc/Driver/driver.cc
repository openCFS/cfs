#include <fstream>
#include <iostream>
#include <string>

//#include "interface_gridlib.hh"
#include "interface_gridcfs.hh"
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
  ptgrid=NULL;

  ptFileType->ReadNumStepsAndTimeSteps(numsteps, dt0);
//  numsteps=30; dt0=2.000000000000E-07;

  ptFileType->ReadOutputOptions(SaveDer1, SaveDer2);
  SaveDer1=FALSE; SaveDer2=FALSE;

  // read from conf-file the type of using grid-library
  std::string libmesh;
  conf->get("mesh_library",libmesh);

  std::cout << libmesh << std::endl;

  // initialize a pointer to an object with an infromation about mesh  
//  if (libmesh =="gridlib") ptgrid=new InterfaceGridlib<Dim>(ptFileType);
//  else
  if (libmesh =="cfsgrid") ptgrid=new GridInterfaceCFS<Dim>(ptFileType);
   else
     Error("Unknown type of mesh_library in conf-file",__FILE__,__LINE__);

  // read initial grid
  ptgrid->Read();

  ptMaterial=aptMaterial;
}

template<class Dim>
void Driver<Dim>::SolveNewmarkMethod(WriteResults * ptOutput)
{
#ifdef TRACE
  (*trace) << "entering Driver :: SolveNewmarkMethod" << std::endl;
#endif
 
/// Save the grid before a uniform refinement in a separate unverg-file 
/*
  OutResultUnverg<Dim> * ptOutputPreGrid=new OutResultUnverg<Dim>("grid_pre"); 
  ptOutputPreGrid->Create(ptgrid,0);  
  if (ptOutputPreGrid) delete ptOutputPreGrid;
  ptgrid->SubdivideUniform(0);
*/
   ptOutput->Init(ptgrid);
   ptOutput->WriteGrid(0);
/*
   WriteResultsGMV<Dim> * ptGMV=new WriteResultsGMV<Dim>("test",ptgrid);
   ptGMV->WriteGrid(0);
   delete ptGMV;
*/
//  Double endtime=1.0;   ////////////////////////////////////////////
  Double t=0;
/*  
  Double epsilon=1e-25; 
  ptAcPDE=new AcousticPDE<Dim>(dt0,ptgrid,0,ptMaterial, ptFileType);

  Integer i;
  for (i=0; i<numsteps; i++) 
    {
      t+=dt0;

      ptAcPDE->SolveNewmarkMethodStatic(t);

      ptAcPDE->WriteResultsInFile();
    }
*/
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

} // end of namespace
