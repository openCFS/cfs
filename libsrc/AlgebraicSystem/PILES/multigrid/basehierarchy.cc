#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include "environment.hh"
#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

BaseHierarchy :: BaseHierarchy(BaseParameter & aparam, Integer asize, Integer aoffset, 
			       Integer anumrhs, Integer adof)
{
#ifdef TRACE
  (*trace) << "entering BaseHierarchy::BaseHierarchy" << endl;
#endif

  param    = &aparam;
  offset   = aoffset;
  numrhs   = anumrhs;
  dof      = adof;
  amglevel = 0;
  gmglevel = 0;

  Integer i;

  for (i=0; i<20; i++)
    {
      alpha[i] = param->GetAlpha();
    }

  coarsesys = param->GetCoarseSystem();
  epsmat    = param->GetEpsMat();
  rsw       = new Integer[asize];

  for (i=0; i<20; i++)
    {
      node[i].sysmat = NULL;
      node[i].auxmat = NULL;
      node[i].vec1   = NULL;
      node[i].vec2   = NULL;
      node[i].vec3   = NULL;
      node[i].vec4   = NULL;
      node[i].syspro = NULL;
      node[i].sysres = NULL;
      node[i].auxpro = NULL;
      node[i].auxres = NULL;
      node[i].coarse = NULL;
      node[i].smooth = NULL;
    }
}
  
BaseHierarchy :: ~BaseHierarchy()
{
#ifdef TRACE
  (*trace) << "entering BaseHierarchy::~BaseHierarchy" << endl;
#endif

  Integer i;

  if (node[offset].vec1 != NULL)
    {
      delete node[offset].vec1;
      node[offset].vec1 = NULL;
    }
  
  if (node[offset].vec2 != NULL)
    {
      delete node[offset].vec2;
      node[offset].vec2 = NULL;
    }
  
  if (node[offset].vec3 != NULL)
    {
      delete node[offset].vec3;
      node[offset].vec3 = NULL;
    }
  
  if (node[offset].vec4 != NULL)
    {
      delete node[offset].vec4;
      node[offset].vec4 = NULL;
    }
  
  if (node[offset].syspro != NULL)
    {
      delete node[offset].syspro;
      node[offset].syspro = NULL;
      node[offset].sysres = NULL;
    }
  
  if (node[offset].auxpro != NULL)
    {
      delete node[offset].auxpro;
      node[offset].auxpro = NULL;
      node[offset].auxres = NULL;
    }
  
  if (node[offset].coarse != NULL)
    {
      delete node[offset].coarse;
      node[offset].coarse = NULL;
    }
  
  if (node[offset].smooth != NULL)
    {
      delete node[offset].smooth;
      node[offset].smooth = NULL;
    }

  for (i=1; i<amglevel; i++)
    {
      level = offset-i;

      if (node[level].sysmat != NULL)
	{
	  delete node[level].sysmat;
	  node[level].sysmat = NULL;
	}

      if (node[level].auxmat != NULL)
	{
	  delete node[level].auxmat;
	  node[level].auxmat = NULL;
	}

      if (node[level].vec1 != NULL)
	{
	  delete node[level].vec1;
	  node[level].vec1 = NULL;
	}
      
      if (node[level].vec2 != NULL)
	{
	  delete node[level].vec2;
	  node[level].vec2 = NULL;
	}

      if (node[level].vec3 != NULL)
	{
	  delete node[level].vec3;
	  node[level].vec3 = NULL;
	}

      if (node[level].vec4 != NULL)
	{
	  delete node[level].vec4;
	  node[level].vec4 = NULL;
	}

      if (node[level].syspro != NULL)
	{
	  delete node[level].syspro;
	  node[level].syspro = NULL;
	  node[level].sysres = NULL;
	}

      if (node[level].auxpro != NULL)
	{
	  delete node[level].auxpro;
	  node[level].auxpro = NULL;
	  node[level].auxres = NULL;
	}

      if (node[level].coarse != NULL)
	{
	  delete node[level].coarse;
	  node[level].coarse = NULL;
	}
      
      if (node[level].smooth != NULL)
	{
	  delete node[level].smooth;
	  node[level].smooth = NULL;
	}
    }
  
  delete [] rsw;
}

void BaseHierarchy :: Calc(BaseMatrix & sysmat, BaseMatrix & auxmat)
{
#ifdef TRACE
  (*trace) << "entering BaseHierarchy::Calc" << endl;
#endif
  
  level     = offset-amglevel;    
  fsize_sys = sysmat.GetSize();
  fnne_sys  = sysmat.GetNNE();

  (*cla) << amglevel+1 << " " << fsize_sys << " " << fnne_sys << endl;

  if (fsize_sys > coarsesys)
    {
      sysmat.SetAuxMatrix(auxmat);

      node[level].sysmat = &sysmat;
      node[level].auxmat = &auxmat;
      
      fsize_aux = auxmat.GetSize();
      fnne_aux  = auxmat.GetNNE();

      /// define the smoother

      MakeSmoother();

      node[level].smooth->Calc(*node[level].sysmat);

      /// define the structure for coarsening

      MakeCoarsening();

      /// setup the neighbours and the other sets and perform a splitting

      /// (*cla) << "after smoother and coarse" << endl;

      node[level].coarse->SetNeighbour(alpha[level], epsmat);
      node[level].coarse->Calc(rsw);
      node[level].coarse->SetSize();

      /// get the sizes

      csize_sys = node[level].coarse->GetSysCSize();
      csize_aux = node[level].coarse->GetAuxCSize();
      cnne_sys  = node[level].coarse->GetSysCNNE();
      cnne_aux  = node[level].coarse->GetAuxCNNE();

      /// construct a temporary prolongation - auxiliary

      MakeAuxTransfer();

      /// construct a temporary matrix - auxiliary
      
      MakeAuxMatrix();

      ///(*cla) << "after aux matrix" << endl;

      /// calculate the prolongation/restriction operator

      node[level].auxpro->Calc(node[level].coarse->GetTopology());
      node[level].auxres->Calc(node[level].coarse->GetTopology());

      /// calculate the auxiliary coarse grid matrix

      node[level-1].auxmat->SetCoarseGraph(node[level].coarse->GetTopology());
      node[level].auxmat->Galerkin(node[level].auxres, node[level].auxpro, *node[level-1].auxmat);

      //node[level-1].auxmat->Print();

      /// topology for edge discretization, number of weights, etc.

      node[level].coarse->CalcTopology();

      /// construct the system prolongation operator

      MakeSysTransfer();

      /// construct the coarse system matrix

      MakeSysMatrix();

      /// calculate the prolongation/restriction operator

      node[level].syspro->Calc(node[level].coarse->GetTopology());
      node[level].sysres->Calc(node[level].coarse->GetTopology());

      /// perform the Galerkin projection

      node[level-1].sysmat->SetCoarseGraph(node[level].coarse->GetTopology());
      node[level].sysmat->Galerkin(node[level].sysres, node[level].syspro, *node[level-1].sysmat);

      //node[level-1].sysmat->Print();

      /// allocate the vectors for MG

      MakeVector();

      /// delete the memory from coarsening

      /// (*cla) << "after make vector" << endl;
      
      DeleteAuxMemory();

      /// recursive call for Calc()

      amglevel++;

      /// (*cla) << "before calc" << endl;

      Calc(*node[level-1].sysmat, *node[level-1].auxmat);
    }
  else
    {
      ///(*cla) << "coarse grid solver" << endl;

      node[level].sysmat = &sysmat;

      cputime = cpuclock.GetTime();
      node[level].sysmat->Factor();
      cputime = cpuclock.GetTime()-cputime;

      MakeVector();

      (*cla) << "### coarse grid solver CPU-time " << cputime << endl;

      amglevel++;
      
      return;
    }
}

}
