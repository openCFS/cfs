#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

BlockSystem :: BlockSystem(Integer anumsys, Integer anumgraph)
  : AlgebraicSystem(anumsys, anumgraph)
{
#ifdef TRACE
  (*trace) << "entering BlockSystem::BlockSystem" << endl;
#endif

}
  
BlockSystem :: ~BlockSystem()
{
#ifdef TRACE
  (*trace) << "entering BlockSystem::~BlockSystem" << endl;
#endif
}

void BlockSystem :: CalcPrecond()
{
#ifdef TRACE
  (*trace) << "entering BlockSystem::CalcPrecond" << endl;
#endif

  Integer i;

  for (i=0; i<numsys; i++)
    {
      sysmat[i]->BuildInDirichlet(*rhs[i]);
      spemat[i]->BuildInDirichlet();
      premat[i]->Calc(*spemat[i], *auxmat[i]);
    }

#ifdef MEMTRACE
  (*memtrace) << "------------------------------------------------------------------" << endl;
  (*memtrace) << "+++ ALLOCATE MEMORY: double                   " << sumdmem << " MB" << endl;
  (*memtrace) << "+++ ALLOCATE MEMORY: integer                  " << sumimem << " MB" << endl;
  (*memtrace) << endl;
#endif
}


void BlockSystem :: CalcPrecond(Integer newprecond, Integer mat)
{
#ifdef TRACE
  (*trace) << "entering BlockSystem::CalcPrecond" << endl;
#endif

  switch (newprecond)
    {
    case 1: // nonlinear stuff
      premat[mat-1]->InitPrecond();
      spemat[mat-1]->Copy(sysmat[mat-1]);
      sysmat[mat-1]->BuildInDirichlet(*rhs[mat-1]);
      sysmat[mat-1]->BuildInDirichlet();
      spemat[mat-1]->BuildInDirichlet();

      premat[mat-1]->Calc(*spemat[mat-1], *auxmat[mat-1]);
     
      break;

    case 2:
      premat[mat-1]->InitPrecond();
      spemat[mat-1]->Copy(sysmat[mat-1]);
      sysmat[mat-1]->BuildInDirichlet(*rhs[mat-1]);
      sysmat[mat-1]->BuildInDirichlet();
      spemat[mat-1]->BuildInDirichlet();

      premat[mat-1]->Calc(*spemat[mat-1], *auxmat[mat-1]);
      break;
      
    case 3:
      sysmat[mat-1]->UpdateDirichletRHS(*rhs[mat-1]);
      break;

    default:
      sysmat[mat-1]->BuildInDirichlet(*rhs[mat-1]);
    }
#ifdef MEMTRACE
  (*memtrace) << "------------------------------------------------------------------" << endl;
  (*memtrace) << "+++ ALLOCATE MEMORY: double                   " << sumdmem << " MB" << endl;
  (*memtrace) << "+++ ALLOCATE MEMORY: integer                  " << sumimem << " MB" << endl;
  (*memtrace) << endl;
#endif
}

void BlockSystem :: Solve()
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::Solve" << endl;
#endif

  Integer i;

  for (i=0; i<numsys; i++)
    {
      scheme[i]->Calc(*sysmat[i], *premat[i], *rhs[i], *sol[i]);
    }
}

void BlockSystem :: Solve(Double * f, Double * u)
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::Solve" << endl;
#endif

  Integer i,dof,size;

  RealVector & v = (RealVector &) *sol[0];

  size = v.GetSize();
  dof  = v.GetNumDoF();
  
  rhs[0]->SetPointer(f);
  sol[0]->SetPointer(u);

  scheme[0]->Calc(*sysmat[0], *premat[0], *rhs[0], *sol[0]);

  for (i=0; i<size*dof; i++)
    {
      f[i] = v.Get(i+1);
    }
}

void BlockSystem :: CreateParameter()
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::CreateParameter" << endl;
#endif

  Integer i;

  for (i=0; i<numsys; i++)
    {
      param[i]  = new BaseParameter();
      param[i]->Set();
    }
}
}
