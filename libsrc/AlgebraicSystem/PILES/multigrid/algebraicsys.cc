#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

AlgebraicSystem :: AlgebraicSystem(Integer anumsys, Integer anumgraph)
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::AlgebraicSystem" << endl;
#endif

  numsys   = anumsys;
  numgraph = anumgraph;
  numrhs   = 1;
  numbases = 0;
  numpart  = 6;
  offset   = 16;
  count    = 0;
}
  
AlgebraicSystem :: ~AlgebraicSystem()
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::~AlgebraicSystem" << endl;
#endif
}

void AlgebraicSystem :: CreateGraph(Integer * size)
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::CreateGraph" << endl;
#endif
  cout << "PILES: numnode;" << size[0] << endl;  
  Integer i;

  for (i=0; i<numgraph; i++)
    {
      graph[i] = new BaseGraph(size[i]);
    }
}

void AlgebraicSystem :: CreateMatrix(Integer * typesys, Integer * typemat, Integer * typegraph, 
				     Integer * dof, Integer * dir, Integer * con)
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::CreateMatrix" << endl;
#endif

  Integer i,j,size,nne;

  for (i=0; i<numsys; i++)
    {
      graph[typegraph[i]-1]->Create();

      nne  = graph[typegraph[i]-1]->GetNumNNE();
      size = graph[typegraph[i]-1]->GetSize();

      switch (typesys[i])
	{
	case 1:
	  /// real sparse 
	  RScalarCreate(i,size, nne, dir[i], dof[i], &typemat[i*numpart]);
	  break;

	case 2:
	  /// complex sparse
	  CScalarCreate(i,size, nne, dir[i], dof[i], &typemat[i*numpart]);
	  break;
	  
	case 3:
	  /// real block
	  RBlockCreate(i,size, nne, dir[i], dof[i], &typemat[i*numpart]);
	  break;

	case 4:
	  /// complex block
	  CBlockCreate(i,size, nne, dir[i], dof[i], &typemat[i*numpart]);
	  break;

	case 5:
	  /// real full
	  RFullCreate(i,size, nne, dir[i], dof[i], &typemat[i*numpart]);
	  break;

	case 6:
	  /// complex full
	  CFullCreate(i,size, nne, dir[i], dof[i], &typemat[i*numpart]);
	  break;

	default:
	  /// no such structure available
	  cout << "no such system type available" << endl;
	}
    }
}

void AlgebraicSystem :: AddGMGLevel()
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::AddGMGLevel" << endl;
#endif

  cout << "adding a new level" << endl;
  cout << "update the hierarchy - premat" << endl;
  cout << "update the solver" << endl;

}

///////////////////////////////////////////////////////////////////////////////////
/// private functions

void AlgebraicSystem :: RScalarCreate(Integer i, Integer size, Integer nne, 
				      Integer dir, Integer dof, Integer * typemat)
{
  Integer j;
  Integer * pos, rs;

  for (j=0; j<numpart; j++)
    {
      if (typemat[j] != 0)
	{
	  sysmat[(typemat[j]-1)*offset+i] = new RScalarMatrix(size, nne, dir);
	}
    }

  if (param[i]->GetAuxiliaryMatrix() == 1)
    {
      auxmat[i] = new RScalarMatrix(size, nne, dir);
    }
  else
    {
      auxmat[i] = sysmat[i];
    }

  if (param[i]->GetSpectralMatrix() == 1)
    {
      spemat[i] = new RScalarMatrix(size, nne, dir);

      if (param[i]->GetAuxiliaryMatrix() == 0)
	{
	  auxmat[i] = spemat[i];
	}
    }
  else
    {
      spemat[i] = sysmat[i];
    }

  if (param[i]->GetOutBackMemory() == 0)
    {
      rhs[i] = new RealVector(size, numrhs, dof);
      sol[i] = new RealVector(size, numrhs, dof);
    }
  else
    {
      rhs[i] = new RealVector(size, numrhs, dof, FALSE);
      sol[i] = new RealVector(size, numrhs, dof, FALSE);
    }


  switch (param[i]->GetPrecond())
    {
    case 1:
      premat[i] = new RIDScalarPrecond(*param[i], size, numrhs);
      break;
      
    case 2:
      premat[i] = new RJACScalarPrecond(*param[i], size, numrhs);
      break;
      
    case 3:
      premat[i] = new RSSORScalarPrecond(*param[i], size, numrhs);
      break;
      
    case 4:
      premat[i] = new RScalarMGPrecond(*param[i], size, numrhs);
      break;
      
    default:
      premat[i] = NULL;
      cout << "no preconditioner available" << endl;
    }
  
  switch (param[i]->GetSolver())
    {
    case 1:
      scheme[i] = new RRSolver(*param[i], size, numrhs, dof);
      break;
      
    case 2:
      scheme[i] = new CGSolver(*param[i], size, numrhs, dof);
      break;
      
    default:
      scheme[i] = NULL;
      cout << "no iteration scheme available" << endl;
    }

  // set the graph for the system matrix
  
  for (j=1; j<=size; j++)
    {
      pos = graph[0]->GetGraphRow(j);
      rs  = graph[0]->GetRowSize(j);
      
      sysmat[i]->SetGraphRow(rs, pos, j);
      
      if (param[i]->GetSpectralMatrix() == 1)
	{
	  spemat[i]->SetGraphRow(rs, pos, j);
	}
    } 
}

void AlgebraicSystem :: CScalarCreate(Integer i, Integer size, Integer nne, 
				      Integer dir, Integer dof, Integer * typemat)
{
  Integer j;
  Integer * pos, rs;

  for (j=0; j<numpart; j++)
    {
      if (typemat[j] != 0)
	{
	  sysmat[(typemat[j]-1)*offset+i] = new CScalarMatrix(size, nne, dir);
	}
    }

  auxmat[i] = new RScalarMatrix(size, nne, dir);

  if (param[i]->GetSpectralMatrix() == 1)
    {
      spemat[i] = new CScalarMatrix(size, nne, dir);
    }
  else
    {
      spemat[i] = sysmat[i];
    }

  if (param[i]->GetOutBackMemory() == 0)
    {
      rhs[i] = new ComplexVector(size, numrhs, dof);
      sol[i] = new ComplexVector(size, numrhs, dof);
    }
  else
    {
      rhs[i] = new ComplexVector(size, numrhs, dof, FALSE);
      sol[i] = new ComplexVector(size, numrhs, dof, FALSE);
    }

  switch (param[i]->GetPrecond())
    {
    case 1:
      //premat[i] = new CIDScalarPrecond(*param[i], size, numrhs);
      break;
      
    case 2:
      //premat[i] = new CJACScalarPrecond(*param[i], size, numrhs);
      ;
      break;
      
    case 3:
      //premat[i] = new CSSORScalarPrecond(*param[i], size, numrhs);
      ;
      break;
      
    case 4:
      premat[i] = new CScalarMGPrecond(*param[i], size, numrhs);
      break;
      
    default:
      premat[i] = NULL;
      cout << "no preconditioner available" << endl;
    }
  
  switch (param[i]->GetSolver())
    {
    case 1:
      //scheme[i] = new CRSolver(*param[i], size, numrhs, dof);
      ;
      break;
      
    case 2:
      scheme[i] = new QMRSolver(*param[i], size, numrhs, dof);
      break;
      
    default:
      scheme[i] = NULL;
      cout << "no iteration scheme available" << endl;
    }

  // set the graph for the system matrix
  
  for (j=1; j<=size; j++)
    {
      pos = graph[0]->GetGraphRow(j);
      rs  = graph[0]->GetRowSize(j);
      
      sysmat[i]->SetGraphRow(rs, pos, j);
      
      if (param[i]->GetSpectralMatrix() == 1)
	{
	  spemat[i]->SetGraphRow(rs, pos, j);
	}
    } 
}

void AlgebraicSystem :: RBlockCreate(Integer i, Integer size, Integer nne, 
				     Integer dir, Integer dof, Integer * typemat)
{
  Integer j;
  Integer * pos, rs;

  for (j=0; j<numpart; j++)
    {
      if (typemat[j] != 0)
	{
	  sysmat[(typemat[j]-1)*offset+i] = new RBlockMatrix(size, nne, dof, dir);
	}
    }

  auxmat[i] = new RScalarMatrix(size, nne, 0);

  if (param[i]->GetSpectralMatrix() == 1)
    {
      spemat[i] = new RBlockMatrix(size, nne, dof, dir);
    }
  else
    {
      spemat[i] = sysmat[i];
    }

  if (param[i]->GetOutBackMemory() == 0)
    {
      rhs[i]    = new RealVector(size, numrhs, dof);
      sol[i]    = new RealVector(size, numrhs, dof);
    }
  else
    {
      rhs[i]    = new RealVector(size, numrhs, dof, FALSE);
      sol[i]    = new RealVector(size, numrhs, dof, FALSE);
    }

  switch (param[i]->GetPrecond())
    {
    case 1:
      ;//premat[i] = new RIDBlockPrecond(*param[i], size, numrhs);
      break;
      
    case 2:
      ;//premat[i] = new RJACBlockPrecond(*param[i], size, numrhs);
      break;
      
    case 3:
      ;//premat[i] = new RSSORBlockPrecond(*param[i], size, numrhs);
      break;
      
    case 4:
      premat[i] = new RBlockMGPrecond(*param[i], size, numrhs, dof);
      break;
      
    default:
      premat[i] = NULL;
      cout << "no preconditioner available" << endl;
    }
  
  switch (param[i]->GetSolver())
    {
    case 1:
      scheme[i] = new RRSolver(*param[i], size, numrhs, dof);
      break;
      
    case 2:
      scheme[i] = new CGSolver(*param[i], size, numrhs, dof);
      break;
      
    default:
      scheme[i] = NULL;
      cout << "no iteration scheme available" << endl;
    }  

  // set the graph for the system matrix
      
  for (j=1; j<=size; j++)
    {
      pos = graph[0]->GetGraphRow(j);
      rs  = graph[0]->GetRowSize(j);
      
      sysmat[i]->SetGraphRow(rs, pos, j);
      
      if (param[i]->GetSpectralMatrix() == 1)
	{
	  spemat[i]->SetGraphRow(rs, pos, j);
	}
    } 
}

void AlgebraicSystem :: CBlockCreate(Integer i, Integer size, Integer nne, 
				     Integer dir, Integer dof, Integer * typemat)
{

}

void AlgebraicSystem :: RFullCreate(Integer i, Integer size, Integer nne, 
				    Integer dir, Integer dof, Integer * typemat)
{
  Integer j;
  Integer * pos, rs;

  for (j=0; j<numpart; j++)
    {
      if (typemat[j] != 0 && param[i]->GetOutBackMemory() == 0)
	{
	  sysmat[(typemat[j]-1)*offset+i] = new RFullMatrix(size, dir);
	}
      else if (typemat[j] != 0 && param[i]->GetOutBackMemory() == 1)
	{
	  sysmat[(typemat[j]-1)*offset+i] = new RFullMatrix(size, dir, FALSE);
	}	
    }

  auxmat[i] = new RScalarMatrix(size, nne, dir);

  if (param[i]->GetSpectralMatrix() == 1)
    {
      spemat[i] = new RFullMatrix(size, dir);
    }
  else
    {
      spemat[i] = sysmat[i];
    }

  if (param[i]->GetOutBackMemory() == 0)
    {
      rhs[i]    = new RealVector(size, numrhs, dof);
      sol[i]    = new RealVector(size, numrhs, dof);
    }
  else
    {
      rhs[i]    = new RealVector(size, numrhs, dof, FALSE);
      sol[i]    = new RealVector(size, numrhs, dof, FALSE);
    }

  switch (param[i]->GetPrecond())
    {
    case 1:
      ;//premat[i] = new RIDBlockPrecond(*param[i], size, numrhs);
      break;
      
    case 2:
      ;//premat[i] = new RJACBlockPrecond(*param[i], size, numrhs);
      break;
      
    case 3:
      ;//premat[i] = new RSSORBlockPrecond(*param[i], size, numrhs);
      break;
      
    case 4:
      premat[i] = new RFullMGPrecond(*param[i], size, numrhs);
      break;
      
    default:
      premat[i] = NULL;
      cout << "no preconditioner available" << endl;
    }

  switch (param[i]->GetSolver())
    {
    case 1:
      scheme[i] = new RRSolver(*param[i], size, numrhs, dof);
      break;
      
    case 2:
      scheme[i] = new CGSolver(*param[i], size, numrhs, dof);
      break;
      
    default:
      scheme[i] = NULL;
      cout << "no iteration scheme available" << endl;
    }  

  // set the graph for the system matrix
  
  for (j=1; j<=size; j++)
    {
      pos = graph[0]->GetGraphRow(j);
      rs  = graph[0]->GetRowSize(j);
      
      auxmat[i]->SetGraphRow(rs, pos, j);
    } 
}

void AlgebraicSystem :: CFullCreate(Integer i, Integer size, Integer nne, 
				    Integer dir, Integer dof, Integer * typemat)
{

}

}
