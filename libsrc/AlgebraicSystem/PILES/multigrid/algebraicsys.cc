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

AlgebraicSystem :: AlgebraicSystem(Integer anumsys, Integer anumgraph)
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::AlgebraicSystem" << endl;
#endif

  numsys   = anumsys;
  numgraph = anumgraph;
  numrhs   = 1;
  numbases = 0;
  numpart  = 5;
  blocksize= 4;
  offset   = blocksize*blocksize;
  count    = 0;

  for (Integer i=0;i<numpart*offset; i++)
    {
      sysmat[i] = NULL;
    }

#ifdef MEMTRACE
  sumdmem = 0;
  sumimem = 0;
#endif
}
  
AlgebraicSystem :: ~AlgebraicSystem()
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::~AlgebraicSystem" << endl;
#endif

  if (sysmat[0] != NULL)
    {
      delete sysmat[0];
      sysmat[0] = NULL;
    }

  if (spemat[0] != NULL)
    {
      delete spemat[0];
      spemat[0] = NULL;
    }

  if (auxmat[0] != NULL)
    {
      delete auxmat[0];
      auxmat[0] = NULL;
    }

  delete graph[0];

  delete premat[0];
  delete param[0];
  delete scheme[0];
  delete sol[0];
  delete rhs[0];
}

void AlgebraicSystem :: CreateGlobal2Local(Integer aglobsize)
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::CreateGlobal2Local" << endl;
#endif

  Integer i;

  globsize = aglobsize;
  glob2loc = new Integer[globsize*numsys];

  for (i=0; i<globsize*numsys; i++)
    {
      glob2loc[i] = 0;
    }
}


void AlgebraicSystem :: CreateMatrix(Integer matrix_row, Integer matrix_col, Integer * matrix_type, 
				     Integer matrix_class, Integer matrix_graph, Integer adof, Integer dir, 
				     Integer con)
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::CreateMatrix" << endl;
#endif

  Integer ind;

  ind = (matrix_row-1)*blocksize+matrix_col-1;

  // sort the graph array
  graph[ind]->Create();
  
  nne  = graph[ind]->GetNumNNE();
  size = graph[ind]->GetSize();
  dof  = adof;
  numrhs = 1;
  
  switch (matrix_class)
    {
    case 1:
      /// real sparse 
      RScalarCreate(matrix_row, ind, size, nne, dir, con, matrix_type);
      break;
      
    case 2:
      /// complex sparse
      CScalarCreate(ind, size, nne, dir, con, matrix_type);
      break;
      
    case 3:
      /// real block
      RBlockCreate(matrix_row, ind, size, nne, dir, dof, con, matrix_type);
      break;
      
    case 4:
      /// complex block
      CBlockCreate(ind, size, nne, dir, dof, con, matrix_type);
      break;
      
    case 5:
      /// real full
      RFullCreate(ind, size, nne, dir, dof, con, matrix_type);
      break;
      
    case 6:
      /// complex full
      CFullCreate(ind, size, nne, dir, dof, con, matrix_type);
      break;
      
    default:
      /// no such structure available
      (*cla) << "no such system type available" << endl;
    }
}

void AlgebraicSystem :: CreateSolver(Integer sys_id)
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::CreateSolver" << endl;
#endif

  switch (param[sys_id-1]->GetSolver())
    {
    case 1:
      scheme[sys_id-1] = new RRSolver(*param[sys_id-1], size, numrhs, dof);
      break;
      
    case 2:
      scheme[sys_id-1] = new CGSolver(*param[sys_id-1], size, numrhs, dof);
      break;
      
    case 3:
      scheme[sys_id-1] = new LanczosSolver(*param[sys_id-1], size, numrhs, dof);
      break; 
      
    default:
      scheme[sys_id-1] = NULL;
      (*cla) << "no iteration scheme available" << endl;
    }  
}

void AlgebraicSystem :: CreatePrecond(Integer matrix_class, Integer sys_id)
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::CreatePrecond" << endl;
#endif

  switch (matrix_class)
    {
    case 1:
      switch (param[sys_id-1]->GetPrecond())
	{
	case 1:
	  premat[sys_id-1] = new RIDScalarPrecond(*param[sys_id-1], size, numrhs);
	  break;
	  
	case 2:
	  premat[sys_id-1] = new RScalarMGPrecond(*param[sys_id-1], size, numrhs);
	  break;
	  
	default:
	  premat[sys_id-1] = NULL;
	  (*cla) << "no preconditioner available" << endl;
	}
      break;

    case 2:
      switch (param[sys_id-1]->GetPrecond())
	{
	case 1:
	  ;//premat[sys_id-1] = new CIDScalarPrecond(*param[sys_id-1], size, numrhs);
	  break;
	  
	case 2:
	  premat[sys_id-1] = new CScalarMGPrecond(*param[sys_id-1], size, numrhs);
	  break;
	  
	default:
	  premat[sys_id-1] = NULL;
	  (*cla) << "no preconditioner available" << endl;
	}
      break; 

    case 3:
      switch (param[sys_id-1]->GetPrecond())
	{
	case 1:
	  ;//premat[sys_id-1] = new RIDBlockPrecond(*param[sys_id-1], size, numrhs);
	  break;
	  
	case 2:
	  premat[sys_id-1] = new RBlockMGPrecond(*param[sys_id-1], size, numrhs, dof);
	  break;
	  
	default:
	  premat[sys_id-1] = NULL;
	  (*cla) << "no preconditioner available" << endl;
	}
      break;

    case 4:
      switch (param[sys_id-1]->GetPrecond())
	{
	case 1:
	  ;//premat[sys_id-1] = new CIDBlockPrecond(*param[sys_id-1], size, numrhs);
	  break;
	  
	case 2:
	  ;//premat[sys_id-1] = new CBlockMGPrecond(*param[sys_id-1], size, numrhs, dof);
	  break;
	  
	default:
	  premat[sys_id-1] = NULL;
	  (*cla) << "no preconditioner available" << endl;
	}
      break;

    case 5:
      switch (param[sys_id-1]->GetPrecond())
	{
	case 1:
	  ;//premat[sys_id-1] = new RIDFullPrecond(*param[sys_id-1], size, numrhs, dof);
	  break;
	  
	case 2:
	  ;//premat[sys_id-1] = new RFullMGPrecond(*param[sys_id-1], size, numrhs, dof);
	  break;
	  
	default:
	  premat[sys_id-1] = NULL;
	  (*cla) << "no preconditioner available" << endl;
	}
      break;

    case 6:
      switch (param[sys_id-1]->GetPrecond())
	{
	case 1:
	  ;//premat[sys_id-1] = new CIDFullPrecond(*param[sys_id-1], size, numrhs, dof);
	  break;
	  
	case 2:
	  ;//premat[sys_id-1] = new CFullMGPrecond(*param[sys_id-1], size, numrhs, dof);
	  break;
	  
	default:
	  premat[sys_id-1] = NULL;
	  (*cla) << "no preconditioner available" << endl;
	}
      break;

    case 7:
      switch (param[sys_id-1]->GetPrecond())
	{
	case 1:
	  ;//premat[sys_id-1] = new IDMixedPrecond(*param[sys_id-1], size, numrhs, dof);
	  break;
	  
	case 2:
	  ;//premat[sys_id-1] = new MixedMGPrecond(*param[sys_id-1], size, numrhs, dof);
	  break;
	  
	default:
	  premat[sys_id-1] = NULL;
	  (*cla) << "no preconditioner available" << endl;
	}
      break;

    default:
      (*cla) << "no matrix class available" << endl;
    }
}

void AlgebraicSystem :: ConstructEffectiveMatrix(Integer matrix_row, Integer matrix_col, Double * matrix_fac)
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::ConstructEffectiveMatrix" << endl;
#endif

  Integer ind;
  BaseMatrix *sys, * mat[4];

  ind = (matrix_row-1)*blocksize+matrix_col-1;

  sys    = sysmat[ind];
  mat[0] = sysmat[1*offset+ind];
  mat[1] = sysmat[2*offset+ind];
  mat[2] = sysmat[3*offset+ind];
  mat[3] = sysmat[4*offset+ind];

  sys->ConstructEffectiveMatrix(mat, matrix_fac);
}

void AlgebraicSystem :: AddGMGLevel()
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::AddGMGLevel" << endl;
#endif

  (*cla) << "adding a new level" << endl;
  (*cla) << "update the hierarchy - premat" << endl;
  (*cla) << "update the solver" << endl;

}

///////////////////////////////////////////////////////////////////////////////////
/// private functions

void AlgebraicSystem :: RScalarCreate(Integer sys_id, Integer ind, Integer size, Integer nne, 
				      Integer dir, Integer con, Integer * matrix_type)
{
#ifdef TRACE
  (*trace) << "entering AlgebraicSystem::RScalarCreate" << endl;
#endif
  Integer i,j;
  Integer * pos, rs;

  for (j=0; j<numpart; j++)
    {
      if (matrix_type[j] != 0)
	{
	  sysmat[(matrix_type[j]-1)*offset+ind] = new RScalarMatrix(size, nne, dir);

	  /// graph for the system matrix

	  for (i=1; i<=size; i++)
	    {
	      pos = graph[ind]->GetGraphRow(i);
	      rs  = graph[ind]->GetRowSize(i);
	      
	      sysmat[(matrix_type[j]-1)*offset+ind]->SetGraphRow(rs, pos, i);
	    } 
	}
      else
	{
	  sysmat[(matrix_type[j]-1)*offset+ind] = NULL;
	}
    }

  if ((ind % blocksize - (sys_id-1)) == 0) // diagonal entry
    {
      if (param[sys_id-1]->GetAuxiliaryMatrix() == 1)
	{
	  auxmat[sys_id-1] = new RScalarMatrix(size, nne, dir);
	}
      else
	{
	  auxmat[sys_id-1] = sysmat[(sys_id-1)*blocksize+sys_id-1];
	}
      
      if (param[sys_id-1]->GetSpectralMatrix() == 1)
	{
	  spemat[sys_id-1] = new RScalarMatrix(size, nne, dir);
	  
	  if (param[sys_id-1]->GetAuxiliaryMatrix() == 0)
	    {
	      auxmat[sys_id-1] = spemat[sys_id-1];
	    }

	  for (j=1; j<=size; j++)
	    {
	      pos = graph[ind]->GetGraphRow(j);
	      rs  = graph[ind]->GetRowSize(j);
	      
	      spemat[sys_id-1]->SetGraphRow(rs, pos, j);
	      
	      if (param[sys_id-1]->GetAuxiliaryMatrix() == 0)
		{
		  auxmat[sys_id-1]->SetGraphRow(rs, pos, j);
		}
	    } 
	}
      else
	{
	  spemat[sys_id-1] = sysmat[(sys_id-1)*blocksize+sys_id-1];
	}
      
      if (param[sys_id-1]->GetOutBackMemory() == 0)
	{
	  rhs[sys_id-1] = new RealVector(size, numrhs, 1);
	  sol[sys_id-1] = new RealVector(size, numrhs, 1);
	}
      else
	{
	  rhs[sys_id-1] = new RealVector(size, numrhs, 1, FALSE);
	  sol[sys_id-1] = new RealVector(size, numrhs, 1, FALSE);
	}
    }
}

void AlgebraicSystem :: CScalarCreate(Integer ind, Integer size, Integer nne, 
				      Integer dir, Integer con, Integer * matrix_type)
{
//   Integer j;
//   Integer * pos, rs;

//   for (j=0; j<numpart; j++)
//     {
//       if (typemat[j] != 0)
// 	{
// 	  sysmat[(typemat[j]-1)*offset+i] = new CScalarMatrix(size, nne, dir);
// 	}
//     }

//   auxmat[i] = new RScalarMatrix(size, nne, dir);

//   if (param[i]->GetSpectralMatrix() == 1)
//     {
//       spemat[i] = new CScalarMatrix(size, nne, dir);
//     }
//   else
//     {
//       spemat[i] = sysmat[i];
//     }

//   if (param[i]->GetOutBackMemory() == 0)
//     {
//       rhs[i] = new ComplexVector(size, numrhs, dof);
//       sol[i] = new ComplexVector(size, numrhs, dof);
//     }
//   else
//     {
//       rhs[i] = new ComplexVector(size, numrhs, dof, FALSE);
//       sol[i] = new ComplexVector(size, numrhs, dof, FALSE);
//     }

//   // set the graph for the system matrix
  
//   for (j=1; j<=size; j++)
//     {
//       pos = graph[0]->GetGraphRow(j);
//       rs  = graph[0]->GetRowSize(j);
      
//       sysmat[i]->SetGraphRow(rs, pos, j);
      
//       if (param[i]->GetSpectralMatrix() == 1)
// 	{
// 	  spemat[i]->SetGraphRow(rs, pos, j);
// 	}
//     } 
}

void AlgebraicSystem :: RBlockCreate(Integer sys_id, Integer ind, Integer size, Integer nne, 
				     Integer dir, Integer dof, Integer con, Integer * matrix_type)
{
  Integer i,j;
  Integer * pos, rs;
  Boolean setgraph = FALSE;

  for (j=0; j<numpart; j++)
    {
      if (matrix_type[j] != 0)
	{
	  sysmat[(matrix_type[j]-1)*offset+ind] = new RBlockMatrix(size, nne, dof, dir);

	  /// graph for the system matrix

	  if (!setgraph)
	    {
	      for (i=1; i<=size; i++)
		{
		  pos = graph[ind]->GetGraphRow(i);
		  rs  = graph[ind]->GetRowSize(i);
		  
		  sysmat[(matrix_type[j]-1)*offset+ind]->SetGraphRow(rs, pos, i);
		} 

	      setgraph = TRUE;
	    }
	  else
	    {
	      sysmat[(matrix_type[j]-1)*offset+ind]->SetGraphPointer(sysmat[ind]);
	    }
	}
      else
	{
	  sysmat[(matrix_type[j]-1)*offset+ind] = NULL;
	}
    }

  if ((ind % blocksize - (sys_id-1)) == 0) // diagonal entry
    {
      auxmat[sys_id-1] = new RScalarMatrix(size, nne, dir);
      
      if (param[sys_id-1]->GetSpectralMatrix() == 1)
	{
	  spemat[sys_id-1] = new RBlockMatrix(size, nne, dof, dir);

	  spemat[sys_id-1]->SetGraphPointer(sysmat[ind]);

// 	  for (j=1; j<=size; j++)
// 	    {
// 	      pos = graph[ind]->GetGraphRow(j);
// 	      rs  = graph[ind]->GetRowSize(j);
	      
// 	      spemat[sys_id-1]->SetGraphRow(rs, pos, j);
// 	    } 
	}
      else
	{
	  spemat[sys_id-1] = sysmat[(sys_id-1)*blocksize+sys_id-1];
	}

      if (param[sys_id-1]->GetOutBackMemory() == 0)
	{
	  rhs[sys_id-1] = new RealVector(size, numrhs, dof);
	  sol[sys_id-1] = new RealVector(size, numrhs, dof);
	}
      else
	{
	  rhs[sys_id-1] = new RealVector(size, numrhs, dof, FALSE);
	  sol[sys_id-1] = new RealVector(size, numrhs, dof, FALSE);
	}
    }
}

void AlgebraicSystem :: CBlockCreate(Integer ind, Integer size, Integer nne, 
				     Integer dir, Integer dof, Integer con, Integer * matrix_type)
{

}

void AlgebraicSystem :: RFullCreate(Integer ind, Integer size, Integer nne, 
				    Integer dir, Integer dof, Integer con, Integer * matrix_type)
{
//   Integer j;
//   Integer * pos, rs;

//   for (j=0; j<numpart; j++)
//     {
//       if (typemat[j] != 0 && param[i]->GetOutBackMemory() == 0)
// 	{
// 	  sysmat[(typemat[j]-1)*offset+i] = new RFullMatrix(size, dir);
// 	}
//       else if (typemat[j] != 0 && param[i]->GetOutBackMemory() == 1)
// 	{
// 	  sysmat[(typemat[j]-1)*offset+i] = new RFullMatrix(size, dir, FALSE);
// 	}	
//     }

//   auxmat[i] = new RScalarMatrix(size, nne, dir);

//   if (param[i]->GetSpectralMatrix() == 1)
//     {
//       spemat[i] = new RFullMatrix(size, dir);
//     }
//   else
//     {
//       spemat[i] = sysmat[i];
//     }

//   if (param[i]->GetOutBackMemory() == 0)
//     {
//       rhs[i]    = new RealVector(size, numrhs, dof);
//       sol[i]    = new RealVector(size, numrhs, dof);
//     }
//   else
//     {
//       rhs[i]    = new RealVector(size, numrhs, dof, FALSE);
//       sol[i]    = new RealVector(size, numrhs, dof, FALSE);
//     }

//   switch (param[i]->GetPrecond())
//     {
//     case 1:
//       ;//premat[i] = new RIDBlockPrecond(*param[i], size, numrhs);
//       break;
      
//     case 2:
//       ;//premat[i] = new RJACBlockPrecond(*param[i], size, numrhs);
//       break;
      
//     case 3:
//       ;//premat[i] = new RSSORBlockPrecond(*param[i], size, numrhs);
//       break;
      
//     case 4:
//       premat[i] = new RFullMGPrecond(*param[i], size, numrhs);
//       break;
      
//     default:
//       premat[i] = NULL;
//       (*cla) << "no preconditioner available" << endl;
//     }

//   // set the graph for the system matrix
  
//   for (j=1; j<=size; j++)
//     {
//       pos = graph[0]->GetGraphRow(j);
//       rs  = graph[0]->GetRowSize(j);
      
//       auxmat[i]->SetGraphRow(rs, pos, j);
//     } 
}

void AlgebraicSystem :: CFullCreate(Integer ind, Integer size, Integer nne, 
				    Integer dir, Integer dof, Integer con, Integer * matrix_type)
{

}

}
