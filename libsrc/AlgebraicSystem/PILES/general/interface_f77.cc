#include <fstream.h>

#include <general.hh>
#include <multigrid.hh>

using namespace CoupledField;

/////////////////////////////////////////////////////////////////////////
/// important definitions
  
enum GraphType {NOGRAPH, NODE, EDGE, FACE, VOLUME};
enum MatrixType {NOTYPE, SYSTEM_MAT, STIFFNESS_MAT, DAMPING_MAT, CONVECTION_MAT, MASS_MAT};
enum MatrixClass {NOCLASS, RSPARSE, CSPARSE, RBLOCK, CBLOCK, RFULL, CFULL, MIXED};
enum PrecondType {NOPRECOND, ID, JACOBI, SSOR, MG};
enum SolverType {NOSOLVER, RICHARDSON, CG};

Integer dir[1],con[1],type[1],graph[1],dof[1],size[1],typemat[6];
Integer up, count = 0, matrix = 1;
Integer pos[50];
Integer in, cn;

Double cputime;
CPUClock cpuclock;

AlgebraicSystem * algsys;
FileSystem filesys("capa-test");

extern "C"
{
void setsize_(int *n, int *bsize, int *esize, int *dsize, int *csize,
              int *mtype, int *sym, int *def, int *dm,  int p[], double q[])
{
  algsys = new BlockSystem(1, 1);

  dof[0]   = *bsize;
  con[0]   = *csize;
  dir[0]   = *dsize;
  size[0]  = *n;
  graph[0] = NODE;

  typemat[0] = SYSTEM_MAT; 
  typemat[1] = NOTYPE; 
  typemat[2] = NOTYPE; 
  typemat[3] = NOTYPE;  
  typemat[4] = NOTYPE; 
  typemat[5] = NOTYPE;  

  switch (*mtype)
    {
    case 1:
      type[0] = RSPARSE;
      break;

    case 2:
      type[0] = RBLOCK;
      break;

    case 3:
      type[0] = RFULL;
      break;

    default:
      type[0] = RSPARSE;
    }

  algsys->CreateParameter();

  algsys->SetOutBackMemory(matrix);
  algsys->SetAccuracy(q[0],matrix);
  //algsys->SetEpsMat(q[3],matrix);
  algsys->SetEpsMat(0,matrix);
  algsys->SetAlpha(q[4],matrix);
  //algsys->SetCoarseSystem(p[3],matrix);
  algsys->SetCoarseSystem(1000,matrix);
  //algsys->SetMaxNumIter(p[5],matrix);
  algsys->SetMaxNumIter(100,matrix);
  //algsys->SetCycle(p[11],matrix);
  algsys->SetCycle(1,matrix);
  algsys->SetSmoothType(p[12],matrix);

  algsys->SetSmoothStepFor(1,matrix);
  algsys->SetSmoothStepBack(1,matrix);

//   algsys->SetSmoothStepFor(p[13],matrix);
//   algsys->SetSmoothStepBack(p[13],matrix);

  algsys->SetPrecond(MG,matrix);
  algsys->SetSolver(CG,matrix);

  if (*dm == 1)
    {
      algsys->SetSpectralMatrix(matrix);
    }

  algsys->CreateGraph(size);

  cout << "### after create graph " << endl;
  cout << "size         " << size[0] << endl;
  cout << "dof per node " << dof[0] << endl;
  cout << "elements     " << *esize << endl;
  cout << "constraints  " << con[0] << endl;
  cout << "dirichlet    " << dir[0] << endl;
}

void setelement_(double v[], int p[], int * q, int * es, int * ep, int * ns)
{
  if (dof[0] != 1)
    {
      Integer i,k;
      
      k = 0;
      
      for (i=0; i<*es; i++)
	{
	  pos[i] = Integer((p[k]-1)/dof[0])+1;
	  k += dof[0];
	}
 
      algsys->SetElementMatrix(v, pos, *es, 0, matrix);
    }
  else
    {
      algsys->SetElementMatrix(v, p, *es, 0, matrix);
    }
}

void setconnection_(int p[], int * q, int * es)
{
  if (dof[0] != 1)
    {
      Integer i,k;
      
      k = 0;
      
      for (i=0; i<*es; i++)
	{
	  pos[i] = Integer((p[k]-1)/dof[0])+1;
	  k += dof[0];
	}
      
      algsys->SetElementPos(pos, *es, matrix);
    }
  else
    {
      algsys->SetElementPos(p, *es, matrix);
    }
}

void creatematrix_()
{
  //algsys->InitRHS(matrix);
  //algsys->InitSol(matrix);
  //cout << "before init matrix " << endl;
  
  cputime = cpuclock.GetTime();

  if (count == 0)
    {
      algsys->CreateMatrix(type, typemat, graph, dof, dir, con);
    }

  algsys->InitMatrix(matrix, 0);

  count++;

  cputime = cpuclock.GetTime() - cputime;

  cout << "### create the matrix structure (memory) " << cputime << " sec" << endl;
}

void constructmatrixhierarchy_()
{
  cputime = cpuclock.GetTime();

  algsys->CalcPrecond(up, matrix);
  up = 0;

  cputime = cpuclock.GetTime() - cputime;

  cout << "### construct matrix hierarchy " << cputime << " sec" << endl;
}

void solvediffrhs_(double f[], double u[])
{
  cputime = cpuclock.GetTime();

  algsys->Solve(f, u);

  cputime = cpuclock.GetTime() - cputime;

  cout << "### solution phase " << cputime << " sec" << endl;
}

void solvediffmat_(double f[], double u[])
{
  cputime = cpuclock.GetTime();

  algsys->Solve(f, u);

  cputime = cpuclock.GetTime() - cputime;

  cout << "### solution phase " << cputime << " sec" << endl;
}

int getnumbiter_()
{
  return algsys->GetNumIter(matrix);
}

void setcount_()
{
  up = 1;
  algsys->SetCounter();
}

void setdirichlet_(int * i, int * ind, int * comp, double *value)
{
  if (dof[0] != 1)
    {
      //cout << "dirichlet " << *i << " " << *ind << " " << *comp << " " << *value << endl;

      in = Integer((*ind-1)/dof[0]+1);
      cn = ((*ind-1)%dof[0])+1;
    }
  else
    {
      in = *ind;
      cn = *comp;
    }

  //cout << "dirichlet " << *i << " " << in << " " << *value << " " << cn << endl;

  algsys->SetDirichlet(*i, in, *value, cn, matrix);
}

void updatedirichlet_(int *i, double *val)
{
  algsys->UpdateDirichlet(*i, *val, matrix);
}

void setconstraint_(int cm[], double cf[])
{
  algsys->SetConstraint(cm, cf, matrix);
}

void setmatrixpointer_(double * v, int * s, int * p)
{
  algsys->SetMatrixPointer(s, p, v, matrix);
}

void setauxelementmatrix_(double * v, int * p, int * es)
{
  algsys->SetAuxElementMatrix(v, p, *es, matrix);
}

void cleanup_()
{
  delete algsys;
}

void initprecond_(int *nrhs, int *nbases)
{
  //fe_code->InitPrecond(*nrhs,*nbases);
} 

void setkrylov_()
{
  //algsys->SetKrylov();
}

}


