#ifndef FILE_INTERFACE_PILES_2002
#define FILE_INTERFACE_PILES_2002

#include <abstractAlgSys.hh>

#ifdef GRIDLIB
{
using namespace std;
#include <general.hh>
#include <multigrid.hh>
}
#else
#include <general.hh>
#include <multigrid.hh>
#endif

namespace CoupledField
{

class AlgSysPILES: public AbstractAlgebraicSys
{
public:
  //!
  AlgSysPILES(){};

  //!
  virtual void InitAlgSys(Integer anumsys, Integer anumgraph)
  { algsys = new BlockSystem(anumsys,anumgraph);}

  //!
  virtual void SetSolverParameter(Integer nsys, Double eps, Double dampiter, Integer maxnumit,
                                  Integer solvertype, Integer precondtype, Integer numeqcoarse,
				  Double coarsealpha)
  {  
    nsys++;
    algsys->CreateParameter();
    algsys->SetAccuracy(eps,nsys);
    algsys->SetMaxNumIter(maxnumit,nsys);
    algsys->SetPrecond(precondtype,nsys);
    algsys->SetSolver(solvertype,nsys);
    algsys->SetDampIter(dampiter,nsys);
    algsys->SetCoarseSystem(numeqcoarse,nsys);
    algsys->SetAlpha(coarsealpha,nsys);
  }

  virtual void InitAlgSysGraph(Integer numnode, Integer matrix_row, Integer matrix_col, 
			       Integer matrix_graphtype)
  { 
    matrix_row ++;
    matrix_col ++;
    algsys->CreateGlobal2Local(numnode);
    algsys->CreateGraph(numnode,matrix_graphtype,matrix_row,matrix_col);
  }

  virtual void SetAlgSysGraph(Integer *pos, Integer elemsize, Integer fe_type, Integer matrix_row, 
			      Integer matrix_col)
  { 
    matrix_row ++;
    matrix_col ++;
    //    algsys->SetElementPos(pos,elemsize,2,matrix_row,matrix_col);
    algsys->SetElementPos(pos,elemsize,fe_type,matrix_row,matrix_col);
  }

  virtual void CreateAlgSysMatrices(Integer matrix_row, Integer matrix_col, Integer *matrixsystype, 
				    Integer matrixtype, Integer graphtype, Integer numdofpernode, 
				    Integer numdirichlets, Integer numconstraints)
  {
    matrix_row ++;
    matrix_col ++;
    algsys->CreateMatrix(matrix_row,matrix_col,matrixsystype,matrixtype,graphtype,numdofpernode,
			 numdirichlets,numconstraints);
    algsys->CreateSolver(matrix_row);
    algsys->CreatePrecond(matrixtype, matrix_row);
  }

  virtual void ResetAlgSys(Integer matrix_row, Integer matrix_col, Integer matrix_id)
  {
    matrix_row ++;
    matrix_col ++;
    algsys->InitRHS(matrix_row);
    algsys->InitSol(matrix_row);
    algsys->InitMatrix(matrix_row,matrix_col,matrix_id);
  }

  virtual void ResetRHS(Integer matrix_row)
  {
    matrix_row ++;
    algsys->InitRHS(matrix_row);
  }

  virtual void ResetMatrix(Integer matrix_row, Integer matrix_col, Integer matrix_id)
  {
    matrix_row ++;
    matrix_col ++;
    algsys->InitMatrix(matrix_row,matrix_col,matrix_id);
  }

  virtual void PutElemMatAlgSys(Double *elemmat, Integer *help, Integer numnodelem, Integer matrix_row,
				Integer matrix_col, Integer matrix_id)
  {
    matrix_row ++;
    matrix_col ++;
    algsys->SetElementMatrix(elemmat, help, numnodelem,  matrix_row, matrix_col, matrix_id);
  }

  //!
  virtual void ComputeSysMatrix(Integer matrix_row, Integer matrix_col, Double * matrix_factor)
  {
    matrix_row ++;
    matrix_col ++;
    algsys->ConstructEffectiveMatrix(matrix_row,matrix_col,matrix_factor);
  }

  //!
  virtual void SetBCDirichlet(Integer restrnum, Integer nodenr, Double restrval, Integer ndof, 
			      Integer matrix_row, Integer matrix_col, Integer matrix_id)
  {
    matrix_row ++;
    matrix_col ++;
    algsys->SetDirichlet(restrnum,nodenr,restrval,ndof,matrix_row,matrix_col,matrix_id);
  }

  //!
  virtual void SetBCDirichletUpdate(Integer restrnum, Integer nodenr, Double restrval, Integer ndof, 
			      Integer matrix_row, Integer matrix_col, Integer matrix_id)
  {
    matrix_row ++;
    matrix_col ++;
    algsys->UpdateDirichlet(restrnum,restrval,matrix_row,matrix_col,matrix_id);
  }

  //!
  virtual void UpdateRHS(Integer matrix_row, Integer matrix_col, Integer matrix_id, Double * vec)
  {
    matrix_row ++;
    matrix_col ++;
    algsys->UpdateRHS(matrix_row,matrix_col,matrix_id,vec);
  }

  //! Add value addval to rhs-vector at position pos
  virtual void AddRHS(Double addval, Integer pos, Integer sys_id)
  {
    sys_id++;
//    algsys->AddRHS(addval,pos,sys_id);
  }


  //!
  virtual void ComputePrecond(Integer job, Integer nsys)
  {
    nsys++;
    algsys->CalcPrecond(job,nsys);
  }

  //!
  virtual void SolveAlgSys(Integer nsys)
  {
    nsys++;
    algsys->Solve();
  }

  //!
  virtual Double * GetSolution(Integer nsys)
  {
    nsys++;
    return algsys->GetSolutionVal(nsys);
  }

  //!
  virtual Double CalcEnergyNorm(Integer matrix_row, Integer matrix_col, Integer matrix_id, Double * u)
  {
    matrix_row ++;
    matrix_col ++;
    return algsys->CalcEnergyNorm(matrix_row,matrix_col,matrix_id,u);
  }

  virtual ~AlgSysPILES()
  {
#ifdef TRACE
    (*trace) << "entering AlgSysPILES::~AlgSysPILES " << std::endl;
#endif
    if (algsys) delete algsys;
  }

private:
  AlgebraicSystem * algsys;

} ;

} // end of namespace

#endif
