#ifndef FILE_ALGEBRAICSYSTEM_PILES
#define FILE_ALGEBRAICSYSTEM_PILES

namespace CoupledField
{

class AlgebraicSystem
{
public:
  ///
  AlgebraicSystem(Integer anumsys, Integer anumgraph);

  ///
  virtual ~AlgebraicSystem();

  ///
  virtual void CalcPrecond() = 0;

  ///
  virtual void CalcPrecond(Integer newprecond, Integer mat) = 0;

  ///
  virtual void Solve() = 0;

  ///
  virtual void Solve(Double * f, Double * u) = 0;

  ///
  virtual void CreateParameter() = 0;

  ////////////////////////////////////////////////////////////////////////////

  ///
  void SetAccuracy(Double aeps, Integer mat)
    {param[mat-1]->SetAccuracy(aeps);}

  ///
  void SetMaxNumIter(Integer amaxnumiter, Integer mat)
    {param[mat-1]->SetMaxNumIter(amaxnumiter);}
  
  ///
  void SetPrecond(Integer aprecond, Integer mat)
    {param[mat-1]->SetPrecond(aprecond);}

  ///
  void SetSolver(Integer asolver, Integer mat)
    {param[mat-1]->SetSolver(asolver);}
  
  ///
  void SetDampIter(Double adampiter, Integer mat)
    {param[mat-1]->SetDampIter(adampiter);}

  ///
  void SetEpsMat(Double aepsmat, Integer mat)
    {param[mat-1]->SetEpsMat(aepsmat);};

  ///
  void SetAlpha(Double aalpha, Integer mat)
    {param[mat-1]->SetAlpha(aalpha);};

  ///
  void SetCoarseSystem(Integer acoarsesystem, Integer mat)
    {param[mat-1]->SetCoarseSystem(acoarsesystem);};
  
  ///
  void SetCycle(Integer acycle, Integer mat)
    {param[mat-1]->SetCycle(acycle);};

  ///
  void SetSmoothType(Integer asmoothingtype, Integer mat)
    {param[mat-1]->SetSmoothType(asmoothingtype);};

  ///
  void SetSmoothStepFor(Integer asmoothingfor, Integer mat)
    {param[mat-1]->SetSmoothStepFor(asmoothingfor);};

  ///
  void SetSmoothStepBack(Integer asmoothingback, Integer mat)
    {param[mat-1]->SetSmoothStepBack(asmoothingback);};

  ///
  void SetSpectralMatrix(Integer mat)
    {param[mat-1]->SetSpectralMatrix();}

  ///
  void SetPenalty(double apenalty, Integer mat)
    {param[mat-1]->SetPenalty(apenalty);};

  ///
  void SetOutBackMemory(Integer mat)
    {param[mat-1]->SetOutBackMemory();}

  /////////////////////////////////////////////////////////////////////////////

  ///
  void SetSimRHS(Integer anumrhs)
    {numrhs = anumrhs;};

  ///
  void SetKrylov(Integer anumbases)
    {numbases = anumbases;};

  ///
  void SetCounter(){count = 0;};

  /////////////////////////////////////////////////////////////////////////////////////////

  ///
  void SetElementMatrix(Double * elemmat, Integer * connect, Integer elemsize, 
			Integer type,Integer mat)
    {sysmat[type*offset+mat-1]->Assemble(elemmat, connect, elemsize);}

  ///
  void SetAuxElementMatrix(Double * elemmat, Integer * connect, Integer elemsize ,Integer mat)
    {auxmat[mat-1]->Assemble(elemmat, connect, elemsize);}

  ///
  void SetElementRHS(Double * elemrhs, Integer * connect, Integer elemsize, Integer mat)
    {rhs[mat-1]->Assemble(elemrhs, connect, elemsize);};

  ///
  void InitRHS(Integer mat)
    {rhs[mat-1]->Init();};

  void InitSol(Integer mat)
    {sol[mat-1]->Init();};

  ///
  void CreateMatrix(Integer * typesys, Integer * typemat, Integer * typegraph, Integer * dof, 
		    Integer * dir, Integer * con);

  ///
  void SetMatrixFactor(Double val, Integer type, Integer mat)
    {sysmat[type*offset+mat-1]->SetMatrixFactor(val);};

  ///
  void MultMat(Integer i, Integer j, Integer type, Double * u, Double * v)
    {sysmat[type*offset+(i-1)*numsys+j-1]->Mult(u,v);};

  ///
  void AddGMGLevel();

  ///
  void InitMatrix(Integer mat, Integer type)
    {sysmat[type*offset+mat-1]->InitMatrix();};

  ///
  void SetMatrixPointer(Integer * astart, Integer * apos, Double * aval, Integer mat)
    {sysmat[mat-1]->SetPointer(astart, apos, aval);};

  /////////////////////////////////////////////////////////////////////////////////////////
  
  ///
  void SetDirichlet(Integer i, Integer num, Double val, Integer comp, Integer mat)
    {sysmat[mat-1]->SetDirichlet(i, num, val, comp);};

  ///
  void UpdateDirichlet(Integer num, Double val, Integer mat)
    {sysmat[mat-1]->UpdateDirichlet(num, val);};

  ///
  void SetConstraint(Integer * cm, Double * cf, Integer mat)
    {sysmat[mat-1]->SetConstraint(cm, cf);};
  
  /////////////////////////////////////////////////////////////////////////////////////////

  ///
  void CreateGraph(Integer * size);

  ///
  void SetGraph(Integer * row, Integer * col, Integer size, Integer mat);

  ///
  void SetElementPos(Integer * connect, Integer elemsize, Integer mat)
  { graph[mat-1]->SetElementPos(connect, elemsize);};

  /////////////////////////////////////////////////////////////////////////////////////////

  ///
  Integer GetNumIter(Integer mat) const 
    {return scheme[mat-1]->GetNumIter();};

  ///
  Double * GetSolution(Integer mat) const {return sol[mat-1]->GetPointer();};

  ///
  void Print(Integer type, Integer mat) const {sysmat[type*offset+mat-1]->Print();};

  /////////////////////////////////////////////////////////////////////////////////////////
  

protected:

  ///
  BaseMatrix * sysmat[80], * spemat[4], * auxmat[16];

  ///
  BaseGraph * graph[4];

  ///
  BaseVector * rhs[4], * sol[4];

  ///
  BaseSolver * scheme[4];

  ///
  BasePrecond * premat[4];

  ///
  BaseParameter * param[4];

  ///
  Integer numsys, numgraph, numrhs, numbases, count, offset, numpart;

private:

  ///
  void RScalarCreate(Integer i, Integer size, Integer nne, Integer dir, Integer dof, Integer * typemat);

  ///
  void CScalarCreate(Integer i, Integer size, Integer nne, Integer dir, Integer dof, Integer * typemat);

  ///
  void RBlockCreate(Integer i, Integer size, Integer nne, Integer dir, Integer dof, Integer * typemat);

  ///
  void CBlockCreate(Integer i, Integer size, Integer nne, Integer dir, Integer dof, Integer * typemat);

  ///
  void RFullCreate(Integer i, Integer size, Integer nne, Integer dir, Integer dof, Integer * typemat);

  ///
  void CFullCreate(Integer i, Integer size, Integer nne, Integer dir, Integer dof, Integer * typemat);
};

}

#endif // FILE_ALGEBRAICSYSTEM_PILES
