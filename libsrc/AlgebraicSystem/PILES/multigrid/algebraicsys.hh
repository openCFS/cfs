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
  virtual void CalcPrecond(Integer newprecond, Integer sys_id) = 0;

  ///
  virtual void Solve() = 0;

  ///
  virtual void Solve(Double * f, Double * u) = 0;

  ///
  virtual void CreateParameter() = 0;

  ////////////////////////////////////////////////////////////////////////////

  ///
  void SetAccuracy(Double aeps, Integer sys_id)
    {param[sys_id-1]->SetAccuracy(aeps);}

  ///
  void SetMaxNumIter(Integer amaxnumiter, Integer sys_id)
    {param[sys_id-1]->SetMaxNumIter(amaxnumiter);}
  
  ///
  void SetPrecond(Integer aprecond, Integer sys_id)
    {param[sys_id-1]->SetPrecond(aprecond);}

  ///
  void SetSolver(Integer asolver, Integer sys_id)
    {param[sys_id-1]->SetSolver(asolver);}
  
  ///
  void SetDampIter(Double adampiter, Integer sys_id)
    {param[sys_id-1]->SetDampIter(adampiter);}

  ///
  void SetEpsMat(Double aepsmat, Integer sys_id)
    {param[sys_id-1]->SetEpsMat(aepsmat);};

  ///
  void SetAlpha(Double aalpha, Integer sys_id)
    {param[sys_id-1]->SetAlpha(aalpha);};

  ///
  void SetCoarseSystem(Integer acoarsesystem, Integer sys_id)
    {param[sys_id-1]->SetCoarseSystem(acoarsesystem);};
  
  ///
  void SetCycle(Integer acycle, Integer sys_id)
    {param[sys_id-1]->SetCycle(acycle);};

  ///
  void SetSmoothType(Integer asmoothingtype, Integer sys_id)
    {param[sys_id-1]->SetSmoothType(asmoothingtype);};

  ///
  void SetSmoothStepFor(Integer asmoothingfor, Integer sys_id)
    {param[sys_id-1]->SetSmoothStepFor(asmoothingfor);};

  ///
  void SetSmoothStepBack(Integer asmoothingback, Integer sys_id)
    {param[sys_id-1]->SetSmoothStepBack(asmoothingback);};

  ///
  void SetSpectralMatrix(Integer sys_id)
    {param[sys_id-1]->SetSpectralMatrix();}

  ///
  void SetPenalty(double apenalty, Integer sys_id)
    {param[sys_id-1]->SetPenalty(apenalty);};

  ///
  void SetOutBackMemory(Integer sys_id)
    {param[sys_id-1]->SetOutBackMemory();}

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
			Integer matrix_row, Integer matrix_col,Integer matrix_id)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->Assemble(elemmat, connect, elemsize);}

  ///
  void SetAuxElementMatrix(Double * elemmat, Integer * connect, Integer elemsize ,Integer sys_id)
    {auxmat[sys_id-1]->Assemble(elemmat, connect, elemsize);}

  ///
  void SetElementRHS(Double * elemrhs, Integer * connect, Integer elemsize, Integer sys_id)
    {rhs[sys_id-1]->Assemble(elemrhs, connect, elemsize);};

  ///
  void InitRHS(Integer sys_id)
    {rhs[sys_id-1]->Init();};

  void InitSol(Integer sys_id)
    {sol[sys_id-1]->Init();};

  ///
  void CreateMatrix(Integer matrix_row, Integer matrix_col, Integer * matrix_type, Integer matrix_class,
		    Integer graph, Integer adof, Integer dir, Integer con);

  ///
  void CreateSolver(Integer sys_id);

  ///
  void CreatePrecond(Integer matrix_class, Integer sys_id);

  ///
  void SetMatrixFactor(Double val, Integer matrix_row, Integer matrix_col, Integer matrix_id)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->SetMatrixFactor(val);};

  ///
  void MultMat(Integer matrix_row, Integer matrix_col, Integer matrix_id, Double * u, Double * v)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->Mult(u,v);};

  ///
  void AddGMGLevel();

  ///
  void InitMatrix(Integer matrix_row, Integer matrix_col, Integer matrix_id)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->InitMatrix();};

  ///
  void SetMatrixPointer(Integer * astart, Integer * apos, Double * aval, 
			Integer matrix_row, Integer matrix_col, Integer matrix_id)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->SetPointer(astart, apos, aval);};

  /////////////////////////////////////////////////////////////////////////////////////////
  
  ///
  void UpdateRHS(Integer matrix_row, Integer matrix_col, Integer matrix_id, Double * fup)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->MultAdd(fup, *rhs[matrix_row-1]);};

  ///
  void SetDirichlet(Integer i, Integer num, Double val, Integer comp, 
		    Integer matrix_row, Integer matrix_col, Integer matrix_id)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->SetDirichlet(i, num, val, comp);};

  ///
  void UpdateDirichlet(Integer num, Double val, Integer matrix_row, Integer matrix_col, Integer matrix_id)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->UpdateDirichlet(num, val);};

  ///
  void SetConstraint(Integer * cm, Double * cf, Integer matrix_row, Integer matrix_col, Integer matrix_id)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->SetConstraint(cm, cf);};
  
  /////////////////////////////////////////////////////////////////////////////////////////

  ///
  void CreateGlobal2Local(Integer aglobsize);

  ///
  void CreateGraph(Integer size, Integer graph_row, Integer graph_col)
    {graph[(graph_row-1)*blocksize+graph_col-1] = new BaseGraph(size);};

  ///
  void SetElementPos(Integer * connect, Integer elemsize, Integer matrix_row, Integer matrix_col)
    {graph[(matrix_row-1)*blocksize+matrix_col-1]->SetElementPos(connect, elemsize);};

  /////////////////////////////////////////////////////////////////////////////////////////

  ///
  void ConstructEffectiveMatrix(Integer matrix_row, Integer matrix_col, Double * matrix_fac);

  ///
  Integer GetNumIter(Integer sys_id) const 
    {return scheme[sys_id-1]->GetNumIter();};

  ///
  Integer GetSizeSystemVec(Integer sys_id) const
    {return sysmat[(sys_id-1)*blocksize+sys_id-1]->GetSize();};

  ///
  Double * GetSolutionVal(Integer sys_id) const 
    {return sol[sys_id-1]->GetPointer();};

  Integer * GetSolutionInd(Integer sys_id) const
    {return NULL;};

  ///
  void Print(Integer matrix_id, Integer matrix_row, Integer matrix_col) const 
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->Print();};

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
  Integer numsys, numgraph, numrhs, numbases, count, offset, numpart, blocksize, size, nne, dof, globsize;

private:

  ///
  void RScalarCreate(Integer sys_id, Integer ind, Integer size, Integer nne, 
		     Integer dir, Integer con, Integer * matrix_type);

  ///
  void CScalarCreate(Integer ind, Integer size, Integer nne, 
		     Integer dir, Integer con, Integer * matrix_type);

  ///
  void RBlockCreate(Integer sys_id, Integer ind, Integer size, Integer nne, 
		    Integer dir, Integer dof, Integer con, Integer * matrix_type);

  ///
  void CBlockCreate(Integer ind, Integer size, Integer nne, 
		    Integer dir, Integer dof, Integer con, Integer * matrix_type);

  ///
  void RFullCreate(Integer ind, Integer size, Integer nne, 
		   Integer dir, Integer dof, Integer con, Integer * matrix_type);

  ///
  void CFullCreate(Integer ind, Integer size, Integer nne, 
		   Integer dir, Integer dof, Integer con, Integer * matrix_type);

  ///
  Integer * glob2loc;
};

}

#endif // FILE_ALGEBRAICSYSTEM_PILES
