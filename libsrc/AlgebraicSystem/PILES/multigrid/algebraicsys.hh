#ifndef FILE_ALGEBRAICSYSTEM_CLA
#define FILE_ALGEBRAICSYSTEM_CLA

namespace CoupledField
{

//! Class for Algebraic System
class AlgebraicSystem
{
public:
  //! Constructor
  /*! \param anumsys number of systems
      \param anumgraph number of needed graphs
  */
  AlgebraicSystem(Integer anumsys, Integer anumgraph);

  //! Destructor
  virtual ~AlgebraicSystem();

  virtual void CalcPrecond() = 0;

  //! calculate the preconditioner (here pure virtual)
  /*!
    \param newprecond job variable: 1..for nonlinear stuff; 2..new computation; 3..just update Dirichlet values
    \param sys_id ID for system in block matrix
  */
  virtual void CalcPrecond(Integer newprecond, Integer sys_id) = 0;

  //! Solve the algebraic system with ID sys_id (here just virtual)
  virtual void Solve() = 0;

  //!
  virtual void Solve(Double * f, Double * u) = 0;

  //! pure virtual (implemented in Class BlockSystem)
  virtual void CreateParameter() = 0;

  ////////////////////////////////////////////////////////////////////////////

  //! set the accuracy of the solver
  /*!
     \param aeps  numerical value of accuracy
     \param sys_id according system in the block structure
  */
  void SetAccuracy(Double aeps, Integer sys_id)
    {param[sys_id-1]->SetAccuracy(aeps);}

  //! set the maximum number of allowed iterations for the solver
  /*!
     \param amaxnumiter value for the maximum number of allowed iterations
     \param sys_id according system in the block structure
  */
  void SetMaxNumIter(Integer amaxnumiter, Integer sys_id)
    {param[sys_id-1]->SetMaxNumIter(amaxnumiter);}
  
  //! set the preconditioner
  /*!
     \param aprecond type of preconditioner
     \param sys_id according system in the block structure
  */
  void SetPrecond(Integer aprecond, Integer sys_id)
    {param[sys_id-1]->SetPrecond(aprecond);}

  //! set the type of solver 
  /*!
     \param asolver type of solver
     \param sys_id according system in the block structure
  */
  void SetSolver(Integer asolver, Integer sys_id)
    {param[sys_id-1]->SetSolver(asolver);}
  
  //! set the maximum number of allowed iterations for the solver
  /*!
     \param adampiter value for damping the iterative solver
     \param sys_id according system in the block structure
  */
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

  //! Put element matrix to global (stiffness, mass, ..) matrix
  /*!
    \param elemmat array of entries of element matrix
    \param connect array of global node numbers
    \param elemsize number of nodes (edges) per element
    \param graph_row  row index in matrix structure
    \param graph_col column index in matrix structure
    \param matrix_id specifies type od matrix (mass, stiffness, ..)
  */
  void SetElementMatrix(Double * elemmat, Integer * connect, Integer elemsize, 
			Integer matrix_row, Integer matrix_col,Integer matrix_id)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->Assemble(elemmat, connect, elemsize);}

  ///
  void SetAuxElementMatrix(Double * elemmat, Integer * connect, Integer elemsize ,Integer sys_id)
    {auxmat[sys_id-1]->Assemble(elemmat, connect, elemsize);}

  ///
  void SetElementRHS(Double * elemrhs, Integer * connect, Integer elemsize, Integer sys_id)
    {rhs[sys_id-1]->Assemble(elemrhs, connect, elemsize);};

  //! initialize right hand side of system with ID sys_id to zero
  void InitRHS(Integer sys_id)
    {rhs[sys_id-1]->Init();};

  //! initialize solution vector  of system with ID sys_id to zero
  void InitSol(Integer sys_id)
    {sol[sys_id-1]->Init();};

  //! Get memory for matrices of defined system and setup according matrix graph
  /*! \param graph_row  row index in matrix structure
      \param graph_col column index in matrix structure
      \param matrix_type array of needed matrices (system, stiffness, mass, ..)
      \param matrix_class type of matrices (real, complex, sparse,full,..)
      \param graph dummy argument (currently not used, will be for specifying node, edge, etc. graph)!!!!
      \param adof number of degree of freedom per node (edge)
      \param dir number od restraints (inhomogeneous dirichlet BC)
      \param con number of constraints
  */
  void CreateMatrix(Integer matrix_row, Integer matrix_col, Integer * matrix_type, 
		    Integer matrix_class, Integer graph, Integer adof, Integer dir, Integer con);

  //! initialize the solver for system with ID sys_id
  void CreateSolver(Integer sys_id);

  //! initialize preconditioner
  /*!
    \param matrix_class type of matrix (real, complex, sparse,full,..)
    \param sys_id ID of system
  */
  void CreatePrecond(Integer matrix_class, Integer sys_id);

  ///
  void SetMatrixFactor(Double val, Integer matrix_row, Integer matrix_col, Integer matrix_id)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->SetMatrixFactor(val);};

  ///
  void MultMat(Integer matrix_row, Integer matrix_col, Integer matrix_id, Double * u, Double * v)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->Mult(u,v);};

  double CalcEnergyNorm(Integer matrix_row, Integer matrix_col, Integer matrix_id, Double * u)
  {
    return sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->CalcEnergyNorm(u);
  };

  ///
  void AddGMGLevel();

  //! Initialize matrix (set to zero)
  /*!  
    \param graph_row  row index in matrix structure
    \param graph_col column index in matrix structure
    \param matrix_id specifies type od matrix (mass, stiffness, ..)
  */
  void InitMatrix(Integer matrix_row, Integer matrix_col, Integer matrix_id)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->InitMatrix();};

  ///
  void SetMatrixPointer(Integer * astart, Integer * apos, Double * aval, 
			Integer matrix_row, Integer matrix_col, Integer matrix_id)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->SetPointer(astart, apos, aval);};

  /////////////////////////////////////////////////////////////////////////////////////////
  
  //! update right hand side (e.g. in transient analysis)
  /*!
    \param graph_row  row index in matrix structure
    \param graph_col column index in matrix structure
    \param matrix_id specifies type od matrix (mass, stiffness, ..)
    \param fup vector, which will be added to the right hand side
  */
  void UpdateRHS(Integer matrix_row, Integer matrix_col, Integer matrix_id, Double * fup)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->MultAdd(fup, *rhs[matrix_row-1]);};

  //! set Dirichlet boundary conditions (using penalty formulation)
  /*!
    \param i number of Dirichlet boundary condition
    \param num node (edge) number)
    \param val value of Dirichlet boundary condition
    \param comp degree of freedom component
    \param graph_row  row index in matrix structure
    \param graph_col column index in matrix structure
    \param matrix_id system ID (for effective mass/stiffness matrix)
  */
  void SetDirichlet(Integer i, Integer num, Double val, Integer comp, 
		    Integer matrix_row, Integer matrix_col, Integer matrix_id)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->SetDirichlet(i, num, val, comp);};

  //! update Dirichlet boundary conditions (e.g. in transient analysis)
  /*!
    \param num number of Dirichlet boundary condition
    \param val value of Dirichlet boundary condition
    \param graph_row  row index in matrix structure
    \param graph_col column index in matrix structure
    \param matrix_id system ID (for effective mass/stiffness matrix)
  */
  void UpdateDirichlet(Integer num, Double val, Integer matrix_row, Integer matrix_col, Integer matrix_id)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->UpdateDirichlet(num, val);};

  ///
  void SetConstraint(Integer * cm, Double * cf, Integer matrix_row, Integer matrix_col, Integer matrix_id)
    {sysmat[(matrix_id-1)*offset+(matrix_row-1)*blocksize+matrix_col-1]->SetConstraint(cm, cf);};
  
  /////////////////////////////////////////////////////////////////////////////////////////

  ///
  void CreateGlobal2Local(Integer aglobsize);

  //! Initialize matrix graph for each system in the matrix structure
  /*! \param size number of elements (nodes, edges)
      \param graph_row  row index in matrix structure
      \param graph_col column index in matrix structure
  */
  void CreateGraph(Integer size, Integer graph_row, Integer graph_col)
    {graph[(graph_row-1)*blocksize+graph_col-1] = new BaseGraph(size);};

  //! Set up neighbor relations of nodes (edges)
  /*! \param connect contains the global node (edge) numbers
      \param elemsize number of nodes (edges) in array connect
      \param graph_row  row index in matrix structure
      \param graph_col column index in matrix structure
  */
  void SetElementPos(Integer * connect, Integer elemsize, Integer matrix_row, Integer matrix_col)
    {graph[(matrix_row-1)*blocksize+matrix_col-1]->SetElementPos(connect, elemsize);};

  /////////////////////////////////////////////////////////////////////////////////////////

  //! compute the system matrix
  /*! \param graph_row  row index in matrix structure
      \param graph_col column index in matrix structure
      \param matrix_fac factors, with which the matrices (mass, stiffness,..) are multiplied
                        before being added to the system matrix                
  */
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

  //! 
  BaseMatrix * sysmat[80], * spemat[4], * auxmat[16];

  //!
  BaseGraph * graph[4];

  //!
  BaseVector * rhs[4], * sol[4];

  //!
  BaseSolver * scheme[4];

  //!
  BasePrecond * premat[4];

  //! Stores the solver parameter for each system in the block structure
  BaseParameter * param[4];

  //!
  Integer numsys, numgraph, numrhs, numbases, count, offset, numpart, blocksize, size, nne, dof, globsize;

private:

  //! get momory for scalar matrices of system with ID sys_id
  //! if system matrix, then get memory for rhs- and solution-vector
  /*!
     \param sys_id system number
     \param ind index in block matrix structure
     \param size number of nodes (edges)
     \param nne number of nonzero entries in matrix
     \param dir number od restraints (inhomogeneous dirichlet BC)
     \param con number of constraints
     \param matrix_type array of needed matrices (system, stiffness, mass, ..)
  */
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

#endif // FILE_ALGEBRAICSYSTEM_CLA
