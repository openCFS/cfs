#ifndef FILE_ABSTRACTAlgSys_2001
#define FILE_ABSTRACTAlgSys_2001

#include "tools.hh"
 
namespace CoupledField
{

//! Abstract class for algebraic system
class AbstractAlgebraicSys
{
public:
  //! Constructor 
  AbstractAlgebraicSys(){};

  //! Deconstructor
  virtual ~AbstractAlgebraicSys(){};

  //! Initialize the algebraic system (matrix block system)
  /*! 
    \param anumsys number of systems
    \param anumgraph number of needed (node, edge) graphs 
  */
  virtual void InitAlgSys(Integer anumsys, Integer anumgraph)=0;

  //! set the solver parameters
  /*!
    \param nsys number of the system within the block structure
    \param eps value for accuracy
    \param dampiter damping factor within the iterative solution
    \param maxnumit maximum number of iterations
    \param solvertype type of solver:  1.. RICHARDSON; 2..CG
    \param precondtype type of preconditioner: ID .. 1; MG..2
    \param numeqcoarse number of equations for coarse algebraic system (when using AMG)
  */
  virtual void SetSolverParameter(Integer nsys, Double eps, Double dampiter, Integer maxnumit,
                                  Integer solvertype, Integer precondtype, Integer numeqcoarse,
				  Double coarsealpha)=0;

  //! Initialize the graph of the nodal (edge) mesh
  /*!
    \param numelements number of nodes (edges) of mesh (specifies the dimension of the graph)
    \param matrix_row  row index in block structure
    \param matrix_col column index in block structure
  */
  virtual void InitAlgSysGraph(Integer numelements, Integer matrix_row, Integer matrix_col,
			       Integer matrix_gaphtype)=0;

  //! Set up the graph of the nodal (edge) mesh
  /*!
    \param pos global node (edge) numbers of one finite element
    \param elemsize number of nodes (edges) per element
    \param matrix_row  row index in block structure
    \param matrix_col column index in block structure
  */
  virtual void SetAlgSysGraph(Integer *pos, Integer elemsize, Integer fe_type, Integer matrix_row, 
			      Integer matrix_col)=0;

  //! Creates all matrices for the algebraic system
  /*!
    \param matrix_row  row index in block structure
    \param matrix_col column index in block structure
    \param matrixtype array of needed matrices:  matrixtype[0]: SYSTEM; matrixtype[1]: STIFFNESS;
           matrixtype[2]: DAMPING; matrixtype[3]: CONVECTION; matrixtype[4]:  MASS;
           set array value to 1, if matrix type is needed         
    \param matrix_class type of matrices:  1..RSPARSE; 2..CSPARSE; 3..RBLOCK; 4..CBLOCK;
           5..RFULL; 6..CFULL; 7..MIXED;
    \param graphtype dummy argument (currently not used, will be for specifying node, edge, etc. graph)!!!!
    \param numdofpernode number of degree of freedom per node (edge)
    \param numdirichlets number od restraints (inhomogeneous dirichlet BC)
    \param constraints number of constraints
  */
  virtual void CreateAlgSysMatrices(Integer maxtrix_row, Integer matrix_col, Integer *matrixsystype, 
				    Integer matrixtype, Integer graphtype, Integer numdofpernode, 
				    Integer numdirichlets, Integer numconstraints)=0;

  //!Set matrices, rhs-vector and solution vector of algebraic system with ID nsys to zero
  /*!
    \param nsys system id
    \param matrix_id type of matrix (mass, stif,..)
  */
  virtual void ResetAlgSys(Integer nsys, Integer nsys2, Integer matrix_id)=0;

  //! Set RHS-vector of system matrix_row to zero
  virtual void ResetRHS(Integer matrix_row)=0;

  //! Set matrix to zero
  /*!
    \param matrix_row  row index in block structure
    \param matrix_col column index in block structure
    \param matrix_id type of matrix: 1..System; 2..Stiffness; 3..Damping; 4..Convective; 5..Mass
  */
  virtual void ResetMatrix(Integer matrix_row, Integer matrix_col, Integer matrix_id)=0;

  //! Put elemtmatrix to global matrix
  //! Set matrix to zero
  /*!
    \param elemt matrix values
    \param global node (edge) numbers
    \param number of nodes (edges) of element
    \param matrix_row  row index in block structure
    \param matrix_col column index in block structure
    \param matrix_id type of matrix (mass,stiff,..)
  */
  virtual void PutElemMatAlgSys(Double *elemmat, Integer *pos, Integer numnodelem, Integer matrix_row,
				Integer matrix_col, Integer matrix_id)=0;

  //! Computes the system matrix
  /*! 
    \param graph_row  row index in matrix structure
    \param graph_col column index in matrix structure
    \param matrix_fac factors, with which the matrices (mass, stiffness,..) are multiplyed
                        before being added to the system matrix                
  */
  virtual void ComputeSysMatrix(Integer matrix_row, Integer matrix_col, Double * matrix_factor)=0;

  //! Set the Dirichlet boundary conditions (penalty formulation)
  /*!
    \param restrnum number of Dirichlet boundary condition
    \param nodenum node (edge) number)
    \param restrval value of Dirichlet boundary condition
    \param ndof degree of freedom component
    \param graph_row  row index in matrix structure
    \param graph_col column index in matrix structure
    \param matrix_id system ID (for effective mass/stiffness matrix)
  */
  virtual void SetBCDirichlet(Integer restrnum, Integer nodenr, Double restrval, Integer ndof, 
			      Integer matrix_row, Integer matrix_col, Integer matrix_id)=0;

  //! Update Dirichlet boundary conditions (e.g. in transient analysis)
  /*!
    \param restrnum number of Dirichlet boundary condition
    \param nodenum node (edge) number)
    \param restrval value of Dirichlet boundary condition
    \param ndof degree of freedom component
    \param graph_row  row index in matrix structure
    \param graph_col column index in matrix structure
    \param matrix_id system ID (for effective mass/stiffness matrix)
  */
  virtual void SetBCDirichletUpdate(Integer restrnum, Integer nodenr, Double restrval, Integer ndof, 
			      Integer matrix_row, Integer matrix_col, Integer matrix_id)=0;

  //! Update the right hand side (e.g. in transient analysis)
  /*!
    \param graph_row  row index in matrix structure
    \param graph_col column index in matrix structure
    \param matrix_id specifies type od matrix (mass, stiffness, ..)
    \param vec vector, which will be added to the right hand side
  */
  virtual void UpdateRHS(Integer matrix_row, Integer matrix_col, Integer matrix_id, Double * vec)=0;

  //! Compute the preconditioner
  /*!
    \param job: 1..for nonlinear stuff; 2..new computation; 3..just update Dirichlet values
    \param sys_id ID for system in block matrix
  */

  //! Add element vector to rhs-vector at position pos
  virtual void PutElemRHS(Double * elemrhs, Integer *pos, Integer size, Integer sys_id)=0;

  virtual void ComputePrecond(Integer job, Integer nsys)=0;

  //! Solve the algebraic system with ID nsys
  virtual void SolveAlgSys(Integer nsys)=0;

  //! Get the solution for algebraic system with ID nsys
  virtual Double * GetSolution(Integer nsys)=0;

  //! Compute the energy norm of defined matrix with vector u
  /*!
    \param graph_row  row index in matrix structure
    \param graph_col column index in matrix structure
    \param matrix_id specifies type od matrix (mass, stiffness, ..)
    \param vec vector, which will teh nergy norm will be computed
  */
  virtual Double CalcEnergyNorm(Integer matrix_row, Integer matrix_col, Integer matrix_id, Double * u)=0;

protected:
  

} ;

} // end of namespace

#endif
