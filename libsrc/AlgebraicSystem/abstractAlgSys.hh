#ifndef FILE_ABSTRACTAlgSys_2001
#define FILE_ABSTRACTAlgSys_2001

#include "tools.hh"
 
namespace CoupledField
{

class AbstractAlgebraicSys
{
public:
  //! Constructor 
  AbstractAlgebraicSys(){};

  //!
  virtual void InitAlgSys(Integer anumsys, Integer anumgraph)=0;

  //!
  virtual void SetSolverParameter(Integer nsys, Double eps, Double dampiter, Integer maxnumit,
                                  Integer solvertype, Integer precondtype)=0;

  //!
  virtual void InitAlgSysGraph(Integer numnode, Integer matrix_row, Integer matrix_col)=0;

  //!
  virtual void SetAlgSysGraph(Integer *pos, Integer elemsize, Integer matrix_row, Integer matrix_col)=0;

  //!
  virtual void CreateAlgSysMatrices(Integer maxtrix_row, Integer matrix_col, Integer *matrixsystype, Integer matrixtype, 
				    Integer graphtype, Integer numdofpernode, Integer numdirichlets, Integer numconstraints)=0;

  //!
  virtual void ResetAlgSys(Integer nsys, Integer nsys, Integer matrix_id)=0;

  //!
  virtual void ResetRHS(Integer matrix_row)=0;

  //!
  virtual void PutElemMatAlgSys(Double *elemmat, Integer *help, Integer numnodelem, Integer matrix_row,
				Integer matrix_col, Integer matrix_id)=0;
  //!
  virtual void ComputeSysMatrix(Integer matrix_row, Integer matrix_col, Double * matrix_factor)=0;

  //!
  virtual void SetBCDirichlet(Integer restrnum, Integer nodenr, Double restrval, Integer ndof, 
			      Integer matrix_row, Integer matrix_col, Integer matrix_id)=0;
  //!
  virtual void SetBCDirichletUpdate(Integer restrnum, Integer nodenr, Double restrval, Integer ndof, 
			      Integer matrix_row, Integer matrix_col, Integer matrix_id)=0;

  //!
  virtual void UpdateRHS(Integer matrix_row, Integer matrix_col, Integer matrix_id, Double * vec)=0;

  //!
  virtual void ComputePrecond(Integer job, Integer nsys)=0;

  //!
  virtual void SolveAlgSys(Integer nsys)=0;

  //!
  virtual Double * GetSolution(Integer nsys)=0;

protected:
  

} ;

} // end of namespace

#endif
