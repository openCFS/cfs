#ifndef FILE_ABSTRACTAlgSys_2001
#define FILE_ABSTRACTAlgSys_2001
 
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
  virtual void SetSolverParameter(Integer nsys, Integer eps, Integer dampiter, Integer maxnumit,
                                  Integer solvertype, Integer precondtype)=0;

  //!
  virtual void InitAlgSysGraph(Integer numnode)=0;

  //!
  virtual void SetAlgSysGraph(Integer *pos, Integer elemsize, Integer nsys)=0;

  //!
  virtual void CreateAlgSysMatrices(Integer matrixtype, Integer *matrixsystype, Integer graphtype, 
                                     Integer numdofpernode, Integer numdirichlets, Integer numconstraints)=0;

protected:
  

} ;

} // end of namespace

#endif
