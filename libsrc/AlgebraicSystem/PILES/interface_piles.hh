#ifndef FILE_INTERFACE_PILES_2002
#define FILE_INTERFACE_PILES_2002

#include <abstractAlgSys.hh>
#include <general.hh>
#include <multigrid.hh>

 
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
  virtual void SetSolverParameter(Integer nsys, Integer eps, Integer dampiter, Integer maxnumit,
                                  Integer solvertype, Integer precondtype)
  {  
    algsys->CreateParameter();
    algsys->SetAccuracy(eps,nsys);
    algsys->SetMaxNumIter(maxnumit,nsys);
    algsys->SetPrecond(precondtype,nsys);
    algsys->SetSolver(solvertype,nsys);
    algsys->SetDampIter(dampiter,nsys);
  }

  virtual void InitAlgSysGraph(Integer numnode)
  { algsys->CreateGraph(&numnode);}

  virtual void SetAlgSysGraph(Integer *pos, Integer elemsize, Integer nsys)
  { algsys->SetElementPos(pos,elemsize,nsys);}

  virtual void CreateAlgSysMatrices(Integer matrixtype, Integer *matrixsystype, Integer graphtype, 
                                     Integer numdofpernode, Integer numdirichlets, Integer numconstraints)
  {
    algsys->CreateMatrix(&matrixtype, matrixsystype, &graphtype, &numdofpernode, &numdirichlets,
                         &numconstraints);
  }


private:
  AlgebraicSystem * algsys;

} ;

} // end of namespace

#endif
