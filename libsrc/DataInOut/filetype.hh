#ifndef FILE_FILETYPE_2001
#define FILE_FILETYPE_2001

#include "vector.hh"
#include <list>

namespace CoupledField
{

  //! base class for reading initial data
  /*! 
    Base class for reading input data. Its function is virtual due to handle the different types of input files and to hide their features from developer of code.
  */

struct NodeRestraint;

class FileType
{

public:
  //! enum for reading boundary condition
  enum nameBound{ numdofs, numconstr, numrestr, numloads, resistors,
		  numspring, bembdry, numflux, numrad, numpress, ncurrdens,
		  regl_itmx, reglr_use, endBound}; // 13   

  //! enum for reading data about general analisys
  enum nameGAn{soltype, analtype, numnode, numgroup, restart, inactdofs,
	       circuit, deactDf, masstype, damptype, nooptimiz, endGAnal}; //11

  //! enum for redaing general info about element
  enum nameGElem {numelem, ielemtyp, isubtype, ielemsave, maxnode,
		  nonlinear, form1, form2, endGElem}; // 8 

  //!
  FileType(const Char * const afilename);

  //!
  virtual ~FileType();

  //!
  virtual  void ReadGeneralAnalChoice(Integer * dataGAnalCh,
				      enum nameGAn first ...)=0;

  //! num - number of group element 
  virtual void ReadGeneralElemChoice(const Integer num, Integer * dataGElemCh,
				     enum nameGElem first ...)=0;

  //!
  virtual void ReadElemConnectionGH(const Integer maxelem, Integer * connect, const Integer maxnode, const Integer numelemgr)=0;

  //!
  virtual void ReadCoordinate(Point3D * const InitNodalCo,                                                     const Integer maxnumNodes)=0;

  //!
  virtual void ReadCoordinate(Point2D * const InitNodalCo,
			      const Integer maxnumNodes)=0;

  //!
  virtual void ReadIntegrationParam(Double &, Double &, Double &)=0;

  //!
  virtual void ReadParabolicParam(Double & gamma_p)=0;
  //!
  virtual void ReadNumStepsAndTimeSteps(Integer &, Double &)=0;

  //!
  virtual  void ReadTimeFunc(Double * const timeTimeFunc,
			     Double * const valIimeFunc,
			     const Integer numTimeFunc)=0;

  //!
  virtual  void ReadNumTimeFunc(Integer & maxnumTimeFunc)=0;    

  //!
  virtual  void ReadInfoTimeFunc(Integer * maxvalTimeFunc,
				 const Integer maxnumTimeFunc)=0;

  //! Read number of nodes at which there is Dirichlet BC
  virtual void ReadNumNodesforDirichletBC(Integer &)=0;

  //! Read array of nodes at which there is Dirichlet BC
  virtual void ReadDirichletBC(Integer *)=0;

  //! Read array of nodes at which there is Dirichlet BC in vec
  virtual void ReadDirichletBC(Vector<Integer> &)=0;

  //! Read output options, concerned saving first, second derivatives of vp
  virtual void ReadOutputOptions(Boolean & Der1, Boolean & Der2)=0;

 //! Read analisys specification data
  virtual void preReadTransAnal(Integer & soltype, Integer & statickey)=0;

  //! Read analisys specification data
  virtual void ReadTransAnalTran(Integer * dataTAnalTran,
                         Double * dataTAnalTranReal)=0;

  //! Read specific info about boundary condition Restrain
  virtual void ReadBoundRestr(Integer ** dataBRestr, Integer numberRestr,
                      Double * factorRestr)=0;

  virtual void ReadBoundRestr(std::list<NodeRestraint> & restr, const Integer numberRestr)=0;

  //! Read step data; special for TransientDriver
  virtual void ReadStepData(Integer & anumstep, Integer & aisavebegin, Integer & aisaveend, Integer & aisaveincr, Double & afirstdt)=0;

protected:

  //! name of input file
  Char * filename;
};

}


#endif // FILE_FILETYPE
