#ifndef FILE_ANSYSFILE_2002
#define FILE_ANSYSFILE_2002

#include "filetype.hh"

namespace CoupledField
{

  //! base class for reading initial data
  /*! 
    Class, that is derived from class FileType for reading input data, which is produced by Ansys interface. 
  */

class AnsysFile: virtual public FileType 
{

public:

  //!
  AnsysFile(const Char * const afilename);

  //!
  virtual ~AnsysFile();

  //!
  virtual void ReadMaxnumnodes(Integer & );

  //!
  virtual void ReadMaxnumelem(Integer & );

  //!
  virtual void ReadGeneralAnalChoice(Integer * dataGAnalCh,
				      enum nameGAn first ...)=0;

  //!
  virtual void ReadGeneralAnal(Integer * dataGAnal)=0;

  //! num - number of group element 
  virtual void ReadGeneralElemChoice(const Integer num, Integer * dataGElemCh,
				     enum nameGElem first ...)=0;

  //!
  virtual void ReadElemConnectionGH(const Integer maxelem, Integer * connect, const Integer maxnode, const Integer numelemgr, const Integer startposinarrayconn)=0;

  //! read maximum number of elements in group number num
  virtual void ReadMaxnumelemGroup(Integer & maxelem, const Integer numgr)=0;

  //!
  virtual void ReadCoordinate(Point3D * const NodesCoord,                                                     const Integer maxnumNodes);

  //!
  virtual void ReadCoordinate(Point2D * const NodesCoord,
			      const Integer maxnumNodes);

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

  virtual void ReadBoundRestr(std::list<NodeRestraint> & restr, Integer & numberRestr)=0;

  //! Read step data; special for TransientDriver
  virtual void ReadStepData(Integer & anumstep, Integer & aisavebegin, Integer & aisaveend, Integer & aisaveincr, Double & afirstdt)=0;

protected:
  //! 
  std::ifstream infile;
  
private:

  std::string::size_type pos_end;

  // dimension of problem
  Integer dim_;

  //!
  void ReadDim();

  // take position in section
  void takePosition(const std::string seekexp, std::string::size_type & pos);

  // transform string for level of boundary condition in number
  Integer TransformInDof(const Char * el);

  // read number of nodes for boundary condition
  void ReadMaxnumnodesbc(Integer & nbc);

};

}

#endif // FILE_FILETYPE
