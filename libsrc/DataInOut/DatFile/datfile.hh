#ifndef FILE_DATFILE
#define FILE_DATFILE

#include "filetype.hh"

namespace CoupledField
{

  //! Derived class from FileType for reading initial data from dat file
  /*!
    This class is derived from class FileType. It reads necessary data from {\it dat} - file, which is produced by CAPAPRE.
  */

class DatFile : virtual public FileType
{
public:
  //! Constructor
  DatFile(const Char * const afilename);

  //! Deconstructor
  ~DatFile();

protected:
  //! Auxialary Variables  
  std::ifstream infile;

  //! usefull enums
  //! enum for reading data about general analisys
//  enum nameGAn{soltype, analtype, numnode, numgroup, restart, inactdofs,
//               circuit, deactDf, masstype, damptype, nooptimiz, endGAnal}; //11

  //! enum for redaing general info about element
//  enum nameGElem {numelem, ielemtyp, isubtype, ielemsave, maxnode,
//                  nonlinear, form1, form2, endGElem}; // 8
public:

  //! Read maximum number of nodes
  virtual void ReadMaxnumnodes(Integer & maxnumnodes);

  //! Read maximum number of elements
  virtual void ReadMaxnumelem(Integer & maxnumelem);

  //!  Read a title of simulation run 
  Boolean ReadTitle(std::string & title);

  //! Print a title of simulation run
  void  PrintTitle (std::ostream * outfile);

  //! Read analysis data from General section
  void ReadGeneralAnal(Integer * dataGAnal);

  //!  Read analysis data from General section by choice
  void ReadGeneralAnalChoice(Integer * dataGAnalCh,
			     enum nameGAn first ...);

  //! Read output options, concerned saving first, second derivatives of vp
  void ReadOutputOptions(Boolean & Der1, Boolean & Der2);   

  //! Read output specification data 
  void ReadOutputData(Integer * dataOut);

  //!
  void ReadOutputDataArr(Integer * dataOutArr, Integer sizeArr,
			 std::string & seekStr);

  //! Read Nodal Coordinates (overloading for 3D and 2D) 
  void ReadCoordinate(Point3D * const InitNodalCo, const Integer maxnumNodes);
  void ReadCoordinate(Point2D * const InitNodalCo, const Integer maxnumNodes);


  //! Read General Info about element 
  void ReadGeneralElem(Integer * dataGElem, const Integer num);

  //! Read General Info about element by choice
  void ReadGeneralElemChoice(const Integer num, Integer * dataGElemCh,
			     enum nameGElem first ...);

  //! Read Info about connection of element
  void ReadElemConnect(Integer elemnum, Integer * connect,
		       Integer maxnode);

  //! Read connection of element in grid hierarchy
  void  ReadElemConnectionGH(const Integer maxelem, Integer * connect,
			     const Integer maxnode, const Integer numelemgr, const Integer startposinarrayconnection);

  //! read maximum number of elements in group number num
  void ReadMaxnumelemGroup(Integer & maxelem, const Integer numgr);

  //! Auxialary functions  
  void ReadElemRecord3d(Integer elemnum, Integer * connect,Integer maxnode,
                        const Integer numelemgr);

  //! Read info about connection elements of group number: numelemgr
  void ReadElem(const Integer maxelem, Integer ** connect,const Integer maxnode,
		const Integer numelemgr);
  //!
  void ReadElemMatCalc(Integer & matnumber, std::string & calc_expr);

  //!
  void ReadThickness(Double & thickness);

  //! Read a number of time functions 
  void ReadNumTimeFunc(Integer & maxnumTimeFunc);
  
  //! Read Max Value for a number of time functions
  void ReadInfoTimeFunc(Integer * maxvalTimeFunc,                                                 const Integer maxnumTimeFunc);

  //! Read time functions with number numTimeFunc
  void ReadTimeFunc(Double * const timeTimeFunc,
                    Double * const valIimeFunc,
                    const Integer numTimeFunc);
  //! Read integration parameters for Newmark method
  void ReadIntegrationParam(Double & alpha, Double & beta, Double & gamma);

  //! Read gamma_parabolic
  void ReadParabolicParam(Double & gamma_p);

  //! Read first time step
  void ReadNumStepsAndTimeSteps(Integer & numsteps, Double & dt);
  
  //! Read analisys specification data 
  void preReadTransAnal(Integer & soltype, Integer & statickey);

  //! Read analisys specification data
  void ReadTransAnalDir(Integer * dataTAnalDir, Double & bemtolD);

  //! Read analisys specification data
  void ReadTransAnalIter(Integer * dataTAnalIt, Double & bemtolIt,
			 Double & convtolIt);
 
  //! Read analisys specification data
  void ReadTransAnalNonl(Integer * dataTAnalNonl,
			 Double * dataTAnalNonlTol);

  //! Read analisys specification data 
  void ReadTransAnalTran(Integer * dataTAnalTran,
			 Double * dataTAnalTranReal);

  //! Read analisys specification data
  void ReadTransAnalStat(Integer & numstepStat,
			 Double & deltatStat );

  //! Read analisys specification data
  void ReadEigenvalAnal(Integer * dataEAnal,
			Double * dataEAnalReal );

  //! Read analisys specification data
  void ReadHarmAnal(Integer * dataHarm,
		    Double * dataHarmReal );

  //! Read general info about boundary conditions 
  void ReadGeneralBound(Integer * dataGBound, Double * dataGBoundD);

  //! Read general info about boundary conditions by choice
  void ReadGeneralBoundChoice(Integer * dataGBoundCh,
			      enum nameBound first ...);

  //! Read number of nodes at which there is Dirichlet boundary condition
  void ReadNumNodesforDirichletBC(Integer &);

  //! Read array of nodes at which there is Dirichlet boundary condition
  void ReadDirichletBC(Integer *);

  //! Read array of nodes at which there is Dirichlet boundary condition in vec
  void ReadDirichletBC(Vector<Integer> & nodes);

  //! Read specific info about boundary condition Dof
  void ReadBoundDof(Integer ** dataBDof, Integer numberdofs,Integer maxrecord);

  //! Read specific info about boundary condition Restrain 
  void ReadBoundRestr(Integer ** dataBRestr, Integer numberRestr,
		      Double * factorRestr); 

  //! Read specific info about boundary condition in list
  void ReadBoundRestr(std::list<NodeRestraint> & restr, Integer & numberRestr);
 
  //! Auxialary function for reading boundary condition
  void ReadMaxRecord(Integer & maxrecord);

  //! Read specific info about boundary condition Constraint
  void ReadBoundConstr(Integer ** dataBDof, Integer numberdofs,
		       Double * factor);

  //! Read specific info about boundary condition Load
  void ReadBoundLoads(Integer ** dataBLoads, Integer number,
		      Double * factor); 

  //! Read step data; special function for TransientDriver
  void ReadStepData(Integer & anumstep, Integer & aisavebegin, Integer & aisaveend, Integer & aisaveincr, Double & afirstdt);

  //! Enum for output specification data 
  enum nameOut{ formsave, nsavenode, savemax, savedofs, savederiv, 
		savederiv2}; // 6

  //! Enums for boundary specification data 
  enum nameTransDir{impexpDir, ncorrectDir, bemcheckDir}; // 3

  //!
  enum nameTransIt{itersoltype, ncorrectIt, bemcheckIt, iprecon,
		   iterlog, vecgmres}; //6

  //!
  enum nameTransNonl{nlinalg, nlcrit, nform, nlmaxit, updgeo}; //5

  //!
  enum nameTransTran{numstep, isavebegin, isaveend, isaveinc}; //4

  //!
  enum nameTransTranReal{deltatTr, alpha, beta, gamma_h, gamma_p}; //5

  //!
  enum nameEigen{method,neigval,maxiter}; //3

  //!
  enum nameEigenReal{flower,fupper,fshift,convtol}; //4

  //!
  enum nameHarm{cwnumf, cwspace,cwout}; //3

  //!
  enum nameHarmReal{cwlower,cwupper};

  //! Enum for part of boundary conditions 
  enum nameBoundD{regler_offset, regler_tol}; // 2

  //! Coding of the nodal degrees of freedom 

  enum nameDf{non, dx, dy, dz, ep, vp, chg, ax, ay, az, esp, anx, any, anz,
              tp, tnp, enp, vnp, mue, sig}; // 20

  //! Axialary function : transform element of enum in std::string 
  std::string PrintnameDf(Integer el);

  //!  Axialary function : transform element of enum in std::string
  std::string PrintnameDf(enum nameDf el);

  //!  Axialary function : transform std::string in element of enum 
  enum nameDf TransformInNameDf(const Char * el);

  //!  Axualary function : transform std::string in Char *  
  void Peel(std::string & buf, Char * init);

  //!  Axualary function : return position in file where seekexp is found first time
  void TakePos(const std::string seekexp,  std::string::size_type & pos);
  void TakePos(const std::string seekexp,  std::string::size_type & pos, std::string & buf_prev);

  //! Return TRUE if there is seekexp in datfile
  Boolean IsThere(const std::string seekexp) ;

  //! Count elements in std::string 
  Integer CountElemInString(std::string::size_type pos);

  //! Count number of character a in std::string buf 
  Integer CountCharInStr(std::string & buf, const Char a);

  //!  Axialary function: For coding number 
  Integer FindDigits(const Integer x);

  //!  Axialary function: For coding number
  Integer Step(Integer n, Integer a);

  //!  Axialary function: For coding number
  Integer Cryptogr(const Integer x);

  //!  Axialary function: For coding number
  Integer Encode(const Integer x);

  std::string::size_type pos_end; 
 
};
}

#endif // FILE_DATFILE
