#ifndef FILE_WRITERESULTS_2002
#define FILE_WRITERESULTS_2002

#include "Domain/grid.hh"
#include "Domain/bcs.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/elemstoresol.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{

 //! base class for printing results
  /*! 
    Base class for printing results. Functions of the class are virtual in order to handle the different types of output files.
    At present {\tt gmv}-format, {\tt unverg}-format and {\tt database} are implemented.
  */ 
class WriteResults
{
public:
   //! constructor
   WriteResults(const Char * const filename, FileType * const aInFile=NULL);

  //! initialization with grid
  //! \param ptgrid pointer to grid object
  //! \param aptbcs pointer to BCs object
  virtual void Init(Grid * aptgrid, BCs * aptbcs)=0;

   //! deconstructor
   virtual ~WriteResults();

  //! write information about grid on level i in file
  /*!
    \param level in: level of grid
  */
  virtual void WriteGrid(const Integer level)=0;
 
  // *************************
  //     TRANSIENT SECTION
  // *************************
  
  //! write element solution vector (transient/static)
  /*!
    \param data vector with data (ex. value of an error for the cell)
    \param step step of calculation
    \param time time of calculation
  */
  virtual void WriteNodeSolutionTransient(const NodeStoreSol<Double>& data, 
					  const Integer step, 
					  const Double time) = 0;

 //! write element solution vector (transient/static)
  /*!
    \param data vector with data (ex. value of an error for the cell)
    \param step step of calculation
    \param time time of calculation
  */
  virtual void WriteElemSolutionTransient(const ElemStoreSol<Double>& data, 
					  const Integer step, 
					  const Double time) = 0;

  //! write node history vector (transient/static)
  /*!
    \param data vector with data (ex. value of an error for the cell)
    \param step step of calculation
    \param time time of calculation
  */
  virtual void WriteNodeHistoryTransient(const NodeStoreSol<Double>& data, 
				 const Integer step, 
				 const Double time);
  

  // *************************
  //     HARMONIC SECTION
  // *************************

  //! write element solution vector (harmonic)
  /*!
    \param data vector with data (ex. value of an error for the cell)
    \param step step of calculation
    \param frequency frequency of exciting function
    \param format format for writing complex solution (real-imag/amplitude-phase)
  */
  virtual void WriteNodeSolutionHarmonic(const NodeStoreSol<Complex>& data, 
					 const Integer step,
					 const Double frequency,
					 const ComplexFormat format) = 0;

 //! write element solution vector (harmonic)
  /*!
    \param data vector with data (ex. value of an error for the cell)
    \param step step of calculation
    \param frequency frequency of exciting function
    \param frequencyStep step of calculation
    \param format format for writing complex solution (real-imag/amplitude-phase)
  */
  virtual void WriteElemSolutionHarmonic(const ElemStoreSol<Complex>& data, 
					 const Integer step,
					 const Double frequency,
					 const ComplexFormat format) = 0;

  //! write nodal history vector (harmonic)
  /*!
    \param data vector with data (ex. value of an error for the cell)
    \param step step of calculation
    \param frequency frequency of exciting function
    \param frequencyStep step of calculation
    \param format format for writing complex solution (real-imag/amplitude-phase)
  */
  void WriteNodeHistoryHarmonic(const NodeStoreSol<Complex>& data, 
				const Integer step,
				const Double frequency,
				const ComplexFormat format);

  //! to open new file for printing results only for GMV
  /*!
    \number number number for output-file (ex. result.gmv001)
  */
  virtual void OpenFile(const Integer number)
  { Error("Not implemented",__FILE__,__LINE__);}

  //! check, is it the gmv-output file
  virtual Boolean IsGMV()=0;


  /// writes a matrix of the form ptCoordX ptCoorY ptCoordZ SolDof1 SolDof2 ...
  /*!
    \param sol solution
    \param nrDofs nr of degrees of freedom in solution
  */
  void WriteSolMatrix(Grid * ptgrid, const Integer level, const Vector<Double> sol, 
		      const std::string matFileName, const Integer nrDofs=1);


protected:
  //! Convertes enum SolutionType to string
  virtual std::string SolutionTypeToString(const SolutionType type) const
  {Error(" Not implemented here", __FILE__, __LINE__);}

  //! name of file for output results
  std::string namefile_;

  //! pointer to Grid
  Grid * ptgrid;

  //! pointer to BCs
  BCs * ptBCs_;

  // ************************************
  //  Section dealing with history files 
  // ************************************

  //! indicator: print history file or not
  Boolean NeedHistory_;

  //! vector with type of output values
  //! for which a history file is written
  StdVector<SolutionType> histQuantities_;
  
  //! history nodes per PDE
  StdVector<StdVector<Integer> > histNodesPerPDE_;
  

  //! pointer to ofstream with history information
  StdVector<StdVector<std::ofstream*> >  historyFiles_;

   //! indicator: format of output: ascii or binary
  Boolean ascii_;

  //! Ptr to input file, needed for reading the save nodes
  FileType * pt2Inputfile_;
  

protected:
  //! intialize the management of history files
  virtual void InitHistoryFiles();

}; 

} // end of namespace
#endif
