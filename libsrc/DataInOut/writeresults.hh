#ifndef FILE_WRITERESULTS_2002
#define FILE_WRITERESULTS_2002

#include "Domain/grid.hh"
#include <Utils/array.hh>

namespace CoupledField
{

 //! base class for printing results
  /*! 
    Base class for printing results. Functions of the class are virtual in order to handle the different types of output files.
    At present {\tt gmv}-format and {\tt unverg}-format are implemented.
  */ 
class WriteResults
{
public:
   //! constructor
   WriteResults(const Char * const filename,Boolean withHistory, FileType * const aInFile=NULL);

   //! initialization with grid
  virtual void Init(Grid * aptgrid)=0;

   //! deconstructor
   virtual ~WriteResults();

  //! write information about grid on level i in file
  /*!
    \param level in: level of grid
  */
  virtual void WriteGrid(const Integer level)=0;

  //! write nodal soltion vector
  /*!
    \param sol array with solution
    \param step step of calculation
    \param time time of calculation
    \param title name for the solution
    \param nrDofs dimension of solution
  */
  virtual void WriteNodeSolution(const Array<Double>& sol, const Integer step, const Double time, const std::string title)=0;

  //! write element solution vector
  /*!
    \param data array with data (ex. value of an error for the cell)
    \param step step of calculation
    \param time time of calculation
    \param title name for the data    
    \param nrDofs dimension of solution
  */
  virtual void WriteElemSolution(const Array<Double> & data, const Integer step, const Double time, const std::string title)=0;

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


  /// read list of history nodes by name
  void ReadSaveNodes();

 
  

protected:
  //! name of file for output results
  Char * namefile_;

  //! pointer to Grid
  Grid * ptgrid;

  //! indicator: print history file or not
  Boolean NeedHistory_;

  //! vector with nodes history
  std::vector<Integer> nodeshist_;

  //! pointer to ofstream with history information
  std::ofstream * historyfile;

  //! add new data in history file
  void AddInHistory(const Double time, const Double val, const Integer ifile);

  //! add new vector valued data in history file
  void AddVecInHistory(const Double time, const std::vector<Double> val,const Integer ifile);

  //! last time step, for which results were printed
  Vector<Double> lastsavetime;

  //! indicator: format of output: ascii or binary
  Boolean ascii_;

  //! Ptr to input file, needed for reading the save nodes
  FileType * pt2Inputfile_;
  

private:
  //! open history-files
  void InitHistoryFiles();

}; 

} // end of namespace
#endif
