#ifndef FILE_WRITERESULTS_2002
#define FILE_WRITERESULTS_2002

#include "Domain/grid.hh"

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
   WriteResults(const Char * const filename,Boolean withHistory);

   //! initialization with grid
  virtual void Init(Grid * aptgrid)=0;

   //! deconstructor
   virtual ~WriteResults();

  //! write information about grid on level i in file
  /*!
    \param level in: level of grid
  */
  virtual void WriteGrid(const Integer level)=0;

  //! write information about the solution
  /*!
    \param sol vector with solution
    \param step step of calculation
    \param time time of calculation
    \param title name for the solution
  */
  virtual void WriteSolution(const Vector<Double> & sol, const Integer step, const Double time, const std::string title)=0;

  //! write cell data
  /*!
    \param data vector with data (ex. value of an error for the cell)
    \param step step of calculation
    \param time time of calculation
    \param title name for the data
  */
  virtual void WriteDataOnCell(const Vector<Double> & data, const Integer step, const Double time, const std::string title)=0;

   //! write cell vector-data
  /*!
    \param data pointer to vector with data
     \param step step of calculation
    \param time time of calculation
    \param title name for the data
  */
  virtual void WriteVecDataOnCell(const Vector<Double> * data, const Integer step, const Double time, const std::string title)
  { Error("Not implemented",__FILE__,__LINE__);}


  //! to open new file for printing results only for GMV
  /*!
    \number number number for output-file (ex. result.gmv001)
  */
  virtual void OpenFile(const Integer number)
  { Error("Not implemented",__FILE__,__LINE__);}

  //! check, is it the gmv-output file
  virtual Boolean IsGMV()=0;

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
  //! last time step, for which results were printed
  Vector<Double> lastsavetime;
  //! indicator: format of output: ascii or binary
  Boolean ascii_;

private:
  //! open history-files
  void InitHistoryFiles();

}; 

} // end of namespace
#endif
