#ifndef FILE_WRITERESULTS_2002
#define FILE_WRITERESULTS_2002

#include "grid.hh"

namespace CoupledField
{

class WriteResults
{
public:
   /// constructor
   WriteResults(const Char * const filename,Boolean withHistory);

   //!
   enum nameSol{fluid, temperature};

   /// initialization with grid
  virtual void Init(Grid * aptgrid)=0;

   /// deconstructor
   virtual ~WriteResults();

  /// write information about grid on level i in file
  virtual void WriteGrid(const Integer level)=0;

  /// write information about the solution
  virtual void WriteSolution(const Vector<Double> & sol, const Integer step, const Double time, const std::string title)=0;

  //! write cell data
  virtual void WriteDataOnCell(const Vector<Double> & data, const Integer step, const Double time, const std::string title)=0;

   //! write cell data
  virtual void WriteVecDataOnCell(const Vector<Double> * data, const Integer step, const Double time, const std::string title)
  { Error("Not implemented",__FILE__,__LINE__);}


  //! open file - only for GMV
  virtual void OpenFile(const Integer number)
  { Error("Not implemented",__FILE__,__LINE__);}

  virtual Boolean IsGMV()=0;

protected:
  //! name of file for output results
  Char * namefile_;
  ///
  Grid * ptgrid;
  ///
  Boolean NeedHistory_;
  ///
  std::vector<Integer> nodeshist_;
  ///
  std::ofstream * historyfile;
  ///
  void AddInHistory(const Double time, const Double val, const Integer ifile);
  ///
  Vector<Double> lastsavetime;
  //!
  Boolean ascii_;

private:
  //!
  void InitHistoryFiles();

}; 

} // end of namespace
#endif
