#ifndef FILE_WRITERESULTS_2002
#define FILE_WRITERESULTS_2002

#include "grid.hh"

namespace CoupledField
{

class WriteResults
{
public:
   /// constructor
   WriteResults(const Char * const filename);

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

}; 

} // end of namespace
#endif
