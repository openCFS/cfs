#ifndef FILE_WRITERESULTS_2002
#define FILE_WRITERESULTS_2002

#include "grid.hh"

namespace CoupledField
{

class WriteResults
{
public:
   /// constructor
   WriteResults();

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

protected:
  ///
  std::ofstream * output;
  ///
  Grid * ptgrid;
  ///
  Boolean NeedHistory_;
  ///
  Integer history_node_;
  ///
  std::ofstream historyfile;
  ///
  void AddInHistory(const Double time, const Double val);

}; 

} // end of namespace
#endif
