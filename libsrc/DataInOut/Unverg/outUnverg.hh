#ifndef FILE_OUTRESULTUNVERG_2001
#define FILE_OUTRESULTUNVERG_2001

#include "grid.hh"
#include "writeresults.hh"

namespace CoupledField
{

//! Class for writing information about grid and results in unverg-format file

class WriteResultsUnverg: virtual public WriteResults
{

public:
  /// constructor with name of a file for results
  WriteResultsUnverg(const Char * const filename); 

  /// deconstructor
  virtual ~WriteResultsUnverg();
  
   /// initialization with grid
  virtual void Init(Grid * aptgrid);

  /// write information about grid on level i in file
  virtual void WriteGrid(const Integer level);

  /// write information about the solution
  virtual void WriteSolution(const Vector<Double> & sol, const Integer step, const Double time, const std::string title);

private:

  ///
  void Dataset666(const Integer level);

  ///
  void Dataset781(const Integer level);   

  ///
  void Dataset780(const Integer level);

  ///
  void Dataset55(const std::string & title, const Vector<Double> & x, const Integer step, const Double time);

  ///
  void Dataset56(); 
};

} // end of namespace
 
#endif
