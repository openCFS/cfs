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
  WriteResultsUnverg(const Char * const filename,Boolean withHistory=FALSE); 

  /// deconstructor
  virtual ~WriteResultsUnverg();
  
   /// initialization with grid
  virtual void Init(Grid * aptgrid);

  /// write information about grid on level i in file
  virtual void WriteGrid(const Integer level);

  /// write information about the solution
  virtual void WriteSolution(const Vector<Double> & sol, const Integer step, const Double time, const std::string title);

   //! write cell data
  virtual void WriteDataOnCell(const Vector<Double> & data, const Integer step, const Double time, const std::string title);
  
  //
 virtual Boolean IsGMV() { return FALSE;}

private:
  ///
  std::ofstream * output;

  ///
  void Dataset666(const Integer level);

  ///
  void Dataset781(const Integer level);   

  ///
  void Dataset780(const Integer level);

  //! for printing nodal results of simulation
  void Dataset55(const std::string & title, const Vector<Double> & x, const Integer step, const Double time);

  //! for printing cell results of simulation
  void Dataset56(const std::string & title, const Vector<Double> & x, const Integer step, const Double time);

};

} // end of namespace
 
#endif
