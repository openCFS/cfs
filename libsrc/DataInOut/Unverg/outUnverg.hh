#ifndef FILE_OUTRESULTUNVERG_2001
#define FILE_OUTRESULTUNVERG_2001

#include "grid.hh"
#include "writeresults.hh"

namespace CoupledField
{

//! Class for writing information about grid and results in unverg-format file

template<class Dim>
class WriteResultsUnverg: virtual public WriteResults<Dim>
{

public:
  /// constructor with name of a file for results
  WriteResultsUnverg(const Char * const filename); 

  /// deconstructor
  virtual ~WriteResultsUnverg();
  
   /// initialization with grid
  virtual void Init(Grid<Dim> * aptgrid);

  /// write information about grid on level i in file
  virtual void WriteGrid(const Integer level);

  /// write information about the solution
  virtual void WriteSolution(const Vector<Double> & sol, const Integer step, const Double time);

  /// write information about first derivatives of the solution
  virtual void WriteFirstDerSolution(const Vector<Double> & sol, const Integer step, const Double time);

  /// write information about second derivatives of the solution
   virtual void WriteSecondDerSolution(const Vector<Double> & sol, const Integer step, const Double time);

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
  ///
  Dim * ptCoordinate;
 
};

#ifdef __GNUC__
template class WriteResultsUnverg<Point2D>;
template class WriteResultsUnverg<Point3D>;
#endif

} // end of namespace
 
#endif
