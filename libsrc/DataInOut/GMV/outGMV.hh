#ifndef FILE_OUTGMV_2002
#define FILE_OUTGMV_2002

#include "grid.hh"
#include "writeresults.hh"

namespace CoupledField
{

//! This class provides an interface for writing files in gmv-format
template<class Dim>
class WriteResultsGMV: virtual public WriteResults<Dim>
{
public:

  //! Constructor
  WriteResultsGMV(const Char * filename);
  
  //! Deconstructor
  virtual ~WriteResultsGMV();
  
  //! initialization with grid
  virtual void Init(Grid<Dim> * aptgrid);

  //! write information about grid with level in file
  virtual void WriteGrid(const Integer level);

  /// write information about the solution
  virtual void WriteSolution(const Vector<Double> & sol, const Integer step, const Double time, const std::string title);

private:
  
  //!
  Grid<Dim> * ptgrid;

  //!
  std::ostream * output;

  //! write header of gmv-file: only ascii is implemented
  void WriteHeader();

  //! write number of nodes and coordinates of them
  void WriteNodes(const Integer level);

  //! write cell description 
  void WriteCells(const Integer level); 

  //! write variable information
  void WriteVariable(const Vector<Double> var, const std::string name, const Integer type);

  //! 
  Integer variablemode_;  

};

#ifdef __GNUC__
template class WriteResultsGMV<Point3D>;
template class WriteResultsGMV<Point2D>;
#endif

} // end of namespace

#endif
