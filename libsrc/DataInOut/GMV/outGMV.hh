#ifndef FILE_OUTGMV_2002
#define FILE_OUTGMV_2002

#include "grid.hh"

namespace CoupledField
{

//! This class provides an interface for writing files in gmv-format
template<class Dim>
class OutGMV
{
public:
  //! Constructor
  OutGMV(const Char * filename, Grid<Dim> * aptgrid);
  
  //! Deconstructor
  virtual ~OutGMV();
  
  //! write information about grid with level in file
  void Write(const Integer level);

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

};

#ifdef __GNUC__
template class OutGMV<Point3D>;
template class OutGMV<Point2D>;
#endif

} // end of namespace

#endif
