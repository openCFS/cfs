#ifndef FILE_OUTRESULTUNVERG_2001
#define FILE_OUTRESULTUNVERG_2001

#include "grid.hh"

namespace CoupledField
{
///
template <class Dim>
class OutResultUnverg
{
public:
  ///
  OutResultUnverg(const Char * const filename); 
  ///
  void Create(Grid<Dim> * ptgrid, const Integer level);
  ///
  void Dataset666(Grid<Dim> * ptgrid, const Integer level);
  ///
  void Dataset781(Grid<Dim> * ptgrid, const Integer level);   
  ///
  void Dataset780(Grid<Dim> * ptgrid, const Integer level);
  ///
  void Dataset55(const std::string & title, const Vector<Double> & x, const Integer step, const Double time);
  void Dataset56(); 
  ///
  ~OutResultUnverg();
 
private:
  ///
  std::ofstream * output;
  ///
  Dim * ptCoordinate;
 
};

#ifdef __GNUC__
template class OutResultUnverg<Point2D>;
template class OutResultUnverg<Point3D>;
#endif

} // end of namespace
 
#endif
