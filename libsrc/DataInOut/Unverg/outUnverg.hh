#ifndef FILE_OUTRESULTUNVERG_2001
#define FILE_OUTRESULTUNVERG_2001

namespace CoupledField
{

template<class T>
class Grid;
///
class OutResultUnverg
{
public:
  ///
  OutResultUnverg(const Char * const filename); 
  ///
  void Create(Grid<Point2D> * ptgrid, const Integer level);
  ///
  void Dataset666(Grid<Point2D> * ptgrid, const Integer level);
  ///
  void Dataset781(Grid<Point2D> * ptgrid, const Integer level);   
  ///
  void Dataset780(Grid<Point2D> * ptgrid, const Integer level);
  ///
  void Dataset55(const string & title, const Vector<Double> & x, const Integer step, const Double time);
  void Dataset56(); 
  ///
  ~OutResultUnverg();
 
private:
  ///
  ofstream * output;
  ///
  Point2D * ptCoordinate;
  ///
//  string tab; 
 
};

//template class outResultUnverg<Point2D>;
//template class outResultUnverg<Point3D>;

} // end of namespace
 
#endif
