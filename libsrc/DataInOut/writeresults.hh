#ifndef FILE_WRITERESULTS_2002
#define FILE_WRITERESULTS_2002

#include "grid.hh"

namespace CoupledField
{

template<class Dim>
class WriteResults
{
public:

   /// initialization with grid
  virtual void Init(Grid<Dim> * aptgrid)=0;

  /// write information about grid on level i in file
  virtual void WriteGrid(const Integer level)=0;

  /// write information about the solution
  virtual void WriteSolution(const Vector<Double> & sol, const Integer step, const Integer time)=0;

  /// write information about first derivatives of the solution
  virtual void WriteFirstDerSolution(const Vector<Double> & sol, const Integer step, const Integer time)=0;

  /// write information about second derivatives of the solution
   virtual void WriteSecondDerSolution(const Vector<Double> & sol, const Integer step, const Integer time)=0;
  
protected:
  ///
  std::ofstream * output;
  ///
  Grid<Dim> * ptgrid;
}; 

#ifdef __GNUC__
template class WriteResults<Point2D>;
template class WriteResults<Point3D>;
#endif

} // end of namespace
#endif
