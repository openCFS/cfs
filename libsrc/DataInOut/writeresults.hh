#ifndef FILE_WRITERESULTS_2002
#define FILE_WRITERESULTS_2002

#include "grid.hh"

namespace CoupledField
{

template<class Dim>
class WriteResults
{
public:
   //!
   enum nameSol{fluid, temperature};

   /// initialization with grid
  virtual void Init(Grid<Dim> * aptgrid)=0;

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
  Grid<Dim> * ptgrid;

}; 

template<class Dim>
WriteResults<Dim>::~WriteResults()
{
#ifdef TRACE
  (*trace) << "entering WriteResults::~WriteResults" << std::endl;
#endif

 if (!output) delete output;

}

#ifdef __GNUC__
template class WriteResults<Point2D>;
template class WriteResults<Point3D>;
#endif

} // end of namespace
#endif
