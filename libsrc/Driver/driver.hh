#ifndef FILE_DRIVER_2001
#define FILE_DRIVER_2001

#include "tools.hh"
#include "filetype.hh"
#include "conffile.hh"
#include "outGMV.hh"
#include "outUnverg.hh"
#include "pde.hh"

namespace CoupledField
{

  /// class where we implement Newmark method
  template<class Dim>
  class Driver
  {
  public:
    ///
    Driver(FileType * const aptFileType, Integer anummesh, Material * ptMaterial);
    ///
    void SolveNewmarkMethod(WriteResults<Dim> * aptUnverg);
    ///
    virtual ~Driver();
  
  protected:
    ///
    Grid<Dim> * ptgrid;
    ///
    AcousticPDE<Dim> * ptAcPDE;
    ///
    FileType * ptFileType; 
    ///
    Material * ptMaterial;
    ///
    Integer numsteps; 
    ///
    Double dt0;
  
  private:
    Boolean SaveDer1;
    Boolean SaveDer2;

  //! Output results(solution, first, second derivatives of solution on current step) in .unverg format

   void WriteResultsInFile(WriteResults<Dim> *, PDE * , const Integer, const Double );

  };

#ifdef __GNUC__
template class Driver<Point2D>;
template class Driver<Point3D>;
#endif

} // end of namespace
#endif // FILE_DRIVER
