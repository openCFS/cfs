#ifndef FILE_DRIVER_2001
#define FILE_DRIVER_2001

#include "filetype.hh"
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
    void SolveNewmarkMethod(OutResultUnverg<Dim> * aptUnverg);
    ///
    virtual ~Driver();
  
  protected:
    ///
    Grid<Dim> * ptgrid;
    ///
    AcousticPDE * ptAcPDE;
    ///
    FileType * ptFileType; 
    ///
    Material * ptMaterial;
    ///
    Integer numsteps; 
    Double dt0;
  
  private:
    Boolean SaveDer1;
    Boolean SaveDer2;
  };

}

#endif // FILE_DRIVER
