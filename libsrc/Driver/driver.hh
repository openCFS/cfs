#ifndef FILE_DRIVER_2001
#define FILE_DRIVER_2001

namespace CoupledField
{

  /// class where we implement Newmark method
  class Driver
  {
  public:
    ///
    Driver(FileType * const aptFileType, Integer anummesh, Material * ptMaterial);
    ///
    void SolveNewmarkMethod(OutResultUnverg * aptUnverg);
    ///
    virtual ~Driver();
  
  protected:
    ///
    Grid<Point2D> * ptgrid;
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
