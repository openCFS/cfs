#ifndef FILE_DOMAIN_2001
#define FILE_DOMAIN_2001

namespace CoupledField
{
class FileType;

/// contain information about calculation domain and according to different meshes create different grids

template<class Dim>
class Domain
{
public:
  ///
  Domain(Integer anumsubdomain, FileType * const aptFileType);

  ///
  virtual ~Domain();

  ///
  void SetSubdomains();

  ///
  void PrintDomain();

protected:

private:
  ///
  Integer numsubdomain;

  ///
  
  Grid<Dim> * grid;

};

template class Domain<Point3D>;
template class Domain<Point2D>;
}

#endif // FILE_DOMAIN
