#ifndef FILE_COMPLEXMATRIX
#define FILE_COMPLEXMATRIX

namespace CoupledField
{

class ComplexMatrix : public BaseMatrix
{
public:
  ///
  ComplexMatrix();

  ///
  virtual ~ComplexMatrix();

private:
  ///
  BaseMatrix * sysmat[2];

};

}

#endif // FILE_COMPLEXMATRIX
