#ifndef FILE_COMPLEXMATRIX_CLA
#define FILE_COMPLEXMATRIX_CLA

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

#endif // FILE_COMPLEXMATRIX_CLA
