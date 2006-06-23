#ifndef FILE_SMOOTHINT
#define FILE_SMOOTHINT

#include <Elements/basefe.hh>
#include <Forms/bdbInt.hh>
#include <Materials/baseMaterial.hh>

namespace CoupledField
{
  
  
  /// base class for calculation of linear elasticity
class SmoothInt : public BDBInt
{
public:

  //! Constructor
  SmoothInt(BaseMaterial* matData, SubTensorType type = FULL,
            bool coordUpdate = false );
  
  //! Destructor
  virtual ~SmoothInt();
  
protected:    
  
  //! returns B - matrix for BDB
  virtual void calcBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord);

  //! calculate the data-matrix for 2D plain-strain
  virtual void calcDMat(Matrix<Double> & dMat, UInt ip, Matrix<Double> & ptCoord);

  //! returns dimension of D matrix
  virtual UInt getDimD() {
    return dimD_; 
  };
  
  //! returns nr. of degrees of freedom
  virtual UInt getNrDofs() {
    return nrDofs_;
  };

  private:

  //dimension of Dmatrix
  UInt dimD_;

  //! number of degrees 
  UInt nrDofs_;
};
  
} //end namespace

#endif // FILE_SMOOTHINT
