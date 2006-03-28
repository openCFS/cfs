#ifndef FILE_BDBINT
#define FILE_BDBINT

#include "baseForm.hh" 

namespace CoupledField {

  //! abstract class for calculation of bdb forms
  class BDBInt : public BaseForm {

  public:

    //! Constructor with pointer to BaseElem
    BDBInt();
    
    //! Constructor with pointer to BaseElem
    BDBInt(BaseFE *aptelem, BaseMaterial* matData);

    //! Constructor with pointer to BaseElem
    BDBInt(BaseMaterial* matData, SubTensorType type = FULL);

    //! Destructor
    virtual ~BDBInt();

    //! Function for calculation bdb matrix 
    virtual void CalcElementMatrix( Matrix<Double> &ptCoord,
                                    Matrix<Double> &elemmat );

    //! \note This memorial is dedicated to the undocumented function
    virtual void CalcComplexElementMatrix(Matrix<Double> &ptCoord,
                                          Matrix<Complex> &elemMat,
                                          Double &beta, Double &omega );

    virtual void GetDMat(Matrix<Double> & dMat);

    virtual void GetBMat(Matrix<Double> & bMat, Matrix<Double> & ptCoord);

  protected:

    //! returns B - matrix for BDB
    virtual void calcBMat( Matrix<Double> &bMat, UInt ip,
                           Matrix<Double> &ptCoord ) = 0;

    //! returns D - matrix for BDB
    virtual void calcDMat( Matrix<Double> &dMat ) {
      Error( "BDBInt::calcDMat(Matrix<Double>&) not correctly overwritten!",
             __FILE__, __LINE__);
    };

    //! returns D - matrix for BDB
    virtual void calcDMaterialMatWithComplexDamping( Matrix<Complex> &dMat,
                                                     Double &beta,
                                                     Double &omega ) {
      (*error) << "BDBInt::calcDMaterialMatWithComplexDamping"
               << "(Matrix<Complex> &dMat, Double &beta, Double &omega) "
               << "not correctly overwritten!";
      Error( __FILE__, __LINE__ );
    };

    //! returns D - matrix for BDB, changes in every integration point
    virtual void calcDMat( Matrix<Double> &dMat, UInt ip,
                           Matrix<Double> &ptCoord ) {
      (*error) << "BDBInt::calcDMat(Matrix<Double>&, int, Matrix<Double>&) "
               << "not correct overwritten!";
      Error( __FILE__, __LINE__ );
    };

    //! returns dimension of D matrix
    virtual UInt getDimD() = 0;

    //! returns nr. of degrees of freedom
    virtual UInt getNrDofs() = 0;

    //! Boolean for signaling that D matrix is non-constant

    //! In some cases, e.g. in non-linear computations, it may be
    //! necessary to compute the D matrix for each integration point
    //! individually. This attribute is used to signal when the latter is
    //! required.
    Boolean updateDMatInEveryIP_;
  };

}

#endif

