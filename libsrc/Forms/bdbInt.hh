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
    BDBInt(BaseFE *aptelem, MaterialData &matData);

    //! Constructor with pointer to BaseElem
    BDBInt(MaterialData &matData);

    //! Destructor
    virtual ~BDBInt();

    //! Function for calculation bdb matrix 
    virtual void CalcElementMatrix(Matrix<Double> &ptCoord, Matrix<Double> & elemmat);

    virtual void CalcComplexElementMatrix(Matrix<Double>& ptCoord, Matrix<Complex> & elemMat,Double & beta, Double & omega);

  protected:

    //! returns B - matrix for BDB
    virtual void calcBMat( Matrix<Double> &bMat, Integer ip,
			   Matrix<Double> &ptCoord ) = 0;

    //! returns D - matrix for BDB
    virtual void calcDMat( Matrix<Double> &dMat ) {
      Error( "BDBInt::calcDMat(Matrix<Double>&) not correctly overwritten!",
	     __FILE__, __LINE__);
    };

    //! returns D - matrix for BDB
    virtual void calcDMaterialMatWithComplexDamping(Matrix<Complex> &dMat,Double &beta,Double &omega){
      Error( "BDBInt::calcDMaterialMatWithComplexDamping(Matrix<Complex> &dMat, Double &beta, Double &omega) not correctly overwritten!",
	     __FILE__, __LINE__);
    };


    //! returns D - matrix for BDB, changes in every integration point
    virtual void calcDMat( Matrix<Double> &dMat, Integer ip,
			   Matrix<Double> &ptCoord ) {
      Error( "BDBInt::calcDMat(Matrix<Double>&, int, Matrix<Double>&) not correct overwritten!",__FILE__,__LINE__);
    };

    //! returns dimension of D matrix
    virtual Integer getDimD()=0;

    //! returns nr. of degrees of freedom
    virtual Integer getNrDofs()=0;

    //! for nonlinear calculations, the d-matrix has to be updated in every
    //! integration point
    Integer updateDMatInEveryIP_;
  };

}

#endif // FILE_MASSINT
