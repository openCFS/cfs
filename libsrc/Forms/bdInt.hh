#ifndef FILE_BDINT
#define FILE_BDINT

#include "baseForm.hh" 

namespace CoupledField {

  //! abstract class for calculation of bd forms
  class BDInt : public BaseForm {

  public:
    //! Constructor with pointer to BaseElem
    BDInt();
    
    //! Constructor with pointer to BaseElem
    BDInt(BaseFE *aptelem, BaseMaterial* matData,std::string geomType,Double timeStep);

    //! Constructor with pointer to BaseElem
    BDInt(BaseMaterial* matData,std::string geomType, Double timeStep);

    //! Destructor
    virtual ~BDInt();

    //! Function for calculation bd vector 
    virtual void calcElementVector(Matrix<Double> &ptCoord, Vector<Double> & resultStressVector,
				   Vector<Double> & fracDerivStress);



  protected:

    //! returns B - matrix for BD
    virtual void calcBMat( Matrix<Double> &bMat, Integer ip,
			   Matrix<Double> &ptCoord );
    //! returns A - matrix for BD
    //! this matrix is needed for the fractional damping in mechanics
    virtual void calcAMat( Matrix<Double> &bMat);

    //! returns Alpha - matrix for BD
    //! this matrix is needed for the fractional damping in mechanics
    virtual void calcAlphaMat( Matrix<Double> &bMat);


    //! returns D - matrix for BD
    virtual void calcDMat( Matrix<Double> &dMat ) {
      Error( "BDInt::calcDMat(Matrix<Double>&) not correctly overwritten!",
	     __FILE__, __LINE__);
    };



    //! returns D - matrix for BDB, changes in every integration point
    virtual void calcDMat( Matrix<Double> &dMat, Integer ip,
			   Matrix<Double> &ptCoord ) {
      Error( "BDInt::calcDMat(Matrix<Double>&, int, Matrix<Double>&) not correct overwritten!",__FILE__,__LINE__);
    };

    //! returns dimension of D matrix
    virtual Integer getDim();

    /// returns dimension of D matrix
    virtual Integer getDimD(){return 3;};

    //! returns nr. of degrees of freedom
    virtual Integer getNrDofs();

    //! for nonlinear calculations, the d-matrix has to be updated in every
    //! integration point
    Integer updateDMatInEveryIP_;
    std::string geomType_;

  private:
    Double dampAlpha_,  fracDeriv_, timeStepPowerFracDeriv_;

  };

}

#endif // FILE_MASSINT
