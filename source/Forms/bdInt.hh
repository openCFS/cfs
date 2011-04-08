// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_BDINT
#define FILE_BDINT

#include "baseForm.hh" 

namespace CoupledField {

  //! abstract class for calculation of bd forms
  class BDInt : public BaseForm {

  public:

    //! Constructor with pointer to BaseElem
    BDInt(BaseMaterial* matData,std::string geomType, Double timeStep);

    //! Destructor
    virtual ~BDInt();

    //! Function for calculation bd vector 
    virtual void calcElementVector(Vector<Double> & resultStressVector,
                                   EntityIterator& ent,
				   Vector<Double> & fracDerivStress);



  protected:

    /** @see BaseForm::CalcBMat() */
    virtual void CalcBMat( Matrix<Double> &bMat, UInt ip, const Matrix<Double> &ptCoord );

    //! returns A - matrix for BD
    //! this matrix is needed for the fractional damping in mechanics
    virtual void calcAMat( Matrix<Double> &bMat);

    //! returns Alpha - matrix for BD
    //! this matrix is needed for the fractional damping in mechanics
    virtual void calcAlphaMat( Matrix<Double> &bMat);


    //! returns D - matrix for BD
    virtual void calcDMat( Matrix<Double> &dMat ) {
      EXCEPTION("BDInt::calcDMat(Matrix<Double>&) not correctly overwritten!");
    };



    //! returns D - matrix for BDB, changes in every integration point
    virtual void calcDMat( Matrix<Double> &dMat, Integer ip,
			   Matrix<Double> &ptCoord ) {
      EXCEPTION("BDInt::calcDMat(Matrix<Double>&, int, Matrix<Double>&) "
          << "not correct overwritten!");
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
