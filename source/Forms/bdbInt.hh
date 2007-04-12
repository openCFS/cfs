// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_BDBINT
#define FILE_BDBINT

#include "baseForm.hh" 

namespace CoupledField {

  //! abstract class for calculation of bdb forms
  class BDBInt : public BaseForm {

  public:

    //! Constructor with pointer to BaseElem
    BDBInt(BaseMaterial* matData, SubTensorType type = FULL,
           bool coordUpdate = false );

    //! Destructor
    virtual ~BDBInt();

    //! Compute element matrix associated to ADB form
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
    
    //! \note This memorial is dedicated to the undocumented function
    void CalcComplexElementMatrix( Matrix<Complex> & elemMat,
                                   EntityIterator& ent1, 
                                   EntityIterator& ent2,
                                   Double & beta, Double & omega );

  

    //! Get B-Matrix of element midpoint
    virtual void calcBMat(EntityIterator it, Matrix<Double> & bMat);
    
    //! Get DB-Matrix of element midpoint
    virtual void calcDBMat( EntityIterator it,
                    Matrix<Double> & bMat );


    //! returns B - matrix for BDB
    virtual void calcBMat( Matrix<Double> &bMat, UInt ip,
                           Matrix<Double> &ptCoord ) = 0;

    //! returns G - matrix for GDG (incompatible modes)
    virtual void calcGMat( Matrix<Double> &bMat, UInt ip,
                           Matrix<Double> &ptCoord ) {
      Error( "BDBInt::calcGMat not implemented!",
             __FILE__, __LINE__);
    };

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

    //! Query material type for \f$D\f$ tensor
    virtual MaterialType getDMaterialType() = 0;

  protected:
    
    //! bool for signaling that D matrix is non-constant

    //! In some cases, e.g. in non-linear computations, it may be
    //! necessary to compute the D matrix for each integration point
    //! individually. This attribute is used to signal when the latter is
    //! required.
    bool updateDMatInEveryIP_;
  };

}

#endif

