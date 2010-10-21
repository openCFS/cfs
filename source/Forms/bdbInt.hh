// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_BDBINT
#define FILE_BDBINT

#include "baseForm.hh"
#include "Optimization/Design/DesignElement.hh" 

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
                            EntityIterator& ent2 ){
      CalcElementMatrix( elemMat, ent1, ent2, DesignElement::NO_DERIVATIVE);
    }

    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2,
                            const DesignElement::Type direction );

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
      EXCEPTION( "BDBInt::calcGMat not implemented!");
    };

    /** Implement this in your form to return your actual physical tensor.
     * Note, BDBInt and ADBInt use calcDMat(Matrix<Double>, const Elem*)
     * which calls overwritten methods of this method but linElastInt,
     * linElecInt and linPiezoCoupling have calcDMat(Matrix<Double>, const Elem*)
     * implementations. This means, that all direct childs of these three
     * classes must provive a own version of calcDMat with the elem parameter!
     * @see calcDMat(Matrix<Double>, Elem*, Double&) */
    virtual void calcDMat(Matrix<Double> &dMat)
    {
      EXCEPTION("not correctly overwritten!");
    };


    /** This is the SIMP version, where the physical tensor [c], [\epsilon], ... is
     * multiplied with the design variable (pseudo density, pseudo polarization).
     * This could be implemented in CalcElementMatrix but this way we handle read
     * and complex case as easy as the mechanical softening implementations.
     * Note, if you want to do SIMP you have to properly implement this method!
     * This is tricky - you are warned!
     * @see linElastInt::calcDMat(Matrix<Double>, const Elem*)
     * @param dMat the output variable
     * @param elem only relevant for SIMP, the D-Mat is in general no element specific! */
    virtual void calcDMat(Matrix<Double> &dMat, const Elem* elem)
    {
      calcDMat(dMat);
    };
    /** This is the ParamMat optimization version, overwrite this to provide Derivatives for the tensor
     * used in parametric material optimization. */
    virtual void calcDMat(Matrix<Double> &dMat, const Elem* elem, DesignElement::Type direction, double force_factor = 0.0)
    {
      calcDMat(dMat, elem); 
    };


    //! returns D - matrix for BDB
    virtual void calcDMaterialMatWithComplexDamping( Matrix<Complex> &dMat,
                                                     Double &beta,
                                                     Double &omega ) {
      EXCEPTION("BDBInt::calcDMaterialMatWithComplexDamping"
               << "(Matrix<Complex> &dMat, Double &beta, Double &omega) "
               << "not correctly overwritten!" );
    };

    /** returns D - matrix for BDB, changes in every integration point
     * @see calcDMat(Matrix<Double>, EntityIterator*) */
    virtual void calcDMat( Matrix<Double> &dMat, UInt ip,
                           Matrix<Double> &ptCoord ) {
      EXCEPTION( "BDBInt::calcDMat(Matrix<Double>&, int, Matrix<Double>&) "
               << "not correct overwritten!" );
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

