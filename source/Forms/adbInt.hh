// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ADBINT
#define FILE_ADBINT

#include <string>

#include "General/defs.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "baseForm.hh" 

namespace CoupledField {


  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class ADBInt
  //!
  //! Base class for classes representing forms of type \f$ A D B \f$
  //! This class it the base class for classes representing bi-linear
  /*! forms that can be written in the form
   \f$[
   x^TA D By,\quad
   A\in\mathcal{K}^{n\times m},\quad
   D\in\mathcal{K}^{m\times q},\quad
   B\in\mathcal{K}^{q\times l}
   \f$]
   with \f$ \mathcal{K}=\mathcal{R} \f$ or \f$ \mathcal{C}\f$.
  
   The central task of this method is to provide a generalised
   implementation of the CalcElementMatrix() method. The latter uses
   the matrices \f$ A \f$, \f$ D \f$ and \f$ B \f$ to compute the local
   element matrix for the corresponding operator by approximating
   the integral over \f$ ADB \f$ using the weights and integration points
   defined by the underlying element and its base functions.
   \note In the special case of \f$ A=B^T \f$ the BDBInt class is the
   appropriate ancestor. */


class BaseMaterial;
class EntityIterator;
struct Elem;
template <class TYPE> class Matrix;

  class ADBInt : public BaseForm {

  public:

    // =======================================================================
    // CONSTRUCTION and DESTRUCTION
    // =======================================================================

    //@{ \name Construction and destruction



    //! Constructor expecting reference to a material data object
    ADBInt( BaseMaterial* matData, SubTensorType type = FULL ) 
      : BaseForm( matData, type ) {
      baseType_ = STIFFNESS;
      name_ = "ADBInt";
    };
    
    //! Destructor
    virtual ~ADBInt() {
    };

    //@}

    //! Compute element matrix associated to ADB form
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );


    // =======================================================================
    // PURE VIRTUAL METHODS
    // =======================================================================

    //! \name Pure virtual methods.
    //! These methods provide information used by the general implementation
    //! of CalcElementMatrix(). They must therefore be implemented by the
    //! derived classes.
    //@{

    //!  Compute matrix \f$ A \f$ at given integration point.
    virtual void calcAMat( Matrix<Double> &aMat, UInt ip,
                           const Matrix<Double> &ptCoord ) = 0;


    /** @see BDBInt::calcDMat(Matrix<Double>) */
    virtual void calcDMat(Matrix<Double> &dMat) 
    {
      EXCEPTION("not correctly overwritten!");
    };

    
    /** The Implement this with your SIMP variant! 
     * @see BDBInt::calcDMat(Matrix<Double>, Elem*, Double&) */
    virtual void calcDMat(Matrix<Double> &dMat, const Elem* elem)
    {
      calcDMat(dMat);
    };

    //! Query number of degrees of freedom for first physical quantity.

    //! This method can be used to query the number of degrees of freedom at
    //! any node of a finite element for the physical quantity associated to
    //! the base functions whose derivatives form the matrix \f$ A \f$.
    virtual UInt getNumDofsA() = 0;

    //! Query number of degrees of freedom for second physical quantity.

    //! This method can be used to query the number of degrees of freedom at
    //! any node of a finite element for the physical quantity associated to
    //! the base functions whose derivatives form the matrix \f$ B \f$.
    virtual UInt getNumDofsB() = 0;

    //! Query dimensions of matrix \f$ D \f$.
    virtual void getDimD( UInt nRows, UInt nCols ) = 0;

    //! Query material type for \f$ D \f$ tensor.
    virtual MaterialType getDMaterialType() = 0;


    //@}

  };


}

#endif

