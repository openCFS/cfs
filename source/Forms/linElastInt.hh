// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LINELASTINT
#define FILE_LINELASTINT

#include "Forms/bdbInt.hh"
#include "General/environment.hh"

#include "General/defs.hh"
#include "Optimization/Design/DesignElement.hh"

namespace CoupledField {
  
  
  //! base class for calculation of linear elasticity
class BaseFE;
class BaseMaterial;
class EntityIterator;
struct Elem;
template <class TYPE> class Matrix;

  class linElastInt : public BDBInt {

  public:

    //! Constructor
    linElastInt( BaseMaterial* matData, SubTensorType type = FULL );

    //! Destructor
    virtual ~linElastInt();

    /** Function for calculation bdb matrix. 
     * In ersatz material optimization we the form \int B D B where D is [c] (lin elast in)
     * D might be multiplied with an ersatz-Material factor. This is handled in the by calcDMat() */
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1,
                            EntityIterator& ent2,
                            const DesignElement::Type direction);

    void CalcElementMatrix( Matrix<Complex>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2,
                            const DesignElement::Type direction);

    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1,
                            EntityIterator& ent2 ){
      CalcElementMatrix( elemMat, ent1, ent2, DesignElement::NO_DERIVATIVE);
    }

    void CalcElementMatrix( Matrix<Complex>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 ){
      CalcElementMatrix( elemMat, ent1, ent2, DesignElement::NO_DERIVATIVE);
    }

    /** Function for calculation bdb matrix using incompatible modes. Uses */ 
    void CalcElementMatrixICM( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );

    //! Function for calculation bdb matrix of shear terms using BK idea 
    void CalcElementMatrixShearBK1( Matrix<Double>& elemMat,
				    EntityIterator& ent1, 
				    EntityIterator& ent2 );

    //! does the reordering that is done in calculation of the B matrix, this is needed by Shape Derivative
    void ReorderBLikeMatrix( Matrix<Double>& in, Matrix<Double>& out, UInt ip, BaseFE* elem, const Matrix<Double> &ptCoord );

    //! Query material type for \f$D\f$ tensor
    MaterialType getDMaterialType() { return MECH_STIFFNESS_TENSOR; }

    /** Calculates the Material data in the SIMP version!
     * Or in the ParamMat version 
     * Note, that most forms implement this method in the non-SIMP variant with no
     * elem parameter but the BDBInt uses calcDMat with the elem parameter.
     * Therefore all direct childs of this class should either give calcDMat(dMat, elem)
     * or provide a variant of that class the non-elem version.
     * Take care! This is error prone!! 
     * For ParamMat case, if direction is NO_DERIVATIVE, we fall back to the 2 parameter 
     * call (in bdbInt::CalcElementMatrix), so that child classes may also only implement 2 paramater version
     * if no ParamMat optimization is needed
     * @param direction if !=  DesignElement::NO_DERIVATIVE calculate derivative instead
     * @param force_factor if set no pseudo density factor is applied but the force_factor (for direction = DENSITY) */
    void calcDMat(Matrix<Double> &dMat, const Elem* elem, const DesignElement::Type direction, double force_factor = 0.0);

    /** see calcDMat(Matrix<Double>, const Elem*, const DesignElement::Type, force_factor) */
    virtual void calcDMat( Matrix<Double> &dMat, const Elem* elem){
      calcDMat( dMat, elem, DesignElement::NO_DERIVATIVE);
    }

  protected:
    
    /** see calcDMat(Matrix<Double>, const Elem*) */
    virtual void calcDMat( Matrix<Complex> &dMat, const Elem* elem);

    /** @see BaseForm::CalcBMat() */
    void CalcBMat( Matrix<Double> &bMat, UInt ip, const Matrix<Double> &ptCoord );

    //! returns G - matrix for GDG (incompatible modes)
    virtual void calcGMat( Matrix<Double> &bMat, UInt ip,
                           Matrix<Double> &ptCoord );

    //! set dimensions
    virtual void SetDimensions(SubTensorType type);

    //! returns dimension of D matrix
    virtual UInt getDimD() {
      return dimD_; 
    };
    
    //! returns nr. of degrees of freedom
    virtual UInt getNrDofs() {
      return nrDofs_;
    };
    
    /** fetch the ErsatzMaterialTensor from the DesignSpace when using ParamMat optimization
     * @param t output variable, returns the Tensor
     * @elem current element
     * @direction if given return derivative w.r.t. this designelement type
     * @return true if the ersatz material tensor has been set */
    bool GetErsatzMaterialTensor(Matrix<double>& t, const Elem* elem, DesignElement::Type direction);
    
  private:

    //dimension of Dmatrix
    UInt dimD_;
    
    //! number of degrees 
    UInt nrDofs_;

    //! subtype of tensor
    SubTensorType subTensorType_;

    //! number of incompatible modes
    UInt nrICModes_;
    
  };
  



} //end namespace

// Explicit template instantiation
//#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
//template void linElastInt::CalcElementMatrix( Matrix<Double>& elemMat,
//    EntityIterator& ent1,
//    EntityIterator& ent2,
//    const DesignElement::Type direction);
//template void linElastInt::CalcElementMatrix( Matrix<Complex>& elemMat,
//    EntityIterator& ent1,
//    EntityIterator& ent2,
//    const DesignElement::Type direction);
//#endif

#endif // FILE_LINELASTINT
