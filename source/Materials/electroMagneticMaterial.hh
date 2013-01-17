// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef ELECTROMAGNETICMATERIAL_DATA
#define ELECTROMAGNETICMATERIAL_DATA

#include <string>

#include "General/defs.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "baseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling electromagnetic material data
  */

class ApproxData;
class ElemList;

  class ElectroMagneticMaterial : public BaseMaterial {

  public:

    //! Default constructor
    ElectroMagneticMaterial();

    //! Destructor
    ~ElectroMagneticMaterial();
    
    //! Trigger finalization of material
    void Finalize();

    //! set a scalar string material parameter
    void SetScalar(const std::string& param, MaterialType matType);

    //! set a scalar real material parameter
    void SetScalar( Double param, MaterialType matType, 
		    Global::ComplexPart dataType );

    //! set a scalar complex material parameter
    void SetScalar( Complex param, MaterialType matType, 
		    Global::ComplexPart dataType );

    //! set a real material tensor
    void SetTensor(const Matrix<Double>& param, MaterialType matType,
		    Global::ComplexPart dataType );

    //! set a complex material tensor
    void SetTensor(const Matrix<Complex>& param, MaterialType matType,
		    Global::ComplexPart dataType );

    //! get a scalar real material parameter
    void GetScalar( Double& param, MaterialType matType, 
		    Global::ComplexPart dataType ) const;

    //! get a scalar complex real material parameter
    void GetScalar( Complex& param, MaterialType matType, 
		    Global::ComplexPart dataType ) const;

    //! get a real material tensor
    void GetTensor( Matrix<Double>& param, MaterialType matType,
		    Global::ComplexPart dataType,
		    SubTensorType = FULL ) const;	

    //! get a complex material tensor
    void GetTensor( Matrix<Complex>& param, MaterialType matType,
		    Global::ComplexPart dataType,
		    SubTensorType = FULL ) const;	

    //Initialize approximations of nonlinear curves
    void InitApproxCurves();

    virtual ApproxData* GetNonlinFnc( MaterialType matType ) {
      return nlinFncBH_;
    }

    virtual StdVector<ApproxData*>& GetNonlinFncs( MaterialType matType ) {     
      if ( matType == MAG_PERMEABILITYCURVES ) {
        return nlinFncBHcurves_;
      }
      else {
        EXCEPTION( "electroMagneticMaterial::GetNonlinFncs currently just for MAG_PERMEABILITYCURVES");
      }
    };

    // get ansiotropic angles
    virtual StdVector<Double>& GetAnisotropicAngles() {     
      return anisotropicAngles_;
    };

    //============================ HYSTERESIS ===================================

    //Initialize hysteresis
    virtual void InitHyst( UInt numElemSD, shared_ptr<ElemList> actSDList,
                           bool isInverse = false, bool computeInverse = false );

    //set values for differential material approach
    virtual void SetPreviousHystVal( UInt nrElem, Vector<Double>& Xval );

    //! compute scalar differential parameter
    virtual Double ComputeScalarDiffVal( UInt nrElem, Vector<Double>& Xval );

    //! compute scalar differential parameters
    virtual void ComputeScalarDiffValues( UInt nrElem, Vector<Double>& in,
                                            Vector<Double>& scalarValues );


    //! computes the scalar hystereis value
    virtual Double ComputeScalarHystVal( UInt nrElem, Vector<Double>& Xval ) {
      EXCEPTION( "ComputeScalarHystVal not implemented" );
      return 1.0; 
   };

    //! compute the vector hysteresis values
    virtual void ComputeVectorHystVal( UInt nrElem, Vector<Double>& Xin, 
                                       Vector<Double>& Yout );

    //!
    virtual void ComputeInverseScalar( UInt idxEl, UInt comp, Double Yin, 
                               Double& Xout );

    //! computes the scalar hystereis value
    virtual Double GetScalarHystVal( UInt nrElem );

    virtual void GetVectorHystVal( UInt nrElem, Vector<Double>& Val );

    Double ComputeMatDiff( Vector<Double>& dX, Vector<Double>& dY, UInt idx );

    /** overloads BaseMaterial::ComputeSubTensor() */
    void ComputeSubTensor(Matrix<Complex>& matMatrix,
			  MaterialType matType, 
			  SubTensorType subTensor) const;

  private:


    
    //! Calculate full tensor from scalar values
    void ComputeFullMuTensor();

    //! isotropic BH curve    
    ApproxData* nlinFncBH_;

    //! anisotropic BH curves
    StdVector<ApproxData*> nlinFncBHcurves_;

    UInt dim_;

    Matrix<Double> vecXprevious_; //! previous Xval in hysteresis
    Matrix<Double> vecYprevious_; //! previous Yval in hysteresis
    Matrix<Double> vecXact_; //! actual Xval in hysteresis (for inverse hysteresis)
    Matrix<Double> vecYact_; //! actual Yval in hysteresis (for inverse hysteresis)

    Vector<Double> matDiffprevious_;

    Double Xsat_, Ysat_;
  };

} // end of namespace

#endif
