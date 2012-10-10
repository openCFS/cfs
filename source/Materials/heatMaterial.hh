// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef HEATMATERIAL_DATA
#define HEATMATERIAL_DATA

#include "General/defs.hh"
#include "General/environment.hh"
#include "baseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling heat material data
  */

class ApproxData;
template <class TYPE> class Matrix;

  class HeatMaterial : public BaseMaterial {

  public:

    //! Default constructor
    HeatMaterial();

    //! Destructor
    ~HeatMaterial();

    //! set a scalar real material parameter
    void SetScalar( Double param, MaterialType matType, 
		    Global::ComplexPart dataType );

    //! set a scalar complex material parameter
    void SetScalar( Complex param, MaterialType matType, 
		    Global::ComplexPart dataType );

    //! set a real material tensor
    void SetTensor(const Matrix<Double>& param, 
                   MaterialType matType, 
                   Global::ComplexPart  dataType );


    //! set a complex material tensor
    // void SetTensor(const Matrix<Complex>& param, MaterialType matType,  Global::ComplexPart dataType )
    // {
    // matTypeNotAllowed( matType, "tensor");
    //}

    //    //! set a complex material tensor
    //    void SetTensor( Matrix<Complex>& param, const MaterialType& matType,
    //		    const Global::ComplexPart& dataType );

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

    //! get correct nonlinear function
    virtual ApproxData* GetNonlinFnc( MaterialType matType ) {
      if ( matType == HEAT_CONDUCTIVITY )
        return nlinFncConductivity_;
      else 
        return  nlinFncCapacity_;
    }

    //Initialize approximations of nonlinear curves
    void InitApproxCurves();

    /** overloads BaseMaterial::ComputeSubTensor() */
    void ComputeSubTensor(Matrix<Complex>& matMatrix,
			  const MaterialType& matType, 
			  const SubTensorType& subTensor) const;
    
  private:

    ApproxData* nlinFncConductivity_;
    ApproxData* nlinFncCapacity_;
  };

} // end of namespace

#endif
