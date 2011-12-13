// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef MAGNETOSTRIMATERIAL_DATA
#define MAGNETOSTRIMATERIAL_DATA

#include <string>

#include "General/defs.hh"
#include "General/environment.hh"
#include "baseMaterial.hh"

namespace CoupledField {

  //! Class for Magnetostrictive Materials
  /*! 
     Class for handling magnetostrictive material data
  */

template <class TYPE> class Matrix;

  class MagStrictMaterial : public BaseMaterial {

  public:

    //! Default constructor
    MagStrictMaterial();

    //! Destructor
    ~MagStrictMaterial();

   //! set a scalar string material parameter
    void SetScalar(const std::string& param, MaterialType matType);

    //! set a scalar real material parameter
    void SetScalar(double param, MaterialType matType, Global::ComplexPart dataType );

    //! set a scalar complex material parameter
    void SetScalar( Complex param, MaterialType matType, Global::ComplexPart dataType );

    //! set a real material tensor
    void SetTensor( const Matrix<Double>& param, MaterialType matType, Global::ComplexPart dataType );

    //! set a complex material tensor
    void SetTensor( const Matrix<Complex>& param, MaterialType matType, Global::ComplexPart dataType );

    //! get a scalar string material parameter
    void GetScalar( std::string& param, MaterialType matType ) const;

     //! get a scalar integer material parameter
    void GetScalar( Integer& param, MaterialType matType) const;

    //! get a scalar real material parameter
    void GetScalar( Double& param, MaterialType matType, Global::ComplexPart dataType ) const;

    //! get a scalar complex real material parameter
    void GetScalar( Complex& param, MaterialType matType, Global::ComplexPart dataType ) const;

    //! get a real material tensor
    void GetTensor( Matrix<Double>& param, MaterialType matType, Global::ComplexPart dataType,
                   SubTensorType = FULL ) const;  

    //! get a complex material tensor
    void GetTensor( Matrix<Complex>& param, MaterialType matType, Global::ComplexPart dataType,
                   SubTensorType = FULL ) const;  

  private:

    //! compute the correct subTensor (3D, AXI, ..)
    void ComputeSubTensor(Matrix<Complex>& matMatrix, MaterialType matType, 
                          SubTensorType subTensor) const;

  };

} // end of namespace

#endif





