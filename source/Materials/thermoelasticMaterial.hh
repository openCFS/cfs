#ifndef THERMOELASTICMATERIAL_DATA
#define THERMOELASTICMATERIAL_DATA

#include "General/defs.hh"
#include "General/environment.hh"
#include "baseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling thermoelastic material data
  */

template <class TYPE> class Matrix;

  class ThermoelasticMaterial : public BaseMaterial {

  public:

    //! Default constructor
    ThermoelasticMaterial();

    //! Destructor
    ~ThermoelasticMaterial();


   //! set a real material tensor
    void SetTensor(const Matrix<Double>& param, MaterialType matType,
                   Global::ComplexPart dataType );


    //! get a real material tensor
    void GetTensor( Matrix<Double>& param, MaterialType matType,
                    Global::ComplexPart dataType, SubTensorType = FULL ) const;	

    //! set a scalar real material parameter
    void SetScalar( Double param, MaterialType matType, 
		    Global::ComplexPart dataType );


    //! get a scalar real material parameter
    void GetScalar( Double& param, MaterialType matType, 
                    Global::ComplexPart dataType ) const;


  private:

    //! compute the correct subTensor (3D, AXI, ..)
    void ComputeSubTensor(Matrix<Complex>& matMatrix,
			  const MaterialType& matType, 
			  const SubTensorType& subTensor) const;

  };

} // end of namespace

#endif
