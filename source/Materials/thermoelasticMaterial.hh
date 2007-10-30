#ifndef THERMOELASTICMATERIAL_DATA
#define THERMOELASTICMATERIAL_DATA

#include "baseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling thermoelastic material data
  */

  class ThermoelasticMaterial : public BaseMaterial {

  public:

    //! Default constructor
    ThermoelasticMaterial();

    //! Destructor
    ~ThermoelasticMaterial();


   //! set a real material tensor
    void SetTensor( Matrix<Double>& param, const MaterialType& matType,
                    const DataType& dataType );


    //! get a real material tensor
    void GetTensor( Matrix<Double>& param, MaterialType& matType,
                    DataType& dataType, SubTensorType = FULL ) const;	

    //! set a scalar real material parameter
    void SetScalar( Double param, MaterialType matType, 
		    DataType dataType );


    //! get a scalar real material parameter
    void GetScalar( Double& param, MaterialType matType, 
                    DataType dataType ) const;


  private:

    //! compute the correct subTensor (3D, AXI, ..)
    void ComputeSubTensor(Matrix<Complex>& matMatrix,
			  const MaterialType& matType, 
			  const SubTensorType& subTensor) const;

  };

} // end of namespace

#endif
