#ifndef ELECTROSTATICMATERIAL_DATA
#define ELECTROSTATICMATERIAL_DATA

#include "BaseMaterial.hh"

namespace CoupledField {


  //! Class for Material Data
  /*! 
     Class for handling electrostatic material data
  */

  class Hysteresis;

  class ElectroStaticMaterial : public BaseMaterial {

  public:

    //! Default constructor
    ElectroStaticMaterial(MathParser* mp,
                          CoordSystem * defaultCoosy);

    //! Destructor
    ~ElectroStaticMaterial();

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

    //! get a string material parameter
    void GetScalar( std::string& param, 
		    MaterialType matType) const;

    //! get a scalar real material parameter
    void GetScalar( Double& param, MaterialType matType, 
		    Global::ComplexPart dataType ) const;

    //! get a scalar complex real material parameter
    void GetScalar( Complex& param, MaterialType matType, 
		    Global::ComplexPart dataType ) const;


    //! get a scalar integer material parameter
    void GetScalar( Integer& param, 
                    MaterialType matType ) const;
    
    //! get a real material tensor
    void GetTensor( Matrix<Double>& param, MaterialType matType,
		    Global::ComplexPart dataType,
		    SubTensorType = FULL ) const;	

    //! get a complex material tensor
    void GetTensor( Matrix<Complex>& param, MaterialType matType,
		    Global::ComplexPart dataType,
		    SubTensorType = FULL ) const;	

    //! Finalize material
    void Finalize();
    
    //======================== for hysteresis =======================================
    //! compute scalar differential parameter
    Double ComputeScalarDiffVal( UInt nrElem, Double Xval );


    //set values for differential material approach
    void SetPreviousHystVal( UInt nrElem, Double Xval );

    Double ComputeScalarHystVal( UInt nrElem, Double Xval );


  private:

    //! compute the correct subTensor (3D, AXI, ..)
    void ComputeSubTensor(Matrix<Complex>& matMatrix,
			  MaterialType matType, 
			  SubTensorType subTensor) const;
    
    //! compute the full permittivity tensor in case of istoropic material
    void ComputeFullEpsTensor();

    Complex scalarPermittivity_;

  };

} // end of namespace

#endif
