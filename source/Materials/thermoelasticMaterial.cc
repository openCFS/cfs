#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "thermoelasticMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  ThermoelasticMaterial::ThermoelasticMaterial() : BaseMaterial() {

    materialDatabaseName_ = "Thermoelastic";

    //set the allowed material parameters
	isAllowed_.insert( THERMAL_EXPANSION_TENSOR );
  }

  ThermoelasticMaterial::~ThermoelasticMaterial() {

    materialDatabaseName_ = "Thermoelastic";

  }


//   void  ThermoelasticMaterial::SetScalar( std::string& param, const MaterialType& matType) {
// 
//     EXCEPTION("SetScalar for 'String' not implemented");
//   }
// 
// 
   void ThermoelasticMaterial::SetScalar( Double param,MaterialType matType, 
                                          Global::ComplexPart dataType ) { 
     EXCEPTION("SetScalar for 'Double' not implemented");}
 

  void ThermoelasticMaterial::SetTensor(const Matrix<Double>& param,
                                        MaterialType matType, 
                                        Global::ComplexPart dataType ) {
    

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "tensor";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );
      if ( dataType == Global::REAL || dataType == Global::IMAG ) {
	if ( tensorParams_[matType].GetNumRows() == 0 ) {
	  tensorParams_[matType].Resize( param.GetNumRows(), param.GetNumCols() );
          tensorParams_[matType].Init();
	}
	if ( tensorParamsOrig_[matType].GetNumRows() == 0 ) {
	  tensorParamsOrig_[matType].Resize( param.GetNumRows(), param.GetNumCols() );
          tensorParamsOrig_[matType].Init();
	}

	tensorParams_[matType].SetPart( dataType, param );
	tensorParamsOrig_[matType].SetPart( dataType, param );

	// to be consistent to old structure
	if ( dataType == Global::REAL ) {
	  scalarParams_[matType] = Complex( param[2][2], 0.0);
	}
	else {
	  scalarParams_[matType] = Complex( 0.0, param[2][2]);
	  isComplex_.insert( matType );
	}
      }
      else {
	std::string msg = "SetTensor-Double";
	dataTypeNotAllowed4SetGet ( dataType, msg );
      }
    }


  }

//   void ThermoelasticMaterial::SetTensor( Matrix<Complex>& param, const MaterialType& matType, 
// 					   const Global::ComplexPart& dataType ) {
//     
//     ENTER_FCN( "ThermoelasticMaterial::SetTensor" );
//     std::string dim = "tensir";
//     matTypeNotAllowed( matType, dim );
//     
//   }


  void ThermoelasticMaterial::GetScalar( Double& param,
                                         MaterialType matType, 
                                         Global::ComplexPart dataType ) const {
    EXCEPTION("GetScalar for 'Double' not implemented");
  }


  void ThermoelasticMaterial::GetTensor( Matrix<Double>& param, 
                                         MaterialType matType, 
                                         Global::ComplexPart dataType,
                                         SubTensorType subTensor) const {


    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    if ( pos == tensorParams_.end() ) {
      std::string dim = "tensor";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      Matrix<Complex> matTensor;
      if ( subTensor == FULL || subTensor == AXI ) {
    	  matTensor = pos->second;
      }
      else {
    	  ComputeSubTensor(matTensor, matType, subTensor);
      }

      if ( dataType == Global::REAL || dataType == Global::IMAG) {
    	  param = matTensor.GetPart( dataType );
      }
      else {
		std::string msg = "GetTensor-Double";
		dataTypeNotAllowed4SetGet( dataType, msg );
      }
    }
  }

  void ThermoelasticMaterial::ComputeSubTensor(Matrix<Complex>& matMatrix,
					       const MaterialType& matType, 
					       const SubTensorType& subTensor) const {

    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    //2D tensor for the plane xy
    matMatrix.Resize(2,2);
    matMatrix.Init();
    pos->second.GetSubMatrix(matMatrix, 0, 0);
  }

}
