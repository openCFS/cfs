#include <iostream>
#include <fstream>

#include "baseForm.hh"

namespace CoupledField
{

  BaseForm::BaseForm(BaseFE * aptelem, MaterialData & matData)
    :ptMaterial(&matData), ptelem(aptelem), isaxi_(FALSE), dofzero_(0), isFracDamping_(FALSE)
  {
    ENTER_FCN( "BaseForm::BaseForm" );
    isSetIntPoint_ = FALSE;
    piezoMatType_ = realMaterialParameter;

  }

  BaseForm::BaseForm(MaterialData & matData)
    :isaxi_(FALSE), dofzero_(0), isFracDamping_(FALSE)
  {
    ENTER_FCN( "BaseForm::BaseForm" );
    isSetIntPoint_ = FALSE;
    ptMaterial = new MaterialData(matData);
    piezoMatType_ = realMaterialParameter;
    
  }

  BaseForm::BaseForm(BaseFE * aptelem)
    :ptelem(aptelem), isaxi_(FALSE),  dofzero_(0), isFracDamping_(FALSE)
  {
    ENTER_FCN( "BaseForm::BaseForm" );
    isSetIntPoint_ = FALSE;
    piezoMatType_ = realMaterialParameter;
  }
 
  BaseForm::BaseForm()
    : isaxi_(FALSE),  dofzero_(0), isFracDamping_(FALSE)
  {
    ENTER_FCN( "BaseForm::BaseForm" );
    isSetIntPoint_ = FALSE;
    piezoMatType_ = realMaterialParameter;
  }
 


  BaseForm::~BaseForm()
  {
    ENTER_FCN( "BaseForm::~BaseForm" );
  }



  void BaseForm::CalcElementMatrix(Matrix<Double>& ptCoord, Matrix<Double>& elemMat) 
  {
    ENTER_FCN( "BaseForm::CalcElemMatrix" );
    Error(" Function BaseForm::CalcElementMatrix is virtual. You can use it for derived classes.",__FILE__,__LINE__);
  }


  void BaseForm::CalcComplexElementMatrix(Matrix<Double>& ptCoord, Matrix<Complex>& elemMat, Double & beta, Double & omega) 
  {
    ENTER_FCN( "BaseForm::CalcComplexElemMatrix" );
    Error(" Function BaseForm::CalcComplexElementMatrix is virtual. You can use it for derived classes.",__FILE__,__LINE__);
  }


  void BaseForm::Print(std::ostream * out, const Matrix<Double> Result) const
  {
    ENTER_FCN( "BaseForm::Print" );
    Error(" Function BaseForm::Print is virtual. You can use it for derived classes.",__FILE__,__LINE__);
  }



}
