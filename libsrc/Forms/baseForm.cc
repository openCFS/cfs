#include <iostream>
#include <fstream>

#include "baseForm.hh"

namespace CoupledField
{

  BaseForm::BaseForm(BaseFE * aptelem, MaterialData & matData)
    :ptMaterial(&matData), ptelem(aptelem), isaxi_(FALSE)
  {
    ENTER_FCN( "BaseForm::BaseForm" );
    isSetIntPoint_ = FALSE;
  }

  BaseForm::BaseForm(MaterialData & matData)
    :isaxi_(FALSE)
  {
    ENTER_FCN( "BaseForm::BaseForm" );
    isSetIntPoint_ = FALSE;
    ptMaterial = new MaterialData(matData);
  }

  BaseForm::BaseForm(BaseFE * aptelem)
    :ptelem(aptelem), isaxi_(FALSE)
  {
    ENTER_FCN( "BaseForm::BaseForm" );
    isSetIntPoint_ = FALSE;
  }
 
  BaseForm::BaseForm()
    : isaxi_(FALSE)
  {
    ENTER_FCN( "BaseForm::BaseForm" );
    isSetIntPoint_ = FALSE;
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




  void BaseForm::Print(std::ostream * out, const Matrix<Double> Result) const
  {
    ENTER_FCN( "BaseForm::Print" );
    Error(" Function BaseForm::Print is virtual. You can use it for derived classes.",__FILE__,__LINE__);
  }



}
