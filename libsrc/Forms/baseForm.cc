#include <iostream>
#include <fstream>

#include "baseForm.hh"

namespace CoupledField
{

  BaseForm::BaseForm(BaseFE * aptelem, MaterialData & matData)
    :ptMaterial(&matData), ptelem(aptelem)
  {
#ifdef TRACE
    (*trace) << "entering BaseForm::BaseForm" << std::endl;
#endif
  }

  BaseForm::BaseForm(MaterialData & matData)
  {
#ifdef TRACE
    (*trace) << "entering BaseForm::BaseForm" << std::endl;
#endif

    ptMaterial = new MaterialData(matData);
  }

  BaseForm::BaseForm(BaseFE * aptelem)
    :ptelem(aptelem)
  {
#ifdef TRACE
    (*trace) << "entering BaseForm::BaseForm" << std::endl;
#endif
  }
 
  BaseForm::BaseForm()
  {
#ifdef TRACE
    (*trace) << "entering BaseForm::BaseForm" << std::endl;
#endif
  }
 


  BaseForm::~BaseForm()
  {
#ifdef TRACE
    (*trace) << "entering BaseForm::~BaseForm" << std::endl;
#endif
  }



  void BaseForm::CalcElementMatrix(Matrix<Double>& ptCoord, Matrix<Double>& elemMat) 
  {
#ifdef TRACE
    (*trace) <<  "entering BaseForm::CalcElemMatrix" << std::endl;
#endif
    Error(" Function BaseForm::CalcElementMatrix is virtual. You can use it for derived classes.",__FILE__,__LINE__);
  }




  void BaseForm::Print(std::ostream * out, const Matrix<Double> Result) const
  {
#ifdef TRACE
    (*trace) <<  "entering BaseForm::Print" << std::endl;
#endif
    Error(" Function BaseForm::Print is virtual. You can use it for derived classes.",__FILE__,__LINE__);
  }



}
