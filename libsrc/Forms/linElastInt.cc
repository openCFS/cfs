#include <iostream>
#include <fstream>

#include "linElastInt.hh"

namespace CoupledField
{

  mechPlainStrainInt::mechPlainStrainInt(BaseFE * aptelem) 
    : BaseForm(aptelem)
  {
#ifdef TRACE
    (*trace) << "entering mechPlainStrainInt::mechPlainStrainInt" << std::endl;
#endif

    ptelem=aptelem;
  }
 

  mechPlainStrainInt::~mechPlainStrainInt()
  {
#ifdef TRACE
    (*trace) << "entering mechPlainStrainInt::~mechPlainStrainInt" << std::endl;
#endif
  }



  void mechPlainStrainInt::CalcElementMatrix(Matrix<Double>& ptCoord, Matrix<Double>& StiffMat) 
  {
#ifdef TRACE
    (*trace) <<  "entering mechPlainStrainInt::CalcElemMatrix" << std::endl;
#endif
    Error(" Function mechPlainStrainIntCalcElementMatrix is virtual. You can use it for derived classes.",__FILE__,__LINE__);
  }




  void mechPlainStrainInt::Print(std::ostream * out, const Matrix<Double> Result) const
  {
#ifdef TRACE
    (*trace) <<  "entering mechPlainStrainInt::Print" << std::endl;
#endif
    Error(" Function mechPlainStrainInt::Print is virtual. You can use it for derived classes.",__FILE__,__LINE__);
  }



}
