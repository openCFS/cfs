#ifndef FILE_LINEARFORM
#define FILE_LINEARFORM

#include "baseform.hh"

namespace CoupledField
{

/// class for calculation right hand side
template <Integer dim>
class LinearForm : public BaseForm<dim> 
{
public:
  ///
  LinearForm(BaseElem * aptelem);

  ///
  virtual ~LinearForm();

  /// Calculation of vector of right hand side 
  void CalcElemVector(Point<dim> * ptCoord, Vector<Double> & Result);

  /// Calculation of vector of RHS such form: int(shape fnc X interpolated fnc in element)
  void CalcElemVector4InterpolatedFnc(Point<dim> * ptCoord, const Integer acom, Vector<Double> & aValueAtNodePoints, Vector<Double> & Result);

protected:

private:

};

#ifdef __GNUC__
template class  LinearForm <2>;
template class  LinearForm <3>;
#endif

}

#endif // FILE_LINEARFORM
