#ifndef FILE_JILES_2004
#define FILE_JILES_2004

#include <list>

#include "Utils/vector.hh"
#include "hysteresis.hh"

namespace CoupledField {

class Jiles : public Hysteresis
{
public:
  Jiles(Integer numElem, Double ysat, Double a, Double alpha, 
	Double k, Double c);

  //!
  virtual ~Jiles();

  //!
  Double computeValue(Double xVal, Integer idxElem);

  //!
  void updateMinMaxList(Double newX, Integer idxElem);

  //! 
  void SetTimeStepVal(Double dt) 
  {dt_ = dt; };

protected:

private:

  Double Ysaturated_;
  Double a_;
  Double alpha_;
  Double k_;
  Double c_;

  Double dt_;
  Integer actElem_;

  Vector<Double> Xold_;
  Vector<Double> YirrOld_;
  Vector<Double> YirrNew_;
  Vector<Double> YirrPrev_;

  Double iterMax_;
  Double err_;
};


} //end of namespace


#endif

