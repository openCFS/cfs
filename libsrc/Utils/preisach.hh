#ifndef FILE_PREISACH_2004
#define FILE_PREISACH_2004

#include <list>

#include "Utils/vector.hh"
#include "hysteresis.hh"

namespace CoupledField {

class Preisach : public Hysteresis
{
public:
  Preisach(Integer numElem, Double xSat, Double ysat, Double xRem,
	   Boolean isVirgin);

  //!
  virtual ~Preisach();

  //!
  Double computeValue(Double xVal, Integer idxElem);

  //!
  void updateMinMaxList(Double newX, Integer idxElem);

  //!
  void wipout(Integer idx);

  //!
  Double everett(Double x1, Double x2);

  //!
  Double normalizeInput(Double xInput);

protected:

private:

  Double Xsaturated_;
  Double YSaturated_;
  Double YRemnant_;

  Boolean isVirgin_;
  Integer actElem_;

  Vector<Double> lastVal_;
  Vector<Double> preisachSum_;

  std::vector<Integer> *isMinMax_;
  std::vector<Double>  *extremaList_;


};


} //end of namespace


#endif

