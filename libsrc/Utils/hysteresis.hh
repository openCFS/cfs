#ifndef FILE_HYSTERESIS
#define FILE_HYSTERESIS

#include <string>
#include "General/environment.hh"

namespace CoupledField {

class Hysteresis
{
public:
  Hysteresis(Integer numElem);

  ///
  virtual ~Hysteresis();

  //!
  virtual Double computeValue(Double xVal, Integer idxElem)
  {
    Error("computeValue not implemented in base-Class");
    return 0.0;
  };

  //!
  void updateMinMaxList(Double newX, Integer idxElem) 
  {
    Error("updateMinMaxList not implemented in base-Class");
  };


protected:

private:

  Integer numElements_;
};


} //end of namespace


#endif

