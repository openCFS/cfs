#ifndef ACOUSTICPDE_BEM_HH
#define ACOUSTICPDE_BEM_HH

#include "SinglePDE.hh"

namespace CoupledField{

  // forward class declaration
  class BaseResult;
  class ResultHandler;
  class LinearFormContext;
  class ElecForceOp;
  class BaseBDBInt;
  class ResultFunctor;
  class CoefFunctionMulti;

  class AcousticPDE_BEM : public SinglePDE
  {

  };
}

#endif
