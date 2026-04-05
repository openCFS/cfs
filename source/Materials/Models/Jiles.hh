#ifndef FILE_JILES_2004
#define FILE_JILES_2004

#include <list>

#include "MatVec/Vector.hh"
#include "Domain/Domain.hh"
#include "Model.hh"


namespace CoupledField {

class MathParser;

class Jiles : public Model {

public:
  //! Constructor
  Jiles();

  //! Destructor
  virtual ~Jiles();

  void Init(std::map<std::string, double> ParameterMap, UInt numElems, UInt dim);

  Double ComputeMaterialParameter(Vector<Double> E, Integer ElemNum);

  Double Evaluate(Double E, Integer idx);

  void RampUp(Integer Nt, Double E, Integer idx);

  void saveValues(bool InstantSave);

private:
  //==============

  Vector<UInt> ElemNum2Idx_;

  UInt numElems_;
  Double MaxE_;

  //current index
  UInt idx_;

  Double Ps_;
  Double a_;
  Double alpha_;
  Double k_;
  Double c_;

  Vector<Double> E0_;
  Vector<Double> E1_;
  Vector<Double> E0it_;

  Vector<Double> P0_;
  Vector<Double> P1_;
  Vector<Double> P0it_;

  Vector<Double> Pi0_;
  Vector<Double> Pi1_;
  Vector<Double> Pi0it_;

  Vector<Double> Pa0_;
  Vector<Double> Pa1_;
  Vector<Double> Pa0it_;

  //! Pointer to math parser instance
  MathParser* mp_;

  Vector<Integer> isFirstTime_;
  bool isFirstTimeFinished_;

  UInt timeStep_;

  UInt globalIter_;
  double isMH_;

  std::string varHandle_;

  Vector<Double> currentDirection_;
  std::map<UInt,Vector<Double>> initialDirection_;
};
} //end of namespace

#endif

