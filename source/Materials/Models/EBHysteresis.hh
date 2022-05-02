#pragma once

#include <list>

#include <boost/utility.hpp>

#include "MatVec/Vector.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Domain/Domain.hh"
#include "Model.hh"


namespace CoupledField {

class EBHysteresis : public Model {

public:
  //! Constructor
  EBHysteresis();

  //! Destructor
  virtual ~EBHysteresis();

  void Init(std::map<std::string, double> ParameterMap, UInt numElems);

  Double ComputeMaterialParameter(Vector<Double> E, Integer ElemNum);
  Matrix<Double> ComputeTensorialMaterialParameter(Vector<Double> E, Integer ElemNum);

  Matrix<Double> EvaluateLocalEpsilon(StdVector<Double> E, StdVector<Double> D);

  Matrix<Double> EvaluateLocalEpsilonDFP(StdVector<Double> E, StdVector<Double> D);

  Vector<Double> Evaluate(Vector<Double> E);

  //  void RampUp(Integer Nt, Double E, Integer idx);

  void saveValues(bool InstantSave);

private:
  //==============

  Vector<UInt> ElemNum2Idx_;

  UInt numElems_;
  Double MaxE_;

  //current index
  UInt idx_;

  Double epsilon0_;


  Double Ps_;
  Double A_;
  Double mu0_;
  Double numS_;
  Double chi_factor_;

  StdVector< StdVector<Double> > E0_;
  StdVector< StdVector<Double> > E1_;

  StdVector< StdVector<Double> > P0_;
  StdVector< StdVector<Double> > P1_;

  StdVector< StdVector<Double> > HxS_prev_;
  StdVector< StdVector<Double> > HyS_prev_;
  StdVector< StdVector<Double> > MxS_prev_;
  StdVector< StdVector<Double> > MyS_prev_;

  StdVector< StdVector<Double> > HxS_n_;
  StdVector< StdVector<Double> > HyS_n_;
  StdVector< StdVector<Double> > MxS_n_;
  StdVector< StdVector<Double> > MyS_n_;

  StdVector< StdVector<Double> > HxS_n_tmp_;
  StdVector< StdVector<Double> > HyS_n_tmp_;
  StdVector< StdVector<Double> > MxS_n_tmp_;
  StdVector< StdVector<Double> > MyS_n_tmp_;

  // epsilon tensor of the previous iteration
  StdVector< Matrix<Double> > epsilon_;

  //! Pointer to math parser instance
  MathParser* mp_;

  Vector<Integer> isFirstTime_;
  bool isFirstTimeFinished_;

  UInt timeStep_;

  UInt globalIter_;
  double isMH_;

  std::string varHandle_;

};
} //end of namespace
