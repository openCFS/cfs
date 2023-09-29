#pragma once

#include <list>

#include <boost/utility.hpp>

#include "MatVec/Vector.hh"
#include "Domain/Domain.hh"
#include "Model.hh"


namespace CoupledField {

class MathParser;

class EBHysteresis : public Model {

public:
  //! Constructor
  EBHysteresis();

  //! Destructor
  virtual ~EBHysteresis();

  void Init(std::map<std::string, double> ParameterMap, shared_ptr<ElemList> entityList, UInt dim);

  Double ComputeMaterialParameter(Vector<Double> E, Integer ElemNum);
  Matrix<Double> ComputeTensorialMaterialParameter(Vector<Double> E, Integer ElemNum);

  // just for the computation of the residual, we do not store anything here
  Vector<Double> GetFluxDensity(Vector<Double> E, Integer ElemNum);

  Matrix<Double> EvaluateLocalMu(StdVector<Double> E, StdVector<Double> D, UInt idx);

  Matrix<Double> EvaluateLocalMuDFP(StdVector<Double> E, StdVector<Double> D, UInt idx);

  Matrix<Double> EvaluateLocalMuGBM(StdVector<Double> E, StdVector<Double> D, UInt idx);

  Matrix<Double> EvaluateLocalMuFiniteDifferences(Vector<Double> E, StdVector<Double> D, UInt idx);

  Matrix<Double> EvaluateLocalMuAnhystersisOnly(Vector<Double> E, StdVector<Double> D, UInt idx);

  StdVector<Double> inv3x3(StdVector<Double> A);

  Vector<Double> Evaluate(Vector<Double> E, bool saveTmpStageVecs, UInt idx);

private:
  //==============
  UInt dim_; 

  std::map<UInt,UInt> ElemNum2Idx_;

  UInt numElems_;

  Double mu_0;

  Double Ps_;
  Double A_;
  Double mu0_;
  Double numS_;
  Double chi_factor_;

  StdVector< StdVector<Double> > H0_;
  StdVector< StdVector<Double> > H1_;

  StdVector< StdVector<Double> > M0_;
  StdVector< StdVector<Double> > M1_;

  StdVector< StdVector<Double> > HxS_prev_;
  StdVector< StdVector<Double> > HyS_prev_;
  StdVector< StdVector<Double> > HzS_prev_;
  StdVector< StdVector<Double> > MxS_prev_;
  StdVector< StdVector<Double> > MyS_prev_;
  StdVector< StdVector<Double> > MzS_prev_;

  StdVector< StdVector<Double> > HxS_n_;
  StdVector< StdVector<Double> > HyS_n_;
  StdVector< StdVector<Double> > HzS_n_;
  StdVector< StdVector<Double> > MxS_n_;
  StdVector< StdVector<Double> > MyS_n_;
  StdVector< StdVector<Double> > MzS_n_;

  StdVector< StdVector<Double> > HxS_n_tmp_;
  StdVector< StdVector<Double> > HyS_n_tmp_;
  StdVector< StdVector<Double> > HzS_n_tmp_;
  StdVector< StdVector<Double> > MxS_n_tmp_;
  StdVector< StdVector<Double> > MyS_n_tmp_;
  StdVector< StdVector<Double> > MzS_n_tmp_;

  // epsilon tensor of the previous iteration
  StdVector< Matrix<Double> > mu_;

  //! Pointer to math parser instance
  MathParser* mp_;

  Vector<Integer> isFirstTime_;
  bool isFirstTimeFinished_;

  UInt timeStep_;

  UInt globalIter_;
  double isMH_;

  std::string varHandle_;

  StdVector<bool> hasElemSolution_;

};
} //end of namespace