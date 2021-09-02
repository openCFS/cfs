#include <iostream>
#include <fstream>
#include <cmath>

#include "DataInOut/Logging/LogConfigurator.hh"

#include "Jiles.hh"
#include "Model.hh"

#include "MatVec/Vector.hh"

#include "Utils/mathParser/mathParser.hh"
#include "Domain/Domain.hh"

namespace CoupledField {

DEFINE_LOG(ja, "Jiles")

Jiles::Jiles(): Model() {
}

Jiles::~Jiles() {
}

void Jiles::Init(std::map<std::string, double> ParameterMap, UInt numElems) {

  numElems_ = numElems;

  ElemNum2Idx_.Resize(numElems);
  ElemNum2Idx_.Init(0);
  for(UInt i = 0; i < numElems;i++){
    ElemNum2Idx_[i]=i;
  }

  if (ParameterMap.size() != 5) {
    EXCEPTION("The model needs 5 parameter!");
  }
  Ps_ = ParameterMap["Ps"];
  alpha_ = ParameterMap["alpha"];
  a_ = ParameterMap["a"];
  k_ = ParameterMap["k"];
  c_ = ParameterMap["c"];

  MaxE_ = 0;

  E0_.Resize(numElems_);
  E0_.Init(0);

  E1_.Resize(numElems_);
  E1_.Init(0);

  Pi0_.Resize(numElems_);
  Pi0_.Init(0);

  Pi1_.Resize(numElems_);
  Pi1_.Init(0);

  Pa0_.Resize(numElems_);
  Pa0_.Init(0);

  Pa1_.Resize(numElems_);
  Pa1_.Init(0);

  P0_.Resize(numElems_);
  P0_.Init(0);

  P1_.Resize(numElems_);
  P1_.Init(0);


  isFirstTime_.Resize(numElems_);
  isFirstTime_.Init(1);

  isFirstTimeFinished_ = false;

  currentDirection_.Resize(3);
  currentDirection_.Init(0);

  initialDirection_.clear();

  timeStep_ = 1;

  mp_ = domain->GetMathParser();

}

Double Jiles::ComputeMaterialParameter(Vector<Double> EVec, const Integer ElemNum) {

  idx_=ElemNum2Idx_[ElemNum-1];

  saveValues(false);

  Double E, P;

  E = EVec.NormL2();
  if(E == 0){
    return 4.028353e-07;
  }

  currentDirection_ = EVec;
  currentDirection_.ScalarDiv(E);

  if (!isFirstTimeFinished_ && (isFirstTime_[idx_])) {

    //Register intial direction
    initialDirection_[idx_]=currentDirection_;

    //Evaluate the anhysteretic curve to give a better starting point for P
    P0_[idx_] = Ps_ * (cosh(E / a_) / sinh(E / a_) - a_ / E);

    //Ramp up, in case signal doesnt start at 0, to provide stable JA-solution.
    RampUp(512, E, idx_);

    //Set first time for element to false
    isFirstTime_[idx_] = 0;

    if (idx_==numElems_-1){

      // Double checking if each element is initalized once. eg each element has its index.
      for(UInt i = 0; i < isFirstTime_.GetSize();i++){
        if(isFirstTime_[i]){
          EXCEPTION("This should not happen!");
        }
      }
      isFirstTimeFinished_=true;
    }
  }

  //Check if direction has changed
//  double innerProduct = 0;
//  innerProduct = currentDirection_.Inner(initialDirection_[idx_]);
//  if (abs(innerProduct)<0.95){
//    EXCEPTION("Dependend Vector rotates to much. This should not occur. This is a scalar hysteresis model.");
//  }

  P = Evaluate(E, idx_);

  Double epsilon, epsilon0;
  epsilon0 = 8.854187e-12;

  if(E == 0 ){
    epsilon = 4.028353e-07;
  } else{
    epsilon = epsilon0 + P / E;
  }

  if(isinf(epsilon) || isnan(epsilon)){
    std::cout << "E0: "<< E0_[idx_] << std::endl;
    std::cout << "E1: "<< E1_[idx_]<< std::endl;
    std::cout << "P0: "<< P0_[idx_] << std::endl;
    std::cout << "P1: "<< P1_[idx_]<< std::endl;
    std::cout << "Pa0: "<< Pa0_[idx_] << std::endl;
    std::cout << "Pa1: "<< Pa1_[idx_]<< std::endl;
    std::cout << "Pi0: "<< Pi0_[idx_] << std::endl;
    std::cout << "Pi1: "<< Pi1_[idx_]<< std::endl;

    EXCEPTION("Epsilon is inifite or NaN...")
  }

  return epsilon;
}

void Jiles::RampUp(Integer Nt, Double E, Integer idx) {
//  LOG_DBG(ja) << "Computing Ramp Up-Polarization for element " << idxElem ;

  Double tmp_E;
  tmp_E = E / Nt;
  for (int i = 1; i < Nt-1; i++) {
    Evaluate(tmp_E * i, idx);
    saveValues(true);
  }
}
double Jiles::Evaluate(Double E, Integer idx) {

  // Currently only forward integration...
  double Pa, Ee, P, Pi, dE, delta;    //, delta2, dPPa;

  dE = E - E0_[idx];

//  std::cout << "E1 = " << E << std::endl;
//  std::cout << "E0 = " << E0_[idx] << std::endl;
  if(dE == E){
    dE = E - 0.99*E;
  }

  if(dE==0){
    delta = 1;
  }else{
    delta = dE / std::abs(dE);
  }
  //dPPa = P0_[idx] - Pa0_[idx];

  //Effective electric field
  Ee = E + alpha_ * P0_[idx];

  //anhysteretic polarisation
  Pa = Ps_ * (cosh(Ee / a_) / sinh(Ee / a_) - a_ / Ee);


  //irreversible polarisation
  Pi = Pi0_[idx] + (Pa - Pi0_[idx]) * dE / (delta * k_ - alpha_ * (Pa - Pi0_[idx]));

  // total polarization
  P = c_ * Pa + (1 - c_) * Pi;

  // this are the latest values...
  Pa1_[idx] = Pa;
  Pi1_[idx] = Pi;
  P1_[idx] = P;
  E1_[idx] = E;


  return P;
}

void Jiles::saveValues(bool InstantSave){
  if((timeStep_ != mp_->GetExprVars(MathParser::GLOB_HANDLER, "step")) || InstantSave){
    Pa0_=Pa1_;
    Pi0_=Pi1_;
    P0_=P1_;
    E0_=E1_;

    Pa1_.Init(0);
    Pi1_.Init(0);
    P1_.Init(0);
    E1_.Init(0);

    timeStep_ = mp_->GetExprVars(MathParser::GLOB_HANDLER, "step");
  }

}

}

