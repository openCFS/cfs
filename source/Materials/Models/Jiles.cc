#include <iostream>
#include <fstream>
#include <cmath>

#include "DataInOut/Logging/LogConfigurator.hh"

#include "Jiles.hh"
#include "Model.hh"

#include "MatVec/Vector.hh"

#include "Utils/mathParser/mathParser.hh"
#include "Domain/Domain.hh"
#include "General/Environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

DEFINE_LOG(ja, "Jiles")

Jiles::Jiles() : Model(),
numElems_{0}, MaxE_{0},
Ps_{0}, a_{0}, alpha_{0}, k_{0},c_{0},
mp_{nullptr},
timeStep_{0}, globalIter_{0},
isMH_{false}
{}

Jiles::~Jiles() {
}

void Jiles::Init(std::map<std::string, double> ParameterMap, UInt numElems,  UInt dim) {

  numElems_ = numElems;

  ElemNum2Idx_.Resize(numElems);
  ElemNum2Idx_.Init(0);
  for(UInt i = 0; i < numElems;i++){
    ElemNum2Idx_[i]=i;
  }

  if (ParameterMap.size() < 5) {
    EXCEPTION("The model needs 5 or more parameter!");
  }
  Ps_ = ParameterMap["Ps"];
  alpha_ = ParameterMap["alpha"];
  a_ = ParameterMap["a"];
  k_ = ParameterMap["k"];
  c_ = ParameterMap["c"];

  isMH_ = ParameterMap["isMH"];
  if(isMH_ == 1.0){
    varHandle_="cacheResult";
  } else {
    varHandle_="step";
  }

  // JilesAtherton keeps shared per-timestep history (RampUp/saveValues operate on all elements)
  // transient JA can have wrong numerical results when in parallel, the PrismHyst_MH can be forced to serial
  // TODO: JA need a refactoring to be truly thread-save
  if(isMH_ != 1.0 && CFS_NUM_THREADS > 1)
    WARN("Material model 'JilesAthertonModel' is not thread-safe for transient analysis");
    
  Ps_ = ParameterMap["Ps"];
  alpha_ = ParameterMap["alpha"];
  a_ = ParameterMap["a"];
  k_ = ParameterMap["k"];
  c_ = ParameterMap["c"];

  MaxE_ = 0;

  E0_.Resize(numElems_);
  E0_.Init(0);

  E0it_.Resize(numElems_);
  E0it_.Init(0);

  E1_.Resize(numElems_);
  E1_.Init(0);

  Pi0_.Resize(numElems_);
  Pi0_.Init(0);

  Pi0it_.Resize(numElems_);
  Pi0it_.Init(0);

  Pi1_.Resize(numElems_);
  Pi1_.Init(0);

  Pa0_.Resize(numElems_);
  Pa0_.Init(0);

  Pa0it_.Resize(numElems_);
  Pa0it_.Init(0);

  Pa1_.Resize(numElems_);
  Pa1_.Init(0);

  P0_.Resize(numElems_);
  P0_.Init(0);

  P0it_.Resize(numElems_);
  P0it_.Init(0);

  P1_.Resize(numElems_);
  P1_.Init(0);


  isFirstTime_.Resize(numElems_);
  isFirstTime_.Init(1);

  timeStep_ = 1;

  mp_ = domain->GetMathParser();
  globalIter_ = 0;
}

Double Jiles::ComputeMaterialParameter(Vector<Double> EVec, const Integer ElemNum) {

  // JilesAtherton is a serial model: the per-iteration / per-timestep bookkeeping
  // (globalIter_ swap, saveValues, first-time RampUp init) mutates shared,
  // all-element state, and the per-element math reads/writes the same shared
  // vectors. It also has lazy init (within first run) - TODO: refactor it
  // 
  // With the critical below we fix double free and the PrismHist_MH runs, other
  // JA tests can produce wrong numerical results when assembled in parallel.
  const UInt idx = ElemNum2Idx_[ElemNum-1];
  const Double E = EVec.NormL2();

  Double epsilon;
  const Double epsilon0 = 8.854187e-12;

  #pragma omp critical (jiles_transition)
  {
    if(globalIter_ != mp_->GetExprVars(MathParser::GLOB_HANDLER, "iterationCounter")){
      globalIter_ = mp_->GetExprVars(MathParser::GLOB_HANDLER, "iterationCounter");
      //if there is a new iteration, save the values from the previous iteration
      LOG_DBG3(ja) << "Trigger new iteration"<< std::endl;
      for(UInt i = 0; i <numElems_; i++){
        LOG_DBG3(ja) << "Overwritting for idx_: " << i << std::endl;
        E0it_[i]=E1_[i];
        P0it_[i]=P1_[i];
        Pi0it_[i]=Pi1_[i];
        Pa0it_[i]=Pa1_[i];
      }
    }

    // advance the timestep (global E0_=E1_ swap + reset inside saveValues)
    saveValues(false);

    if(E == 0){
      epsilon = 4.028353e-07;
    } else {
      if (isFirstTime_[idx]) {

        //Evaluate the anhysteretic curve to give a better starting point for P
        P0_[idx] = Ps_ * (cosh(E / a_) / sinh(E / a_) - a_ / E);

        //Ramp up, in case signal doesnt start at 0, to provide stable JA-solution.
        RampUp(512, E, idx);

        // Set first time for element to false. Each element initializes itself exactly once
        // note that with parallel assembly the order of element evaluation is not deterministic
        isFirstTime_[idx] = 0;
      }

      //if timestep == 0 -> new iteration, reset the values, only for multiharmonic
      if(timeStep_==0){
        E0_[idx]=E0it_[idx];
        P0_[idx]=P0it_[idx];
        Pi0_[idx]=Pi0it_[idx];
        Pa0_[idx]=Pa0it_[idx];
      }

      Double P = Evaluate(E, idx);
      epsilon = epsilon0 + P / E;
    }
  }

  if(std::isinf(epsilon) || std::isnan(epsilon) ){
    std::cout << "E0: "<< E0_[idx] << std::endl;
    std::cout << "E1: "<< E1_[idx]<< std::endl;
    std::cout << "dE: " <<  E1_[idx] - E0_[idx] << std::endl;
    std::cout << "P0: "<< P0_[idx] << std::endl;
    std::cout << "P1: "<< P1_[idx]<< std::endl;
    std::cout << "Pa0: "<< Pa0_[idx] << std::endl;
    std::cout << "Pa1: "<< Pa1_[idx]<< std::endl;
    std::cout << "Pi0: "<< Pi0_[idx] << std::endl;
    std::cout << "Pi1: "<< Pi1_[idx]<< std::endl;

    EXCEPTION("Epsilon is inifite or NaN...")
  }

//  if(idx_==0 && (timeStep_==0 || timeStep_== 127)){
//    LOG_DBG3(ja) << "Iteration: " << globalIter_ << std::endl;
//    LOG_DBG3(ja) << "Timestep = " << timeStep_ << std::endl;
//    LOG_DBG3(ja) << "For element with intern idx_: " << idx_ << std::endl;
//    LOG_DBG3(ja) << "Last Polarization (prev. It)      = " << P0_[idx_] << std::endl;
//    LOG_DBG3(ja)<< "Current Polarization (curr. It)   = " << P1_[idx_] << std::endl;
//    LOG_DBG3(ja) << "Last elecField (prev. It)         = " << E0_[idx_] << std::endl;
//    LOG_DBG3(ja) << "Current elecField (curr. It)      = " << E1_[idx_] << std::endl;
//    LOG_DBG3(ja) << "------------------------------------" << std::endl;
//  }
  return epsilon;
}

void Jiles::RampUp(Integer Nt, Double E, Integer idx) {

  LOG_DBG(ja) << "Computing Ramp Up-Polarization for element " << idx << std::endl ;

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

  //varHandle_ is different for transient/mh analysis
  if((timeStep_ != mp_->GetExprVars(MathParser::GLOB_HANDLER, varHandle_)) || InstantSave){
    Pa0_=Pa1_;
    Pi0_=Pi1_;
    P0_=P1_;
    E0_=E1_;

    Pa1_.Init(0);
    Pi1_.Init(0);
    P1_.Init(0);
    E1_.Init(0);

    timeStep_ = mp_->GetExprVars(MathParser::GLOB_HANDLER, varHandle_);
  }
}
}


