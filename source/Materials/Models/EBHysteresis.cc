#include <iostream>
#include <fstream>
#include <cmath>

#include "DataInOut/Logging/LogConfigurator.hh"

#include "EBHysteresis.hh"
#include "Model.hh"

#include "MatVec/Vector.hh"

#include "Utils/mathParser/mathParser.hh"
#include "Domain/Domain.hh"

namespace CoupledField {

DEFINE_LOG(eb, "EBHysteresis")

    EBHysteresis::EBHysteresis() : Model(),
    numElems_{0}, MaxE_{0},  idx_{0},
    Ps_{0}, A_{0}, mu0_{0}, numS_{0},chi_factor_{0},
    mp_{nullptr}, isFirstTimeFinished_{0},
    timeStep_{0}, globalIter_{0},
    isMH_{false}
    {}

    EBHysteresis::~EBHysteresis() {
    }

    void EBHysteresis::Init(std::map<std::string, double> ParameterMap, UInt numElems) {

      numElems_ = numElems;

      ElemNum2Idx_.Resize(numElems);
      ElemNum2Idx_.Init(0);
      for(UInt i = 0; i < numElems;i++){
        ElemNum2Idx_[i]=i;
      }

      if (ParameterMap.size() < 5) {
        EXCEPTION("The model needs 5 parameters!");
      }
      Ps_ = ParameterMap["Ps"];
      A_ = ParameterMap["A"];
      mu0_ = ParameterMap["mu0"];
      numS_ = ParameterMap["numS"];
      chi_factor_ = ParameterMap["chi_factor"];

      isMH_ = ParameterMap["isMH"];
      if(isMH_ == 1.0){
        varHandle_="cacheResult";
      } else {
        varHandle_="step";
      }
      epsilon0_ = 1.0/(4 * M_PI * 1e-7 * 299792458);

      Ps_ = ParameterMap["Ps"];
      A_ = ParameterMap["A"];
      mu0_ = ParameterMap["mu0"];
      numS_ = ParameterMap["numS"];
      chi_factor_ = ParameterMap["chi_factor"];

      E0_.Resize(numElems_, StdVector<Double>(2));
      E1_.Resize(numElems_, StdVector<Double>(2));
      P0_.Resize(numElems_, StdVector<Double>(2));
      P1_.Resize(numElems_, StdVector<Double>(2));

      HxS_prev_.Resize(numElems_, StdVector<Double>(numS_));
      HyS_prev_.Resize(numElems_, StdVector<Double>(numS_));
      MxS_prev_.Resize(numElems_, StdVector<Double>(numS_));
      MyS_prev_.Resize(numElems_, StdVector<Double>(numS_));

      HxS_n_.Resize(numElems_, StdVector<Double>(numS_));
      HyS_n_.Resize(numElems_, StdVector<Double>(numS_));
      MxS_n_.Resize(numElems_, StdVector<Double>(numS_));
      MyS_n_.Resize(numElems_, StdVector<Double>(numS_));

      HxS_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));
      HyS_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));
      MxS_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));
      MyS_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));

      epsilon_.Resize(numElems_, Matrix<Double>(2,2));

      isFirstTime_.Resize(numElems_);
      isFirstTime_.Init(1);

      isFirstTimeFinished_ = false;

      timeStep_ = 1;

      mp_ = domain->GetMathParser();
      globalIter_ = 0;
    }


    Double EBHysteresis::ComputeMaterialParameter(Vector<Double> EVec, const Integer ElemNum) {
      EXCEPTION("Not implemented here, this is a vector hysteresis model!");
    }



    Matrix<Double> EBHysteresis::ComputeTensorialMaterialParameter(Vector<Double> EVec, const Integer ElemNum) {

      if(globalIter_ != mp_->GetExprVars(MathParser::GLOB_HANDLER, "iterationCounter")){
        globalIter_ = mp_->GetExprVars(MathParser::GLOB_HANDLER, "iterationCounter");
        //if there is a new iteration, save the values from the previous iteration
        LOG_DBG3(eb) << "Trigger new iteration"<< std::endl;
      }

      idx_=ElemNum2Idx_[ElemNum-1];


      saveValues(false);
      Matrix<Double> epsilon(2,2);
      epsilon.InitValue(0.0);
      if(timeStep_== 1){
        std::cout<< "EVec = "<<EVec.ToString()<<std::endl;
        epsilon[0][0] = 40 * M_PI; //epsilon0_;
        epsilon[0][1] = 0.0;
        epsilon[1][0] = 0.0;
        epsilon[1][1] = 40 * M_PI; //epsilon0_;
        //E1_[idx_][0] = EVec[0];
        //E1_[idx_][1] = EVec[1];
        //E1_[idx_][0] = 0.0;
        //E1_[idx_][1] = 0.0;
        //P1_[idx_][0] = 0.0;
        //P1_[idx_][1] = 0.0;
        E0_[idx_][0] = 0.0;
        E0_[idx_][1] = 0.0;
        P0_[idx_][0] = 0.0;
        P0_[idx_][1] = 0.0;
        epsilon_[idx_] = epsilon;
        if((EVec[0] == 0) && (EVec[1] == 0)){
          return epsilon;
        }
      }



      Vector<Double> P, E;
      P = Evaluate(EVec);



      // get the values for E and P from the last timestep
      StdVector<Double>& E0 = E0_[idx_];
      StdVector<Double>& P0 = P0_[idx_];
      StdVector<Double> D_k(2);
      StdVector<Double> D_k_0(2);


      StdVector<Double> delta_E(2);
      StdVector<Double> delta_D(2);
      E1_[idx_].Resize(2, 0.0);
      P1_[idx_].Resize(2, 0.0);
      for(int i = 0; i < 2; i++){
        D_k[i] = epsilon0_ * EVec[i] + P[i];
        D_k_0[i] = epsilon0_ * E0[i] + P0[i];
        delta_E[i] = EVec[i] - E0[i];
        delta_D[i] = D_k[i] - D_k_0[i];
        E1_[idx_][i] = EVec[i];
        P1_[idx_][i] = P[i];
      }

      // std::cout << "EVec = " << EVec.ToString() << std::endl;
      // std::cout << "E0_[idx_] = " << E0_[idx_].ToString() << std::endl;
      // std::cout << "P0_[idx_] = " << P0_[idx_].ToString() << std::endl;
      // std::cout << "E1_[idx_] = " << E1_[idx_].ToString() << std::endl;
      // std::cout << "P1_[idx_] = " << P1_[idx_].ToString() << std::endl;
      // std::cout << "D_k = " << D_k.ToString() << std::endl;
      // std::cout << "D_k_0 = " << D_k_0.ToString() << std::endl;
      // std::cout << "delta_E = " << delta_E.ToString() << std::endl;
      // std::cout << "delta_D = " << delta_D.ToString() << std::endl;

      //epsilon = EvaluateLocalEpsilon(delta_E, delta_D);
      epsilon = EvaluateLocalEpsilonDFP(delta_E, delta_D);

      Double eps = 1e-8;
      // if( (std::abs(epsilon[0][0]) < eps) && (std::abs(epsilon[0][1]) < eps)
      //     && (std::abs(epsilon[1][0]) < eps) && (std::abs(epsilon[1][1]) < eps) ){
      //   epsilon[0][0] = epsilon0_;
      //   epsilon[0][1] = 0.0;
      //   epsilon[1][0] = 0.0;
      //   epsilon[1][1] = epsilon0_;
      // }
      
      //std::cout << "epsilon = " << epsilon.ToString() << std::endl;

      //P0_[idx_] = P1_[idx_];
      //E0_[idx_] = E1_[idx_];

      epsilon_[idx_] = epsilon;
      return epsilon;
    }

    Matrix<Double> EBHysteresis::EvaluateLocalEpsilon(StdVector<Double> dE, StdVector<Double> dD){
      Matrix<Double> epsilon;
      epsilon.Resize(2,2);

      // Double eps = 1e-14;
      // if( std::abs(dE[1]) < eps){
      //   epsilon[0][1] = 0.0;
      //   epsilon[1][1] = 0.0;
      // }else{
        epsilon[0][1] = 0.0; //dD[0]/dE[1];
        epsilon[1][1] = std::abs(dD[1])/std::abs(dE[1]);
      // }
      // if( std::abs(dE[0]) < eps){
      //   epsilon[1][0] = 0.0;
      //   epsilon[0][0] = 0.0;
      // }else{
        epsilon[1][0] = 0.0; //dD[1]/dE[0];
        epsilon[0][0] = std::abs(dD[0])/std::abs(dE[0]);
      //}

      return epsilon;
    }

    Matrix<Double> EBHysteresis::EvaluateLocalEpsilonDFP(StdVector<Double> dE, StdVector<Double> dD){
      Matrix<Double> epsilon;
      int dim = dE.size();
      epsilon.Resize(dim,dim);

      Matrix<Double> dD_dyadic_dEtransposed(dim, dim);
      Matrix<Double> I(dim, dim);
      for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
          dD_dyadic_dEtransposed[i][j] = dD[i] * dE[j];
        }
        I[i][i] = 1.0;
      }

      Matrix<Double> dE_dyadic_dDtransposed(dim, dim);
      for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
          dE_dyadic_dDtransposed[i][j] = dE[i] * dD[j];
        }
      }

      Matrix<Double> dD_dyadic_dDtransposed(dim, dim);
      for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
          dD_dyadic_dDtransposed[i][j] = dD[i] * dD[j];
        }
      }

      Vector<Double> deltaE(dim, 0.0);
      Vector<Double> deltaD(dim, 0.0);
      for(int i = 0; i < dim; i++){
        deltaE[i] = dE[i];
        deltaD[i] = dD[i];
      }
      Double dDtransposed_in_dE = deltaE.Inner(deltaD);
      
      Matrix<Double> leftM(dim, dim);
      //leftM = (I - dD_dyadic_dEtransposed / dDtransposed_in_dE);
      
      Matrix<Double> leftM_times_epsilon_km1(dim, dim);
      //leftM_times_epsilon_km1 = leftM * epsilon_[idx_];
      
      Matrix<Double> rightM(dim, dim);
      //rightM = (I - dE_dyadic_dDtransposed / dDtransposed_in_dE);
      
      Matrix<Double> A;
      //A = leftM_times_epsilon_km1 * rightM;

      //epsilon = A + dD_dyadic_dDtransposed / dDtransposed_in_dE;

      // std::cout << "dD_dyadic_dEtransposed = " << dD_dyadic_dEtransposed.ToString() << std::endl;
      // std::cout << "dE_dyadic_dDtransposed = " << dE_dyadic_dDtransposed.ToString() << std::endl;
      // std::cout << "dD_dyadic_dDtransposed = " << dD_dyadic_dDtransposed.ToString() << std::endl;
      // std::cout << "leftM = " << leftM.ToString() << std::endl;
      // std::cout << "leftM_times_epsilon_km1 = " << leftM_times_epsilon_km1.ToString() << std::endl;
      // std::cout << "rightM = " << rightM.ToString() << std::endl;
      // std::cout << "A = " << A.ToString() << std::endl;

      for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
          if (std::isinf(epsilon[i][j]) || std::isnan(epsilon[i][j])){
            Matrix<Double> e = epsilon_[idx_];
            epsilon[i][j] = e[i][j];
          }
        }
      }

      return epsilon;
    }

    Vector<Double> EBHysteresis::Evaluate(Vector<Double> Hn) {

      StdVector<Double> weight(numS_);
      weight.Init(1.0/numS_);

      StdVector<Double> chi(numS_);
      for(int i = 0; i < numS_; i++){
        chi[i] = chi_factor_ * (i-1) / (numS_ - 1);
      }

      StdVector<Double> error, dir, HrxS_sol, HryS_sol, MxS_sol, MyS_sol;
      StdVector<UInt> numIter;
      error.Resize(numS_, 0.0);
      dir.Resize(numS_, 0.0);
      numIter.Resize(numS_, 0);
      HrxS_sol.Resize(numS_, 0.0);
      HryS_sol.Resize(numS_, 0.0);
      MxS_sol.Resize(numS_, 0.0);
      MyS_sol.Resize(numS_, 0.0);

      Double Hex_x = Hn[0];
      Double Hex_y = Hn[1];

      Double HrxS_prev, HryS_prev, MxSprev, MySprev, phi, err, ux, uy,
             uabs, ux1, uy1, ux2, uy2, Man, Man1, g1, g2, phiNew, HrS, Px, Py;
      UInt iter;
      StdVector<Double>& HxS_prev = HxS_prev_[idx_];
      StdVector<Double>& HyS_prev = HyS_prev_[idx_];
      StdVector<Double>& MxS_prev = MxS_prev_[idx_];
      StdVector<Double>& MyS_prev = MyS_prev_[idx_];


      for(int k = 0; k < numS_; k++){
        HrxS_prev = HxS_prev[k];
        HryS_prev = HyS_prev[k];
        MxSprev = MxS_prev[k];
        MySprev = MyS_prev[k];
        if(std::sqrt( std::pow(Hex_x - HrxS_prev, 2) + std::pow(Hex_y - HryS_prev, 2) ) < chi[k]){
          HrxS_sol[k] = HrxS_prev;
          HryS_sol[k] = HryS_prev;
        }else{
          phi = std::atan2( (HryS_prev - Hex_y), (HrxS_prev - Hex_x) ); // already in rad
          err = 1.0;
          iter = 0;
          while(err > 1.0e-2 && iter < 10){
            ux = Hex_x + chi[k] * std::cos(phi);
            uy = Hex_y + chi[k] * std::sin(phi);
            uabs = std::sqrt( std::pow(ux, 2) + std::pow(uy, 2));
            ux1 = -chi[k] * std::sin(phi);
            uy1 = chi[k] * std::cos(phi);
            ux2 = -chi[k] * std::cos(phi);
            uy2 = -chi[k] * std::sin(phi);
            Man = 2*Ps_ / M_PI * std::atan(uabs/A_);
            Man1 = 2*Ps_ / (A_ * M_PI) * (1.0/(1.0+std::pow(uabs/A_, 2)));
            g1 = (Man/uabs) * (ux * ux1 + uy * uy1) - MxSprev*ux1 - MySprev*uy1;
            g2 = (uabs*Man1 - Man)/(std::pow(uabs,3)) * std::pow(ux*ux1 + uy*uy1, 2) + Man/uabs *
                (std::pow(ux1, 2) + std::pow(uy1, 2)) + Man/uabs * (ux*ux2 + uy*uy2) - MxSprev*ux2 - MySprev*uy2;
            phiNew = phi - g1/g2;
            err = std::abs(phiNew- phi);
            phi = phiNew;
            iter++;
          }
          error[k] = err;
          dir[k] = phi;
          numIter[k] = iter;
          HrxS_sol[k] = Hex_x +  chi[k]*std::cos(phi);
          HryS_sol[k] = Hex_y +  chi[k]*std::sin(phi);
        }
        HrS = std::sqrt( std::pow(HrxS_sol[k], 2) + std::pow(HryS_sol[k], 2));
        if( std::abs(HrS) > 1.0e-2){
          MxS_sol[k] = (2.0 * Ps_/M_PI) * std::atan(HrS/A_) * HrxS_sol[k]/HrS;
          MyS_sol[k] = (2.0 * Ps_/M_PI) * std::atan(HrS/A_) * HryS_sol[k]/HrS;
        }else{
          MxS_sol[k] = 0.0;
          MyS_sol[k] = 0.0;
        }
      }
      Px = 0.0;
      Py = 0.0;
      for(int k = 0; k < numS_; k++){
        Px += weight[k] * MxS_sol[k];
        Py += weight[k] * MyS_sol[k];
      }

      HxS_n_tmp_[idx_] = HrxS_sol;
      HyS_n_tmp_[idx_] = HryS_sol;
      MxS_n_tmp_[idx_] = MxS_sol;
      MyS_n_tmp_[idx_] = MyS_sol;

      Vector<Double> ret;
      ret.Push_back(Px);
      ret.Push_back(Py);
      return ret;
    }

    void EBHysteresis::saveValues(bool InstantSave){


      //varHandle_ is different for transient/mh analysis
      if((timeStep_ != mp_->GetExprVars(MathParser::GLOB_HANDLER, varHandle_)) || InstantSave){
        // std::cout << "======= BEFORE ================================"<<std::endl;
        // std::cout << "E0_ = " << E0_.ToString() << std::endl;
        // std::cout << "P0_ = " << P0_.ToString() << std::endl;
        // std::cout << "E1_ = " << E1_.ToString() << std::endl;
        // std::cout << "P1_ = " << P1_.ToString() << std::endl;
        // std::cout << "=======================================\n"<<std::endl;

        P0_ = P1_;
        E0_ = E1_;

        // std::cout << "======= AFTER ================================"<<std::endl;
        // std::cout << "E0_ = " << E0_.ToString() << std::endl;
        // std::cout << "P0_ = " << P0_.ToString() << std::endl;
        // std::cout << "E1_ = " << E1_.ToString() << std::endl;
        // std::cout << "P1_ = " << P1_.ToString() << std::endl;
        // std::cout << "=======================================\n"<<std::endl;



        HxS_n_ = HxS_n_tmp_;
        HyS_n_ = HyS_n_tmp_;
        MxS_n_ = MxS_n_tmp_;
        MyS_n_ = MyS_n_tmp_;

        //P1_.Clear(true);
        //E1_.Clear(true);

        timeStep_ = mp_->GetExprVars(MathParser::GLOB_HANDLER, varHandle_);
      }
    }
}


