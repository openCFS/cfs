#include <iostream>
#include <fstream>
#include <cmath>

#include "DataInOut/Logging/LogConfigurator.hh"

#include "EBHysteresis.hh"
#include "Model.hh"

#include "MatVec/Vector.hh"

#include "Utils/mathParser/mathParser.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"

namespace CoupledField {

DEFINE_LOG(eb, "EBHysteresis")

    EBHysteresis::EBHysteresis() : Model(),
    numElems_{0}, 
    Ps_{0}, A_{0}, mu0_{0}, numS_{0},chi_factor_{0},
    mp_{nullptr}, isFirstTimeFinished_{0},
    timeStep_{0}, globalIter_{0},
    isMH_{false}
    {}

    EBHysteresis::~EBHysteresis() {
    }

    void EBHysteresis::Init(std::map<std::string, double> ParameterMap, shared_ptr<ElemList> entityList, UInt dim) {

      dim_ = dim;

      EntityIterator it = entityList->GetIterator();
      numElems_ = entityList->GetSize();
      ElemNum2Idx_.clear();
      UInt num = 0;
      for (it.Begin(); !it.IsEnd(); it++) {
        UInt eNum = it.GetElem()->elemNum;
        ElemNum2Idx_[eNum] = num;
        num++;
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
      mu_0 = 4 * M_PI * 1e-7;

      Ps_ = ParameterMap["Ps"];
      A_ = ParameterMap["A"];
      mu0_ = ParameterMap["mu0"];
      numS_ = ParameterMap["numS"];
      chi_factor_ = ParameterMap["chi_factor"];

      H0_.Resize(numElems_, StdVector<Double>(dim_));
      H1_.Resize(numElems_, StdVector<Double>(dim_));
      M0_.Resize(numElems_, StdVector<Double>(dim_));
      M1_.Resize(numElems_, StdVector<Double>(dim_));

      HxS_prev_.Resize(numElems_, StdVector<Double>(numS_));
      HyS_prev_.Resize(numElems_, StdVector<Double>(numS_));
      HzS_prev_.Resize(numElems_, StdVector<Double>(numS_));
      MxS_prev_.Resize(numElems_, StdVector<Double>(numS_));
      MyS_prev_.Resize(numElems_, StdVector<Double>(numS_));
      MzS_prev_.Resize(numElems_, StdVector<Double>(numS_));

      HxS_n_.Resize(numElems_, StdVector<Double>(numS_));
      HyS_n_.Resize(numElems_, StdVector<Double>(numS_));
      HzS_n_.Resize(numElems_, StdVector<Double>(numS_));
      MxS_n_.Resize(numElems_, StdVector<Double>(numS_));
      MyS_n_.Resize(numElems_, StdVector<Double>(numS_));
      MzS_n_.Resize(numElems_, StdVector<Double>(numS_));

      HxS_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));
      HyS_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));
      HzS_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));
      MxS_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));
      MyS_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));
      MzS_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));

      mu_.Resize(numElems_, Matrix<Double>(dim_,dim_));

      isFirstTime_.Resize(numElems_);
      isFirstTime_.Init(1);

      isFirstTimeFinished_ = false;

      timeStep_ = 1;

      mp_ = domain->GetMathParser();
      globalIter_ = 0;

      hasElemSolution_.Resize(numElems_, false);
    }


    Double EBHysteresis::ComputeMaterialParameter(Vector<Double> HVec, const Integer ElemNum) {
      EXCEPTION("Not implemented here, this is a vector hysteresis model!");
    }


    Vector<Double> EBHysteresis::GetFluxDensity(Vector<Double> HVec, const Integer ElemNum) {
      Vector<Double> B(dim_);
      UInt idx = ElemNum2Idx_[ElemNum];

      LOG_DBG3(eb) << "\n\t HVec = " << HVec.ToString();
      Vector<Double> M;
      M = Evaluate(HVec, false, idx);
      LOG_DBG3(eb) << "\n\t M = " << M.ToString();
      
      for(UInt i = 0; i < dim_; ++i){
        B[i] = mu_0 * (HVec[i] + M[i]);
      }

      LOG_DBG3(eb) << "\n\t B = " << B.ToString();
      return B;
    }



    Matrix<Double> EBHysteresis::ComputeTensorialMaterialParameter(Vector<Double> HVec, const Integer ElemNum) {
#pragma omp critical
{
      if(globalIter_ != mp_->GetExprVars(MathParser::GLOB_HANDLER, "iterationCounter")){
        globalIter_ = mp_->GetExprVars(MathParser::GLOB_HANDLER, "iterationCounter");
        //if there is a new iteration, save the values from the previous iteration
        LOG_DBG3(eb) << "\n\t Trigger new iteration"<< std::endl;
        hasElemSolution_.Init(false);
      }
      if( (timeStep_ != mp_->GetExprVars(MathParser::GLOB_HANDLER, varHandle_))){
        hasElemSolution_.Init(false);
        HxS_n_ = HxS_n_tmp_;
        HyS_n_ = HyS_n_tmp_;
        HzS_n_ = HzS_n_tmp_;
        MxS_n_ = MxS_n_tmp_;
        MyS_n_ = MyS_n_tmp_;
        MzS_n_ = MzS_n_tmp_;
        timeStep_ = mp_->GetExprVars(MathParser::GLOB_HANDLER, varHandle_);
      }
}
      UInt idx = ElemNum2Idx_[ElemNum];

      // since the material tensor is computed per element and one element has 
      // several integration points, there are several evaluations of this
      // ComputeTensorialMaterialParameter function. But after every evaluation
      // of this function we store the H1_ and M1_ and mu_ variables for the particular
      // element and re-evaluating it again, would result in wrong delta values for the
      // Broyden or DFP formula. Therefore, we check if this particular element was already
      // evaluated in this iteration step. If so, we return the already computed values.
      if(hasElemSolution_[idx] == true){
        return mu_[idx];
      }

      Matrix<Double> mu(dim_,dim_);
      mu.InitValue(0.0);

      LOG_DBG3(eb) << "\n\t HVec = " << HVec.ToString();
      if(timeStep_== 1 && globalIter_ == 1){
        mu[0][0] = mu_0;
        mu[0][1] = 0.0;
        mu[1][0] = 0.0;
        mu[1][1] = mu_0;
        if(dim_ == 3){
          mu[2][0] = 0.0;
          mu[2][1] = 0.0;
          mu[0][2] = 0.0;
          mu[1][2] = 0.0;
          mu[2][2] = mu_0;
        }
        mu_[idx] = mu;
        return mu;
      }

      Vector<Double> M, H;
      M = Evaluate(HVec, true, idx);

      // get the values for H and M from the last timestep
      StdVector<Double>& H0 = H0_[idx];
      StdVector<Double>& M0 = M0_[idx];
      StdVector<Double> B_k(dim_);
      StdVector<Double> B_k_0(dim_);


      StdVector<Double> delta_H(dim_);
      StdVector<Double> delta_B(dim_);
      for(UInt i = 0; i < dim_; i++){
        B_k[i] = mu_0 * (HVec[i] + M[i]);
        B_k_0[i] = mu_0 * (H0[i] + M0[i]);
        delta_H[i] = HVec[i] - H0[i];
        delta_B[i] = B_k[i] - B_k_0[i];
        H1_[idx][i] = HVec[i];
        M1_[idx][i] = M[i];
      }



      LOG_DBG3(eb)<< "AFTER EVALUATION" 
                  << "\n\t HVec =  \t" << HVec.ToString()
                  << "\n\t H0_["<<idx<<"] = \t"<<H0_[idx].ToString()
                  << "\n\t M0_["<<idx<<"] = \t"<<M0_[idx].ToString()
                  << "\n\t H1_["<<idx<<"] = \t"<<H1_[idx].ToString()
                  << "\n\t M1_["<<idx<<"] = \t"<<M1_[idx].ToString()
                  << "\n\t B_k = \t\t"<<B_k.ToString()
                  << "\n\t B_k_0 = \t"<<B_k_0.ToString()
                  << "\n\t delta_H = \t"<<delta_H.ToString()
                  << "\n\t delta_B = \t"<<delta_B.ToString();

      //mu = EvaluateLocalMu(delta_H, delta_B);
      //mu = EvaluateLocalMuDFP(delta_H, delta_B);
      mu = EvaluateLocalMuGBM(delta_H, delta_B, idx);


      LOG_DBG3(eb) << "\n\t mu = " << mu.ToString();

      mu_[idx] = mu;
      M0_[idx] = M1_[idx];
      H0_[idx] = H1_[idx];

      // mark this element as computed
      hasElemSolution_[idx] = true;
      return mu;
    }

    Matrix<Double> EBHysteresis::EvaluateLocalMu(StdVector<Double> dH, StdVector<Double> dB, UInt idx){
      Matrix<Double> mu;
      mu.Resize(dim_,dim_);
      mu[0][1] = 0.0; //dB[0]/dH[1];
      mu[1][1] = std::abs(dB[1])/std::abs(dH[1]);
      mu[1][0] = 0.0; //dB[1]/dH[0];
      mu[0][0] = std::abs(dB[0])/std::abs(dH[0]);
      if(dim_ == 3){
        mu[2][0] = 0.0;
        mu[2][1] = 0.0;
        mu[0][2] = 0.0;
        mu[1][2] = 0.0;
        mu[2][2] = std::abs(dB[2])/std::abs(dH[2]);
      }

      return mu;
    }

    Matrix<Double> EBHysteresis::EvaluateLocalMuDFP(StdVector<Double> dH, StdVector<Double> dB, UInt idx){
      Matrix<Double> mu;
      mu.Resize(dim_,dim_);

      Matrix<Double> dD_dyadic_dEtransposed(dim_, dim_);
      Matrix<Double> dE_dyadic_dDtransposed(dim_, dim_);
      Matrix<Double> dD_dyadic_dDtransposed(dim_, dim_);
      Matrix<Double> I(dim_, dim_);
      Matrix<Double> leftM(dim_, dim_);
      Matrix<Double> rightM(dim_, dim_);
      Matrix<Double> A(dim_,dim_);
      
      if(dim_ == 2){
        I[0][0] = 1.0;
        I[1][1] = 1.0;

        dD_dyadic_dEtransposed[0][0] = dB[0] * dH[0];
        dD_dyadic_dEtransposed[1][0] = dB[1] * dH[0];
        dD_dyadic_dEtransposed[0][1] = dB[0] * dH[1];
        dD_dyadic_dEtransposed[1][1] = dB[1] * dH[1];

        dE_dyadic_dDtransposed[0][0] = dH[0] * dB[0];
        dE_dyadic_dDtransposed[1][0] = dH[1] * dB[0];
        dE_dyadic_dDtransposed[0][1] = dH[0] * dB[1];
        dE_dyadic_dDtransposed[1][1] = dH[1] * dB[1];

        dD_dyadic_dDtransposed[0][0] = dB[0] * dB[0];
        dD_dyadic_dDtransposed[1][0] = dB[1] * dB[0];
        dD_dyadic_dDtransposed[0][1] = dB[0] * dB[1];
        dD_dyadic_dDtransposed[1][1] = dB[1] * dB[1];

        Double dDtransposed_in_dE = dH[0]*dB[0] + dH[1]*dB[1];
        

        leftM[0][0] = 1.0 - dD_dyadic_dEtransposed[0][0] / dDtransposed_in_dE;
        leftM[1][0] = dD_dyadic_dEtransposed[1][0] / dDtransposed_in_dE;
        leftM[0][1] = dD_dyadic_dEtransposed[0][1] / dDtransposed_in_dE;
        leftM[1][1] = 1.0 - dD_dyadic_dEtransposed[1][1] / dDtransposed_in_dE;
        
        Matrix<Double> leftM_times_epsilon_km1(leftM);
        leftM.MultT(mu_[idx], leftM_times_epsilon_km1);
        
        rightM[0][0] = 1.0 - dE_dyadic_dDtransposed[0][0] / dDtransposed_in_dE;
        rightM[1][0] = dE_dyadic_dDtransposed[1][0] / dDtransposed_in_dE;
        rightM[0][1] = dE_dyadic_dDtransposed[0][1] / dDtransposed_in_dE;
        rightM[1][1] = 1.0 - dE_dyadic_dDtransposed[1][1] / dDtransposed_in_dE;
      
        leftM_times_epsilon_km1.MultT(rightM, A);

        mu[0][0] = A[0][0] + dD_dyadic_dDtransposed[0][0]/dDtransposed_in_dE;
        mu[1][0] = A[1][0] + dD_dyadic_dDtransposed[1][0]/dDtransposed_in_dE;
        mu[0][1] = A[0][1] + dD_dyadic_dDtransposed[0][1]/dDtransposed_in_dE;
        mu[1][1] = A[1][1] + dD_dyadic_dDtransposed[1][1]/dDtransposed_in_dE;


        LOG_DBG3(eb)<< "\n\t dD_dyadic_dEtransposed = " << dD_dyadic_dEtransposed.ToString()
                    << "\n\t dE_dyadic_dDtransposed = " << dE_dyadic_dDtransposed.ToString()
                    << "\n\t dD_dyadic_dDtransposed = " << dD_dyadic_dDtransposed.ToString()
                    << "\n\t dDtransposed_in_dE = " << dDtransposed_in_dE
                    << "\n\t leftM = " << leftM.ToString()
                    << "\n\t leftM_times_epsilon_km1 = " << leftM_times_epsilon_km1.ToString()
                    << "\n\t rightM = " << rightM.ToString()
                    << "\n\t epsilon_km1 = " << mu_[idx].ToString()
                    << "\n\t A = " << A.ToString();

        mu[0][1] = 0.0;
        mu[1][0] = 0.0;
        for (UInt i = 0; i < dim_; i++){
            if (std::isinf(mu[i][i]) || std::isnan(mu[i][i])){
              Matrix<Double> e = mu_[idx];
                mu[i][i] = mu_0; //e[i][i]; //mu_0;
            } 
        }
      }else{
        
        // 3D version
        I[0][0] = 1.0;
        I[1][1] = 1.0;
        I[2][2] = 1.0;


        dD_dyadic_dEtransposed[0][1] = dB[0] * dH[1];
        dD_dyadic_dEtransposed[1][0] = dB[1] * dH[0];
        dD_dyadic_dEtransposed[0][2] = dB[0] * dH[2];
        dD_dyadic_dEtransposed[2][0] = dB[2] * dH[0];
        dD_dyadic_dEtransposed[1][2] = dB[1] * dH[2];
        dD_dyadic_dEtransposed[2][1] = dB[2] * dH[1];
        dD_dyadic_dEtransposed[0][0] = dB[0] * dH[0];        
        dD_dyadic_dEtransposed[1][1] = dB[1] * dH[1];
        dD_dyadic_dEtransposed[2][2] = dB[2] * dH[2];

        dE_dyadic_dDtransposed[1][0] = dH[1] * dB[0];
        dE_dyadic_dDtransposed[0][1] = dH[0] * dB[1];
        dE_dyadic_dDtransposed[0][2] = dH[0] * dB[2];
        dE_dyadic_dDtransposed[2][0] = dH[2] * dB[0];
        dE_dyadic_dDtransposed[1][2] = dH[1] * dB[2];
        dE_dyadic_dDtransposed[2][1] = dH[2] * dB[1];
        dE_dyadic_dDtransposed[0][0] = dH[0] * dB[0];
        dE_dyadic_dDtransposed[1][1] = dH[1] * dB[1];
        dE_dyadic_dDtransposed[2][2] = dH[2] * dB[2];

        dD_dyadic_dDtransposed[1][0] = dB[1] * dB[0];
        dD_dyadic_dDtransposed[0][1] = dB[0] * dB[1];        
        dD_dyadic_dDtransposed[0][2] = dB[0] * dB[2];
        dD_dyadic_dDtransposed[2][0] = dB[2] * dB[0];
        dD_dyadic_dDtransposed[1][2] = dB[1] * dB[2];
        dD_dyadic_dDtransposed[2][1] = dB[2] * dB[1];
        dD_dyadic_dDtransposed[0][0] = dB[0] * dB[0];
        dD_dyadic_dDtransposed[1][1] = dB[1] * dB[1];
        dD_dyadic_dDtransposed[2][2] = dB[2] * dB[2];


        Double dDtransposed_in_dE = dH[0]*dB[0] + dH[1]*dB[1] + dH[2]*dB[2];

        leftM[1][0] = dD_dyadic_dEtransposed[1][0] / dDtransposed_in_dE;
        leftM[0][1] = dD_dyadic_dEtransposed[0][1] / dDtransposed_in_dE;
        leftM[0][2] = dD_dyadic_dEtransposed[0][2] / dDtransposed_in_dE;
        leftM[2][0] = dD_dyadic_dEtransposed[2][0] / dDtransposed_in_dE;
        leftM[1][2] = dD_dyadic_dEtransposed[1][2] / dDtransposed_in_dE;
        leftM[2][1] = dD_dyadic_dEtransposed[2][1] / dDtransposed_in_dE;
        leftM[0][0] = 1.0 - dD_dyadic_dEtransposed[0][0] / dDtransposed_in_dE;
        leftM[1][1] = 1.0 - dD_dyadic_dEtransposed[1][1] / dDtransposed_in_dE;
        leftM[2][2] = 1.0 - dD_dyadic_dEtransposed[2][2] / dDtransposed_in_dE;

        Matrix<Double> leftM_times_epsilon_km1(leftM);
        leftM.MultT(mu_[idx], leftM_times_epsilon_km1);

        rightM[1][0] = dE_dyadic_dDtransposed[1][0] / dDtransposed_in_dE;
        rightM[0][1] = dE_dyadic_dDtransposed[0][1] / dDtransposed_in_dE;
        rightM[0][2] = dE_dyadic_dDtransposed[0][2] / dDtransposed_in_dE;
        rightM[2][0] = dE_dyadic_dDtransposed[2][0] / dDtransposed_in_dE;
        rightM[1][2] = dE_dyadic_dDtransposed[1][2] / dDtransposed_in_dE;
        rightM[2][1] = dE_dyadic_dDtransposed[2][1] / dDtransposed_in_dE;
        rightM[0][0] = 1.0 - dE_dyadic_dDtransposed[0][0] / dDtransposed_in_dE;
        rightM[1][1] = 1.0 - dE_dyadic_dDtransposed[1][1] / dDtransposed_in_dE;
        rightM[2][2] = 1.0 - dE_dyadic_dDtransposed[2][2] / dDtransposed_in_dE;

        leftM_times_epsilon_km1.MultT(rightM, A);

        mu[1][0] = A[1][0] + dD_dyadic_dDtransposed[1][0]/dDtransposed_in_dE;
        mu[0][1] = A[0][1] + dD_dyadic_dDtransposed[0][1]/dDtransposed_in_dE;
        mu[0][2] = A[0][2] + dD_dyadic_dDtransposed[0][2]/dDtransposed_in_dE;
        mu[2][0] = A[2][0] + dD_dyadic_dDtransposed[2][0]/dDtransposed_in_dE;
        mu[1][2] = A[1][2] + dD_dyadic_dDtransposed[1][2]/dDtransposed_in_dE;
        mu[2][1] = A[2][1] + dD_dyadic_dDtransposed[2][1]/dDtransposed_in_dE;
        mu[0][0] = A[0][0] + dD_dyadic_dDtransposed[0][0]/dDtransposed_in_dE;
        mu[1][1] = A[1][1] + dD_dyadic_dDtransposed[1][1]/dDtransposed_in_dE;
        mu[2][2] = A[2][2] + dD_dyadic_dDtransposed[2][2]/dDtransposed_in_dE;

        LOG_DBG3(eb)<< "\n\t dD_dyadic_dEtransposed = " << dD_dyadic_dEtransposed.ToString()
                    << "\n\t dE_dyadic_dDtransposed = " << dE_dyadic_dDtransposed.ToString()
                    << "\n\t dD_dyadic_dDtransposed = " << dD_dyadic_dDtransposed.ToString()
                    << "\n\t dDtransposed_in_dE = " << dDtransposed_in_dE
                    << "\n\t leftM = " << leftM.ToString()
                    << "\n\t leftM_times_epsilon_km1 = " << leftM_times_epsilon_km1.ToString()
                    << "\n\t rightM = " << rightM.ToString()
                    << "\n\t epsilon_km1 = " << mu_[idx].ToString()
                    << "\n\t A = " << A.ToString();

        mu[0][1] = 0.0;
        mu[1][0] = 0.0;
        mu[2][0] = 0.0;
        mu[0][2] = 0.0;
        mu[1][2] = 0.0;
        mu[2][1] = 0.0;
        for (UInt i = 0; i < dim_; i++){
            if (std::isinf(mu[i][i]) || std::isnan(mu[i][i])){
              Matrix<Double> e = mu_[idx];
                mu[i][i] = e[i][i];//mu_0; //e[i][i]; //mu_0;
            } 
        }
      }
      return mu;
    }



    Matrix<Double> EBHysteresis::EvaluateLocalMuGBM(StdVector<Double> dH, StdVector<Double> dB, UInt idx){
      Matrix<Double> mu;
      mu.Resize(dim_,dim_);

      Vector<Double> dD_minus_epsilon_km1_times_dE(dim_);
      Vector<Double> dHvec(dim_);
      Matrix<Double> rightM(dim_,dim_);

      if(dim_ == 2){
        dD_minus_epsilon_km1_times_dE[0] = dB[0] - (mu_[idx][0][0]*dH[0] + mu_[idx][0][1]*dH[1]);
        dD_minus_epsilon_km1_times_dE[1] = dB[1] - (mu_[idx][1][0]*dH[0] + mu_[idx][1][1]*dH[1]);

        dHvec[0] = dH[0];
        dHvec[1] = dH[1];
        dD_minus_epsilon_km1_times_dE.ScalarDiv(dHvec.NormL2_squared());
        
        rightM.DyadicMult(dD_minus_epsilon_km1_times_dE, dHvec);
        
        mu[0][0] = mu_[idx][0][0] + rightM[0][0];
        mu[1][0] = mu_[idx][1][0] + rightM[1][0];
        mu[0][1] = mu_[idx][0][1] + rightM[0][1];
        mu[1][1] = mu_[idx][1][1] + rightM[1][1];
                
        LOG_DBG3(eb)<< "\n\t dD_minus_epsilon_km1_times_dE = " << dD_minus_epsilon_km1_times_dE.ToString()
                    << "\n\t dHvec = " << dHvec.ToString()
                    << "\n\t rightM = " << rightM.ToString();

        mu[0][1] = 0.0;
        mu[1][0] = 0.0;
        for (UInt i = 0; i < dim_; i++){
            if (std::isinf(mu[i][i]) || std::isnan(mu[i][i])){
              Matrix<Double> e = mu_[idx];
                mu[i][i] = e[i][i]; //mu_0; //e[i][i]; //mu_0;
            } 
        }
      }else{
        
        // 3D version
        dD_minus_epsilon_km1_times_dE[0] = dB[0] - (mu_[idx][0][0]*dH[0] + mu_[idx][0][1]*dH[1] + mu_[idx][0][2]*dH[2]);
        dD_minus_epsilon_km1_times_dE[1] = dB[1] - (mu_[idx][1][0]*dH[0] + mu_[idx][1][1]*dH[1] + mu_[idx][1][2]*dH[2]);
        dD_minus_epsilon_km1_times_dE[2] = dB[2] - (mu_[idx][2][0]*dH[0] + mu_[idx][2][1]*dH[1] + mu_[idx][2][2]*dH[2]);

        dHvec[0] = dH[0];
        dHvec[1] = dH[1];
        dHvec[2] = dH[2];
        dD_minus_epsilon_km1_times_dE.ScalarDiv(dHvec.NormL2_squared());
        
        rightM.DyadicMult(dD_minus_epsilon_km1_times_dE, dHvec);
        
        mu[1][0] = mu_[idx][1][0] + rightM[1][0];
        mu[0][1] = mu_[idx][0][1] + rightM[0][1];
        mu[0][2] = mu_[idx][0][2] + rightM[0][2];
        mu[2][0] = mu_[idx][2][0] + rightM[2][0];
        mu[1][2] = mu_[idx][1][2] + rightM[1][2];
        mu[2][1] = mu_[idx][2][1] + rightM[2][1];
        mu[0][0] = mu_[idx][0][0] + rightM[0][0];
        mu[1][1] = mu_[idx][1][1] + rightM[1][1];
        mu[2][2] = mu_[idx][2][2] + rightM[2][2];

        LOG_DBG3(eb)<< "\n\t dD_minus_epsilon_km1_times_dE = " << dD_minus_epsilon_km1_times_dE.ToString()
                    << "\n\t dHvec = " << dHvec.ToString()
                    << "\n\t rightM = " << rightM.ToString();

        mu[0][1] = 0.0;
        mu[1][0] = 0.0;
        mu[2][0] = 0.0;
        mu[0][2] = 0.0;
        mu[1][2] = 0.0;
        mu[2][1] = 0.0;
        for (UInt i = 0; i < dim_; i++){
            if (std::isinf(mu[i][i]) || std::isnan(mu[i][i])){
              Matrix<Double> e = mu_[idx];
                mu[i][i] = e[i][i]; //e[i][i]; //mu_0;
            } 
        }
      }
      return mu;
    }

    StdVector<Double> EBHysteresis::inv3x3(StdVector<Double> A){
      double detA = 0;
      StdVector<Double> adjA(9); adjA.Init(0.0);
      StdVector<Double> invA(9); invA.Init(0.0);
      detA = A[0]*(A[8]*A[4]-A[5]*A[7]) - (A[1]*(A[8]*A[3]-A[5]*A[6])) + A[2]*(A[7]*A[3]-A[4]*A[6]);
      adjA[0] = A[8]*A[4] - A[5]*A[7];
      adjA[1] = -(A[8]*A[3]-A[5]*A[6]);
      adjA[2] = A[7]*A[3]-A[4]*A[6];
      adjA[3] = -(A[8]*A[1]-A[2]*A[7]);
      adjA[4] = A[8]*A[0]-A[2]*A[6];
      adjA[5] = -(A[7]*A[0]-A[1]*A[6]);
      adjA[6] = A[5]*A[1]-A[2]*A[4];
      adjA[7] = -(A[5]*A[0]-A[2]*A[3]);
      adjA[8] = A[4]*A[0]-A[1]*A[3];
      invA[0] = (1/detA)*adjA[0];
      invA[1] = (1/detA)*adjA[1];
      invA[2] = (1/detA)*adjA[2];
      invA[3] = (1/detA)*adjA[3];
      invA[4] = (1/detA)*adjA[4];
      invA[5] = (1/detA)*adjA[5];
      invA[6] = (1/detA)*adjA[6];
      invA[7] = (1/detA)*adjA[7];
      invA[8] = (1/detA)*adjA[8];
      return invA;
    }


    Vector<Double> EBHysteresis::Evaluate(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx) {

      StdVector<Double> weight(numS_);
      weight.Init(1.0/numS_);
      Vector<Double> ret;

      StdVector<Double> chi(numS_);
      for(UInt i = 0; i < numS_; i++){
        chi[i] = chi_factor_ * i / (numS_ - 1) + 0.001;
      }

      if(dim_ == 2){
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
        StdVector<Double>& HxS_prev = HxS_n_[idx];
        StdVector<Double>& HyS_prev = HyS_n_[idx];
        StdVector<Double>& MxS_prev = MxS_n_[idx];
        StdVector<Double>& MyS_prev = MyS_n_[idx];


        for(int k = 0; k < numS_; k++){
          HrxS_prev = HxS_prev[k];
          HryS_prev = HyS_prev[k];
          MxSprev = MxS_prev[k];
          MySprev = MyS_prev[k];
          //if(std::sqrt( std::pow(Hex_x - HrxS_prev, 2) + std::pow(Hex_y - HryS_prev, 2) ) <= chi[k]){
          if(( std::pow((Hex_x - HrxS_prev)/chi[k], 2) + std::pow((Hex_y - HryS_prev)/chi[k], 2) ) <= 1.0){
            HrxS_sol[k] = HrxS_prev;
            HryS_sol[k] = HryS_prev;
          }else{
            StdVector<Double> i_phi_initial(2);
            i_phi_initial.Init(0.0);
            // dissipation now reached (just arrived at the border of the sphere)
            if(k == 0){
              i_phi_initial[0] = 0.0;
              i_phi_initial[1] = 0.0;
            }else{
              // use the direction of the vector play model as initial direction
              i_phi_initial[0] = (1.0/chi[k] * (HrxS_prev - Hex_x)) / (std::pow((Hex_x - HrxS_prev)/chi[k],2) + std::pow((Hex_y - HryS_prev)/chi[k],2));
              i_phi_initial[1] = (1.0/chi[k] * (HryS_prev - Hex_y)) / (std::pow((Hex_x - HrxS_prev)/chi[k],2) + std::pow((Hex_y - HryS_prev)/chi[k],2));
            }
            phi = std::atan2( i_phi_initial[1], i_phi_initial[0] ); // already in rad
            err = 1.0;
            iter = 0;
            while(err > 1.0e-3 && iter < 10){
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
              if(std::sqrt(std::pow(g2,2)) < 1e-8){
                phiNew = 0;
              }else{
                phiNew = phi - g1/g2;
              }
              err = std::sqrt(std::pow(phiNew-phi,2));
              phi = phiNew;
              iter++;
            }
            error[k] = err;
            dir[k] = phi;
            numIter[k] = iter;
            HrxS_sol[k] = Hex_x +  chi[k]*std::cos(phi);
            HryS_sol[k] = Hex_y +  chi[k]*std::sin(phi);
          }
          HrS = std::sqrt(std::pow(HrxS_sol[k], 2) + std::pow(HryS_sol[k], 2));
          if( std::sqrt(std::pow(HrS,2)) > 1.0e-12){
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
        if(saveTmpStageVecs){
          HxS_n_tmp_[idx] = HrxS_sol;
          HyS_n_tmp_[idx] = HryS_sol;
          MxS_n_tmp_[idx] = MxS_sol;
          MyS_n_tmp_[idx] = MyS_sol;
        }
        ret.Push_back(Px);
        ret.Push_back(Py);

      }else{   
        // 3D Version
        // volumetric weight for each pseudo particle
        StdVector<Double> weight(numS_);
        weight.Init(1.0/numS_);
        // dissipation scalar parameter (can also be a tensor if anisotropy is needed) 
        StdVector<Double> chi(numS_);
        Double offset = 1e-22;
        for(UInt i = 0; i < numS_; i++){
          chi[i] = chi_factor_ * i / (numS_ - 1) + offset;
        }
        StdVector<Double> error, dir_phi, dir_theta, HrxS_sol, HryS_sol,HrzS_sol, MxS_sol, MyS_sol, MzS_sol;
        StdVector<UInt> numIter;
        error.Resize(numS_, 0.0);
        dir_phi.Resize(numS_, 0.0);
        dir_theta.Resize(numS_, 0.0);
        numIter.Resize(numS_, 0);
        HrxS_sol.Resize(numS_, 0.0);
        HryS_sol.Resize(numS_, 0.0);
        HrzS_sol.Resize(numS_, 0.0);
        MxS_sol.Resize(numS_, 0.0);
        MyS_sol.Resize(numS_, 0.0);
        MzS_sol.Resize(numS_, 0.0);

        Double Hex_x = Hn[0];
        Double Hex_y = Hn[1];
        Double Hex_z = Hn[2];
        // all needed variables
        Double HrxS_prev, HryS_prev, HrzS_prev, MxSprev, MySprev, MzSprev, phi, err, ux, uy, uz, uabs,
            du_dphi_x, du_dphi_y, du_dphi_z, du_dtheta_x, du_dtheta_y, du_dtheta_z,
            d2u_dphi2_x, d2u_dphi2_y, d2u_dphi2_z, d2u_dtheta2_x, d2u_dtheta2_y, d2u_dtheta2_z,
            d2u_dthetadphi_x, d2u_dthetadphi_y, d2u_dthetadphi_z,
            Man, Man_strich, theta,
            g_strich_x, g_strich_y,
            phi_initial, theta_initial,
            g_strich_strich_11, g_strich_strich_12, g_strich_strich_21, g_strich_strich_22,
            thetaNew, phiNew, det_g_strich_strich, inv_det_g_strich_strich, normHrev, Px, Py, Pz, c1,c2,
            i_correct_x, i_correct_y, i_correct_z, Hirr_x, Hirr_y, Hirr_z,
            condition1, norm_i_initial;
        UInt iter;
        StdVector<Double>& HxS_prev = HxS_n_[idx];
        StdVector<Double>& HyS_prev = HyS_n_[idx];
        StdVector<Double>& HzS_prev = HzS_n_[idx];
        StdVector<Double>& MxS_prev = MxS_n_[idx];
        StdVector<Double>& MyS_prev = MyS_n_[idx];
        StdVector<Double>& MzS_prev = MzS_n_[idx];
        StdVector<Double> invCHI(9); invCHI.Init(0.0);
        StdVector<Double> CHI(9); CHI.Init(0.0);
        Px = 0.0;
        Py = 0.0;
        Pz = 0.0;
        // #############################################################
        //LOOP OVER ALL PARTICLES
        //#############################################################
        for(int k = 0; k < numS_; k++){
          HrxS_prev = HxS_prev[k];
          HryS_prev = HyS_prev[k];
          HrzS_prev = HzS_prev[k];
          MxSprev = MxS_prev[k];
          MySprev = MyS_prev[k];
          MzSprev = MzS_prev[k];
          // checks what the state of the dissipation is
          CHI[0] = chi[k]; CHI[1] = offset; CHI[2] = offset; CHI[3] = offset; CHI[4] = chi[k]; CHI[5] = offset; CHI[6] = offset; CHI[7] = offset; CHI[8] = chi[k];
          condition1 = (std::pow((Hex_x-HrxS_prev)/chi[k],2) + std::pow((Hex_y-HryS_prev)/chi[k],2) + std::pow((Hex_z-HrzS_prev)/chi[k],2));
          if( condition1 <= 1 ){
            HrxS_sol[k] = HrxS_prev;
            HryS_sol[k] = HryS_prev;
            HrzS_sol[k] = HrzS_prev;
          }else if( condition1 > 1.0){
            StdVector<Double> i_initial(3); i_initial.Init(0.0);
            invCHI = inv3x3(CHI);
            if(k > 0){
              // use the direction of the vector play model as initial direction
/*                 i_initial[0] = invCHI[0]*(HrxS_prev-Hex_x)+invCHI[1]*(HryS_prev-Hex_y)+invCHI[2]*(HrzS_prev-Hex_z) / std::sqrt(std::pow(invCHI[0]*(HrxS_prev-Hex_x)+invCHI[1]*(HryS_prev-Hex_y)+invCHI[2]*(HrzS_prev-Hex_z),2)+std::pow(invCHI[3]*(HrxS_prev-Hex_x)+invCHI[4]*(HryS_prev-Hex_y)+invCHI[5]*(HrzS_prev-Hex_z),2)+std::pow(invCHI[6]*(HrxS_prev-Hex_x)+invCHI[7]*(HryS_prev-Hex_y)+invCHI[8]*(HrzS_prev-Hex_z),2));
                i_initial[1] = invCHI[3]*(HrxS_prev-Hex_x)+invCHI[4]*(HryS_prev-Hex_y)+invCHI[5]*(HrzS_prev-Hex_z) / std::sqrt(std::pow(invCHI[0]*(HrxS_prev-Hex_x)+invCHI[1]*(HryS_prev-Hex_y)+invCHI[2]*(HrzS_prev-Hex_z),2)+std::pow(invCHI[3]*(HrxS_prev-Hex_x)+invCHI[4]*(HryS_prev-Hex_y)+invCHI[5]*(HrzS_prev-Hex_z),2)+std::pow(invCHI[6]*(HrxS_prev-Hex_x)+invCHI[7]*(HryS_prev-Hex_y)+invCHI[8]*(HrzS_prev-Hex_z),2));
                i_initial[2] = invCHI[6]*(HrxS_prev-Hex_x)+invCHI[7]*(HryS_prev-Hex_y)+invCHI[8]*(HrzS_prev-Hex_z) / std::sqrt(std::pow(invCHI[0]*(HrxS_prev-Hex_x)+invCHI[1]*(HryS_prev-Hex_y)+invCHI[2]*(HrzS_prev-Hex_z),2)+std::pow(invCHI[3]*(HrxS_prev-Hex_x)+invCHI[4]*(HryS_prev-Hex_y)+invCHI[5]*(HrzS_prev-Hex_z),2)+std::pow(invCHI[6]*(HrxS_prev-Hex_x)+invCHI[7]*(HryS_prev-Hex_y)+invCHI[8]*(HrzS_prev-Hex_z),2)); */
                i_initial[0] = invCHI[0]*(HrxS_prev-Hex_x)+invCHI[1]*(HryS_prev-Hex_y)+invCHI[2]*(HrzS_prev-Hex_z);
                i_initial[1] = invCHI[3]*(HrxS_prev-Hex_x)+invCHI[4]*(HryS_prev-Hex_y)+invCHI[5]*(HrzS_prev-Hex_z);
                i_initial[2] = invCHI[6]*(HrxS_prev-Hex_x)+invCHI[7]*(HryS_prev-Hex_y)+invCHI[8]*(HrzS_prev-Hex_z);
                norm_i_initial = std::sqrt(std::pow(invCHI[0]*(HrxS_prev-Hex_x)+invCHI[1]*(HryS_prev-Hex_y)+invCHI[2]*(HrzS_prev-Hex_z),2)+std::pow(invCHI[3]*(HrxS_prev-Hex_x)+invCHI[4]*(HryS_prev-Hex_y)+invCHI[5]*(HrzS_prev-Hex_z),2)+std::pow(invCHI[6]*(HrxS_prev-Hex_x)+invCHI[7]*(HryS_prev-Hex_y)+invCHI[8]*(HrzS_prev-Hex_z),2));
                i_initial[0] = i_initial[0]/norm_i_initial;
                i_initial[1] = i_initial[1]/norm_i_initial;
                i_initial[2] = i_initial[2]/norm_i_initial;
            }
            else{
              i_initial[0] = 0.0;
              i_initial[1] = 0.0;
              i_initial[2] = 0.0;
            }
            if (k > 0){
                phi_initial = std::atan2(i_initial[1],i_initial[0]);
                theta_initial = std::acos(i_initial[2]);
            }else{
                phi_initial = 0.0;
                theta_initial = 0.0;
            }
            
            // NEWTON SCHEME
            err = 1.0;
            phi = phi_initial;
            theta = theta_initial;
            iter = 0;
            while(err > 1e-3 && iter<10){
                // calc all needed stuff for the Newton step: x_k+1 = x_k - inv(H)*f'(x_k)
                // H ... Hessian, f' ... first derivative
                ux = Hex_x + chi[k]*std::cos(phi)*std::sin(theta)+offset*std::sin(phi)*std::sin(theta)+offset*std::cos(theta);
                uy = Hex_y + offset*std::cos(phi)*std::sin(theta)+chi[k]*std::sin(phi)*std::sin(theta)+offset*std::cos(theta);
                uz = Hex_z + offset*std::cos(phi)*std::sin(theta)+offset*std::sin(phi)*std::sin(theta)+chi[k]*std::cos(theta);
                uabs = std::sqrt( std::pow(ux, 2) + std::pow(uy, 2) + std::pow(uz, 2));

                du_dphi_x = chi[k]*(-1)*std::sin(phi)*std::sin(theta)+offset*std::cos(phi)*std::sin(theta)+offset*0;
                du_dphi_y = offset*(-1)*std::sin(phi)*std::sin(theta)+chi[k]*std::cos(phi)*std::sin(theta)+offset*0;
                du_dphi_z = offset*(-1)*std::sin(phi)*std::sin(theta)+offset*std::cos(phi)*std::sin(theta)+chi[k]*0;
                du_dtheta_x = chi[k]*std::cos(phi)*std::cos(theta)+offset*std::sin(phi)*std::cos(theta)+offset*(-1)*std::sin(theta);
                du_dtheta_y = offset*std::cos(phi)*std::cos(theta)+chi[k]*std::sin(phi)*std::cos(theta)+offset*(-1)*std::sin(theta);
                du_dtheta_z = offset*std::cos(phi)*std::cos(theta)+offset*std::sin(phi)*std::cos(theta)+chi[k]*(-1)*std::sin(theta);

                d2u_dphi2_x = chi[k]*(-1)*cos(phi)*sin(theta)+offset*(-1)*sin(phi)*sin(theta)+ offset*0;
                d2u_dphi2_y = offset*(-1)*cos(phi)*sin(theta)+chi[k]*(-1)*sin(phi)*sin(theta)+ offset*0;
                d2u_dphi2_z = offset*(-1)*cos(phi)*sin(theta)+offset*(-1)*sin(phi)*sin(theta)+ chi[k]*0;
                d2u_dtheta2_x = chi[k]*(-1)*cos(phi)*sin(theta)+offset*(-1)*sin(phi)*sin(theta)+offset*(-1)*cos(theta);
                d2u_dtheta2_y = offset*(-1)*cos(phi)*sin(theta)+chi[k]*(-1)*sin(phi)*sin(theta)+offset*(-1)*cos(theta);
                d2u_dtheta2_z = offset*(-1)*cos(phi)*sin(theta)+offset*(-1)*sin(phi)*sin(theta)+chi[k]*(-1)*cos(theta);

                d2u_dthetadphi_x = chi[k]*(-1)*std::sin(phi)*std::cos(theta)+offset*std::cos(phi)*std::cos(theta)+offset*0;
                d2u_dthetadphi_y = offset*(-1)*std::sin(phi)*std::cos(theta)+chi[k]*std::cos(phi)*std::cos(theta)+offset*0;
                d2u_dthetadphi_z = offset*(-1)*std::sin(phi)*std::cos(theta)+offset*std::cos(phi)*std::cos(theta)+chi[k]*0;

                //anhysteretic function + first derivative
                Man = (2*Ps_/M_PI) * std::atan(uabs/A_);
                Man_strich = (2*Ps_/(A_*M_PI)) * (1.0/(1.0+std::pow(uabs/A_, 2)));
                // first derivative of objective function after phi and theta (now vector-valued)
                c1 = Man/uabs;
                c2 = (uabs*Man_strich-Man)/(std::pow(uabs,2));
                g_strich_x = (c1*ux-MxSprev)*du_dphi_x + (c1*uy-MySprev)*du_dphi_y + (c1*uz-MzSprev)*du_dphi_z;
                g_strich_y =  (c1*ux-MxSprev)*du_dtheta_x + (c1*uy-MySprev)*du_dtheta_y + (c1*uz-MzSprev)*du_dtheta_z; 
                
                
                // second derivative (jacobi/hessian-matrix)
                g_strich_strich_11 = c1*(ux*d2u_dphi2_x+uy*d2u_dphi2_y+uz*d2u_dphi2_z) + 
                                            c1*(du_dphi_x*du_dphi_x+du_dphi_y*du_dphi_y+du_dphi_z*du_dphi_z) + 
                                            c2*((ux/uabs)*du_dphi_x+(uy/uabs)*du_dphi_y+(uz/uabs)*du_dphi_z)*(ux*du_dphi_x+uy*du_dphi_y+uz*du_dphi_z) - 
                                            (MxSprev*d2u_dphi2_x+MySprev*d2u_dphi2_y+MzSprev*d2u_dphi2_z); 
                
                g_strich_strich_12 = c1*(ux*d2u_dthetadphi_x+uy*d2u_dthetadphi_y+uz*d2u_dthetadphi_z) + 
                                            c1*(du_dphi_x*du_dtheta_x+du_dphi_y*du_dtheta_y+du_dphi_z*du_dtheta_z) + 
                                            c2*((ux/uabs)*du_dtheta_x+(uy/uabs)*du_dtheta_y+(uz/uabs)*du_dtheta_z)*(ux*du_dphi_x+uy*du_dphi_y+uz*du_dphi_z) - 
                                            (MxSprev*d2u_dthetadphi_x+MySprev*d2u_dthetadphi_y+MzSprev*d2u_dthetadphi_z);
                
                g_strich_strich_21 = c1*(ux*d2u_dthetadphi_x+uy*d2u_dthetadphi_y+uz*d2u_dthetadphi_z) + 
                                            c1*(du_dphi_x*du_dtheta_x+du_dphi_y*du_dtheta_y+du_dphi_z*du_dtheta_z) + 
                                            c2*((ux/uabs)*du_dphi_x+(uy/uabs)*du_dphi_y+(uz/uabs)*du_dphi_z)*(ux*du_dtheta_x+uy*du_dtheta_y+uz*du_dtheta_z) - 
                                            (MxSprev*d2u_dthetadphi_x+MySprev*d2u_dthetadphi_y+MzSprev*d2u_dthetadphi_z);
                
                g_strich_strich_22 = c1*(ux*d2u_dtheta2_x+uy*d2u_dtheta2_y+uz*d2u_dtheta2_z) + 
                                            c1*(du_dtheta_x*du_dtheta_x+du_dtheta_y*du_dtheta_y+du_dtheta_z*du_dtheta_z) + 
                                            c2*((ux/uabs)*du_dtheta_x+(uy/uabs)*du_dtheta_y+(uz/uabs)*du_dtheta_z)*(ux*du_dtheta_x+uy*du_dtheta_y+uz*du_dtheta_z) - 
                                            (MxSprev*d2u_dtheta2_x+MySprev*d2u_dtheta2_y+MzSprev*d2u_dtheta2_z);
                if( std::sqrt(std::pow(g_strich_x,2)+std::pow(g_strich_y,2)) < 1e-8){
                    phiNew = phi;
                    thetaNew = theta;
                }else{
                    // full Newton step
                    det_g_strich_strich = (g_strich_strich_11 * g_strich_strich_22) - (g_strich_strich_12*g_strich_strich_21);
                    inv_det_g_strich_strich = 1.0/det_g_strich_strich; // maybe check if det(J) is zero? --> error

                    phiNew = phi - inv_det_g_strich_strich*(g_strich_strich_22*g_strich_x - g_strich_strich_12*g_strich_y);
                    thetaNew = theta - inv_det_g_strich_strich*((-1)*g_strich_strich_21*g_strich_x + g_strich_strich_11*g_strich_y);
                }
                err = std::sqrt(std::pow(phiNew-phi,2)+std::pow(thetaNew-theta,2));
                phi = phiNew;
                theta = thetaNew;
                iter++;
            }
            i_correct_x = std::cos(phi)*std::sin(theta);
            i_correct_y = std::sin(phi)*std::sin(theta);
            i_correct_z = std::cos(theta);
            Hirr_x = chi[k]*i_correct_x+offset*i_correct_y+offset*i_correct_z;
            Hirr_y = offset*i_correct_x+chi[k]*i_correct_y+offset*i_correct_z;
            Hirr_z = offset*i_correct_x+offset*i_correct_y+chi[k]*i_correct_z;
            HrxS_sol[k] = Hex_x + Hirr_x;
            HryS_sol[k] = Hex_y + Hirr_y;
            HrzS_sol[k] = Hex_z + Hirr_z;
            }else{
            HrxS_sol[k] = Hex_x;
            HryS_sol[k] = Hex_y;
            HrzS_sol[k] = Hex_z;
            }
        normHrev = std::sqrt(std::pow(HrxS_sol[k],2)+ std::pow(HryS_sol[k],2)+std::pow(HrzS_sol[k],2));
        if( normHrev > 1e-20){
            MxS_sol[k] = (2.0 * Ps_/M_PI) * std::atan(normHrev/A_) * HrxS_sol[k]/normHrev;
            MyS_sol[k] = (2.0 * Ps_/M_PI) * std::atan(normHrev/A_) * HryS_sol[k]/normHrev;
            MzS_sol[k] = (2.0 * Ps_/M_PI) * std::atan(normHrev/A_) * HrzS_sol[k]/normHrev;
        }else{
            MxS_sol[k] = 0.0;
            MyS_sol[k] = 0.0;
            MzS_sol[k] = 0.0;
        }
        Px += weight[k] * MxS_sol[k];
        Py += weight[k] * MyS_sol[k];
        Pz += weight[k] * MzS_sol[k];
        }
        LOG_DBG3(eb)<< "HrxS_sol" << HrxS_sol.ToString();
        LOG_DBG3(eb)<< "HryS_sol" << HryS_sol.ToString();
        LOG_DBG3(eb)<< "HrzS_sol" << HrzS_sol.ToString();
        LOG_DBG3(eb)<< "MxS_sol" << MxS_sol.ToString();
        LOG_DBG3(eb)<< "MyS_sol" << MyS_sol.ToString();
        LOG_DBG3(eb)<< "MzS_sol" << MzS_sol.ToString();
        LOG_DBG3(eb)<< "weight" << weight.ToString();
        if(saveTmpStageVecs){
          HxS_n_tmp_[idx] = HrxS_sol;
          HyS_n_tmp_[idx] = HryS_sol;
          HzS_n_tmp_[idx] = HrzS_sol;
          MxS_n_tmp_[idx] = MxS_sol;
          MyS_n_tmp_[idx] = MyS_sol;
          MzS_n_tmp_[idx] = MzS_sol;
        }
        // output polarization
        ret.Push_back(Px);
        ret.Push_back(Py);
        ret.Push_back(Pz);
      }
    return ret; 
    }
  }


