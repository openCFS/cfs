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
      mu_0 = 1.256637061e-06;

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
      // This method gives the rho_art for the h-form.

      Double penaltyParameter = 1e-6; // default value for n

      /* Matrix<Double> mu = ComputeTensorialMaterialParameter(HVec, ElemNum);
      Double num1 = mu[0][0]; Double num2 = mu[1][1]; Double num3 = mu[2][2];
      if (num1 >= num2 && num1 >= num3) {
        return num1*std::pow(10,penaltyParameter/2);
      } else if (num2 >= num1 && num2 >= num3) {
        return num2*std::pow(10,penaltyParameter/2);
      } else {
        return num3*std::pow(10,penaltyParameter/2);
      } */

      /* Double normH; normH = std::sqrt(HVec[0]*HVec[0] + HVec[1]*HVec[1] + HVec[2]*HVec[2]);
      Double mu0;Double material_scalar; Double chi;
      mu0 = 1.256637061e-06;;
      chi = (2*A_*Ps_)/(M_PI*(normH*normH+A_*A_));
      material_scalar = mu0*(1+chi); 
      return material_scalar*std::pow(10,penaltyParameter/2); */

      // SVD approach
      Matrix<Double> A = ComputeTensorialMaterialParameter(HVec, ElemNum);
      Matrix<Double> AtA(dim_,dim_);
      for(UInt idx=0;idx<3;idx++){
        for(UInt jdx=0;jdx<3;jdx++){
            // Compute the (i, j)-th entry of AtA
            for (int kdx = 0; kdx < 3; ++kdx) {
              AtA[idx][jdx] += A[kdx][idx]*A[kdx][jdx];
            }
        }
      }
      Double a,b,c;
      // Compute the coefficients of the characteristic polynomial directly
      a = 1.0;
      b = -1.0 * (A[0][0] + A[1][1] + A[2][2]);
      c = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]) +
          A[1][1] * A[2][2] * A[0][0] +
          A[2][2] * A[0][0] * A[1][1] -
          A[0][1] * A[1][0] * A[2][2] -
          A[1][2] * A[2][1] * A[0][0] -
          A[2][0] * A[0][2] * A[1][1];

      Double discriminant = b * b - 4 * a * c;
      Double lambda1 = (-b + std::sqrt(discriminant)) / (2 * a);
      Double lambda2 = (-b - std::sqrt(discriminant)) / (2 * a);
      Double lambda3 = (A[0][0] + A[1][1] + A[2][2]) - lambda1 - lambda2;

      Double svd1 = std::sqrt(std::abs(lambda1)); Double svd2 = std::sqrt(std::abs(lambda2)); Double svd3 = std::sqrt(std::abs(lambda3));

      if (svd1 >= svd2 && svd1 >= svd3) {
          return svd1*std::pow(10,penaltyParameter/2);
      } else if (svd2 >= svd1 && svd2 >= svd3) {
          return svd2*std::pow(10,penaltyParameter/2);
      } else {
          return svd3*std::pow(10,penaltyParameter/2);
      }
    }

    Vector<Double> EBHysteresis::GetFluxDensity(Vector<Double> HVec, const Integer ElemNum) {
      Vector<Double> B(dim_);
      UInt idx = ElemNum2Idx_[ElemNum];

      LOG_DBG3(eb) << "\n\t HVec = " << HVec.ToString();
      Vector<Double> M;
      M = Evaluate(HVec, false, idx);
      LOG_DBG3(eb) << "\n\t M = " << M.ToString();
      
      for(UInt i = 0; i < dim_; ++i){
        B[i] = mu_0*(HVec[i] + M[i]);
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
      if(hasElemSolution_[idx] == true){
        return mu_[idx];
      }

      Vector<Double> M;
      Vector<Double> M_dummy;
      Matrix<Double> mu(dim_,dim_); mu.InitValue(0.0);

      M = Evaluate(HVec, false, idx);

      // get the values for H and M from the last timestep and get deltaH and deltaB
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

   
      if (numS_ > 1 ){ // hysteretic case
        mu = EvaluateLocalMuFiniteDifferences(HVec, B_k, idx);
      } else { // nonlinear case (only anhysteresis)
        mu = EvaluateLocalMuAnhystersisOnly(HVec, idx);
      }
      mu_[idx] = mu;

      M_dummy = Evaluate(HVec, true, idx);

      // mark this element as computed
      hasElemSolution_[idx] = true;
      return mu;
    }

    Matrix<Double> EBHysteresis::EvaluateLocalMuAnhystersisOnly(Vector<Double> HVec, UInt idx){
      Matrix<Double> mu; mu.Resize(dim_,dim_);
      Matrix<Double> chi; chi.Resize(dim_,dim_);
      Matrix<Double> T2; T2.Resize(dim_,dim_);
      Matrix<Double> identity; identity.Resize(dim_,dim_); 

      Double t0, t1, t3, factor1,factor2,factor3;
      Double mu0 = 1.256637061e-06;;

      // 2D case (x,y)
      if (dim_ == 2){
        t0 = std::sqrt(HVec[0]*HVec[0] + HVec[1]*HVec[1]);
        t1 = t0/A_;
        t3 = std::atan(t1);
        factor1 = (2*Ps_/(t1*t1+1))/(M_PI*t0*t0*A_);
        factor2 = (2*Ps_*t3)/(M_PI*t0*t0*t0);
        factor3 = (2*Ps_*t3)/(M_PI*t0);

        for (int i = 0; i < dim_; ++i){
          for (int j = 0; j < dim_; ++j){
              identity[i][j] = (i == j) ? 1 : 0;
              T2[i][j] = HVec[i]*HVec[j];
              chi[i][j] = factor1*T2[i][j] - factor2*T2[i][j] + factor3*identity[i][j];
              if (std::isnan(chi[i][j])){
                chi[i][j] = 0;
              }
              mu[i][j] = mu0*(identity[i][j] + chi[i][j]);
          }
        }
      }else{ // 3D case (x,y,z)
        HVec[0] = (HVec[0] == 0) ? 1e-12 : HVec[0];
        HVec[1] = (HVec[1] == 0) ? 1e-12 : HVec[1];
        HVec[2] = (HVec[2] == 0) ? 1e-12 : HVec[2];
        t0 = std::sqrt(HVec[0]*HVec[0] + HVec[1]*HVec[1] + HVec[2]*HVec[2]);
        t1 = t0/A_;
        t3 = std::atan(t1);
        factor1 = (2*Ps_/(t1*t1+1))/(M_PI*t0*t0*A_);
        factor2 = (2*Ps_*t3)/(M_PI*t0*t0*t0);
        factor3 = (2*Ps_*t3)/(M_PI*t0);

        for (int i = 0; i < dim_; ++i){
          for (int j = 0; j < dim_; ++j){
              identity[i][j] = (i == j) ? 1 : 0;
              T2[i][j] = HVec[i]*HVec[j];
              chi[i][j] = factor1*T2[i][j] - factor2*T2[i][j] + factor3*identity[i][j];
              if (std::isnan(chi[i][j])){
                chi[i][j] = 0;
              }
              mu[i][j] = mu0*(identity[i][j] + chi[i][j]);
          }
        }
      } 
        /* if (idx == 1){
          std::cout << "Hx: " << HVec[0] << ", Hy: " << HVec[1] << ", Hz: " << HVec[2] << std::endl;
          std::cout << "mu[0][0]: " << mu[0][0] << ", mu[0][1]: " << mu[0][1] << ", mu[0][2]: " << mu[0][2] << std::endl;
          std::cout << "mu[1][0]: " << mu[1][0] << ", mu[1][1]: " << mu[1][1] << ", mu[1][2]: " << mu[1][2] << std::endl;
          std::cout << "mu[2][0]: " << mu[2][0] << ", mu[2][1]: " << mu[2][1] << ", mu[2][2]: " << mu[2][2] << std::endl;
        } */
      return mu;
    }

    Matrix<Double> EBHysteresis::EvaluateLocalMuFiniteDifferences(Vector<Double> HVec, StdVector<Double> B_k, UInt idx){
      Matrix<Double> mu;
      mu.Resize(dim_,dim_);

      Double h = 1*10^(-7);
      Vector<Double> eh_x(dim_);
      Vector<Double> eh_y(dim_);
      Vector<Double> eh_z(dim_);
      Vector<Double> H_incx(dim_);
      Vector<Double> H_incy(dim_);
      Vector<Double> H_incz(dim_);

      Vector<Double> M_incx(dim_);
      Vector<Double> M_incy(dim_);
      Vector<Double> M_incz(dim_);
      Vector<Double> B_incx(dim_);
      Vector<Double> B_incy(dim_);
      Vector<Double> B_incz(dim_);

      // 2D case (x,y)
      if (dim_ == 2){

        eh_x[0] = h; eh_x[1] = 0;
        eh_y[0] = 0; eh_y[1] = h;
        for(UInt i = 0; i < dim_; i++){
          H_incx[i] = HVec[i] + eh_x[i];
          H_incy[i] = HVec[i] + eh_y[i];
        }
        // Evaluate EB Hysteresis Model for the increment AND without updating the stage values!!
        M_incx = Evaluate(H_incx, false, idx);
        M_incy = Evaluate(H_incy, false, idx);
        for(UInt i = 0; i < dim_; i++){
          B_incx[i] = mu_0 * (H_incx[i] + M_incx[i]);
          B_incy[i] = mu_0 * (H_incy[i] + M_incy[i]);
        }
        // calculate the mu tesnor (dB/dH) for every entry
        mu[0][0] = (B_incx[0] - B_k[0]) / h; mu[0][1] = 0;
        mu[1][0] = 0; mu[1][1] = (B_incy[1] - B_k[1]) / h;

      }else{ // 3D case (x,y,z)

          eh_x[0] = h; eh_x[1] = 0; eh_x[2] = 0;
          eh_y[0] = 0; eh_y[1] = h; eh_y[2] = 0;
          eh_z[0] = 0; eh_z[1] = 0; eh_z[2] = h;
          for(UInt i = 0; i < dim_; i++){
            H_incx[i] = HVec[i] + eh_x[i];
            H_incy[i] = HVec[i] + eh_y[i];
            H_incz[i] = HVec[i] + eh_z[i];
          }
          // Evaluate EB Hysteresis Model for the increment AND without updating the stage values!!
          M_incx = Evaluate(H_incx, false, idx);
          M_incy = Evaluate(H_incy, false, idx);
          M_incz = Evaluate(H_incz, false, idx);
          for(UInt i = 0; i < dim_; i++){
            B_incx[i] = mu_0 * (H_incx[i] + M_incx[i]);
            B_incy[i] = mu_0 * (H_incy[i] + M_incy[i]);
            B_incz[i] = mu_0 * (H_incz[i] + M_incz[i]);
          }
          // calculate the mu tensor (dB/dH) for every entry
          mu[0][0] = (B_incx[0] - B_k[0]) / h; mu[0][1] = 0;                        mu[0][2] = 0;
          mu[1][0] = 0;                        mu[1][1] = (B_incy[1] - B_k[1]) / h; mu[1][2] = 0;
          mu[2][0] = 0;                        mu[2][1] = 0;                        mu[2][2] = (B_incz[2] - B_k[2]) / h;

      } 
        if (idx == 1){
          std::cout << "FDHx: " << HVec[0] << ", FDHy: " << HVec[1] << ", FDHz: " << HVec[2] << std::endl;
          std::cout << "mu[0][0]: " << mu[0][0] << ", mu[0][1]: " << mu[0][1] << ", mu[0][2]: " << mu[0][2] << std::endl;
          std::cout << "mu[1][0]: " << mu[1][0] << ", mu[1][1]: " << mu[1][1] << ", mu[1][2]: " << mu[1][2] << std::endl;
          std::cout << "mu[2][0]: " << mu[2][0] << ", mu[2][1]: " << mu[2][1] << ", mu[2][2]: " << mu[2][2] << std::endl;
        }
      return mu;
    }

    Matrix<Double> EBHysteresis::EvaluateLocalMu(StdVector<Double> dH, StdVector<Double> dB, UInt idx){
      Matrix<Double> mu;
      mu.Resize(dim_,dim_);
      double offset = 1e-2;

      /* if(dim_ == 2){
        dH[0] = dH[0] + offset; dH[1] = dH[1] + offset;
        dB[0] = dB[0] + offset; dB[1] = dB[1] + offset;
      } */
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
        for (UInt i = 0; i < dim_; i++){
            if (std::isinf(mu[i][i]) || std::isnan(mu[i][i])){
              Matrix<Double> e = mu_[idx];
                mu[i][i] = e[i][i]; //e[i][i]; //mu_0;
            } 
        }
      if (idx == 32){
        printf("mu_simpleFD = %f\n",mu[0][0]);
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
                mu[i][i] = e[i][i]; //e[i][i]; //mu_0;
            } 
        }
      if (idx == 32){
        printf("mu_DFPxx = %f\n",mu[0][0]);
        printf("mu_DFPxy = %f\n",mu[0][1]);
        printf("mu_DFPyx = %f\n",mu[1][0]);
        printf("mu_DFPyy = %f\n",mu[1][1]);
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

        /* mu[0][1] = 0.0;
        mu[1][0] = 0.0;
        mu[2][0] = 0.0;
        mu[0][2] = 0.0;
        mu[1][2] = 0.0;
        mu[2][1] = 0.0; */
        for (UInt i = 0; i < dim_; i++){
            if (std::isinf(mu[i][i]) || std::isnan(mu[i][i])){
              Matrix<Double> e = mu_[idx];
                mu[i][i] = e[i][i];//mu_0; //e[i][i]; //mu_0;
            } 
        }
      }
      if (idx == 32){
        printf("mu_DFP = %f\n",mu[0][0]);
      }
      return mu;
    }

    Matrix<Double> EBHysteresis::EvaluateLocalMuGBM(StdVector<Double> dH, StdVector<Double> dB, UInt idx){
      Matrix<Double> mu;
      mu.Resize(dim_,dim_);

      Vector<Double> dD_minus_epsilon_km1_times_dE(dim_);
      Vector<Double> dHvec(dim_);
      Matrix<Double> rightM(dim_,dim_);
      double offset = 1e-13;

      if(dim_ == 2){
        dH[0] = dH[0] + offset; dH[1] = dH[1] + offset;
        dB[0] = dB[0] + offset; dB[1] = dB[1] + offset;
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

        /* mu[0][1] = 0.0;
        mu[1][0] = 0.0; */
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
      //################### just for checking some things #################
      if (idx == 15){
        printf("mu_Broyden = %f\n",mu[0][0]);
      }
      //############### delete as soon as it works #######################
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
        chi[i] = chi_factor_ * i / (numS_ - 1) + 1e-22;
        chi[i] = chi_factor_ * i / (numS_ - 1) + 0.001;
        chi[i] = chi_factor_ * i / (numS_ - 1) + 1e-22;
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

      }else{   // 3D Version (VPM)
        // volumetric weight for each pseudo particle
        //EXCEPTION("NO 3D!!!")
        StdVector<Double> error, dir, HrxS_sol, HryS_sol, HrzS_sol, MxS_sol, MyS_sol, MzS_sol;
        StdVector<UInt> numIter;
        error.Resize(numS_, 0.0);
        dir.Resize(numS_, 0.0);
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

        Double HrxS_prev, HryS_prev, HrzS_prev, MxSprev, MySprev, MzSprev, phi, err, ux, uy,
              uabs, ux1, uy1, ux2, uy2, Man, Man1, g1, g2, phiNew, HrS, Px, Py, Pz,
              condition1, theta, i_correct_x,i_correct_y,i_correct_z,Hirr_x,Hirr_y,Hirr_z,coth_La,coth_Lb,J_an;
        UInt iter;
        StdVector<Double>& HxS_prev = HxS_n_[idx];
        StdVector<Double>& HyS_prev = HyS_n_[idx];
        StdVector<Double>& HzS_prev = HzS_n_[idx];
        StdVector<Double>& MxS_prev = MxS_n_[idx];
        StdVector<Double>& MyS_prev = MyS_n_[idx];
        StdVector<Double>& MzS_prev = MzS_n_[idx];


        for(int k = 0; k < numS_; k++){
          HrxS_prev = HxS_prev[k];
          HryS_prev = HyS_prev[k];
          HrzS_prev = HzS_prev[k];
          MxSprev = MxS_prev[k];
          MySprev = MyS_prev[k];
          MzSprev = MzS_prev[k];
          condition1 = ( std::pow((Hex_x - HrxS_prev)/chi[k], 2) + std::pow((Hex_y - HryS_prev)/chi[k], 2) + std::pow((Hex_z - HrzS_prev)/chi[k], 2) );
          if( condition1 <= 1.0){
            HrxS_sol[k] = HrxS_prev;
            HryS_sol[k] = HryS_prev;
            HrzS_sol[k] = HrzS_prev;
          }else{
            // dissipation now reached (just arrived at the border of the sphere)
            if(k == 0){
              HrxS_sol[k] = Hex_x;
              HryS_sol[k] = Hex_y;
              HrzS_sol[k] = Hex_z;
            }else{
              // use the direction of the vector play model as initial direction
              Hirr_x = chi[k]* (Hex_x - HrxS_prev) / std::sqrt((std::pow((Hex_x - HrxS_prev),2) + std::pow((Hex_y - HryS_prev),2)) + std::pow((Hex_z - HrzS_prev),2));
              Hirr_y = chi[k]* (Hex_y - HryS_prev) / std::sqrt((std::pow((Hex_x - HrxS_prev),2) + std::pow((Hex_y - HryS_prev),2)) + std::pow((Hex_z - HrzS_prev),2));
              Hirr_z = chi[k]* (Hex_z - HrzS_prev) / std::sqrt((std::pow((Hex_x - HrxS_prev),2) + std::pow((Hex_y - HryS_prev),2)) + std::pow((Hex_z - HrzS_prev),2));
              HrxS_sol[k] = Hex_x - Hirr_x;
              HryS_sol[k] = Hex_y - Hirr_y;
              HrzS_sol[k] = Hex_z - Hirr_z;
            }
          }  
          HrS = std::sqrt(std::pow(HrxS_sol[k], 2) + std::pow(HryS_sol[k], 2) + std::pow(HrzS_sol[k], 2));
          if( std::sqrt(std::pow(HrS,2)) > 1.0e-12){
            MxS_sol[k] = (2.0 * Ps_/M_PI) * std::atan(HrS/A_) * HrxS_sol[k]/HrS;
            MyS_sol[k] = (2.0 * Ps_/M_PI) * std::atan(HrS/A_) * HryS_sol[k]/HrS;
            MzS_sol[k] = (2.0 * Ps_/M_PI) * std::atan(HrS/A_) * HrzS_sol[k]/HrS;
          }else{
            MxS_sol[k] = 0.0;
            MyS_sol[k] = 0.0;
            MzS_sol[k] = 0.0;
          }
        }
        Px = 0.0;
        Py = 0.0;
        Pz = 0.0;
        for(int k = 0; k < numS_; k++){
          Px += weight[k] * MxS_sol[k];
          Py += weight[k] * MyS_sol[k];
          Pz += weight[k] * MzS_sol[k];
        }
        if(saveTmpStageVecs){
          HxS_n_tmp_[idx] = HrxS_sol;
          HyS_n_tmp_[idx] = HryS_sol;
          HzS_n_tmp_[idx] = HrzS_sol;
          MxS_n_tmp_[idx] = MxS_sol;
          MyS_n_tmp_[idx] = MyS_sol;
          MzS_n_tmp_[idx] = MzS_sol;
        }
        /* Px = 0.0;
        Py = 0.0;
        Pz = 0.0;
        Double normH;
        normH = std::sqrt(Hn[0]*Hn[0] + Hn[1]*Hn[1] + Hn[2]*Hn[2]);
        Px = (2.0*Ps_/M_PI)*std::atan(normH/A_) * Hn[0]/normH;
        Py = (2.0*Ps_/M_PI)*std::atan(normH/A_) * Hn[1]/normH;
        Pz = (2.0*Ps_/M_PI)*std::atan(normH/A_) * Hn[2]/normH;
        if (std::isnan(Px)){
          Px = 0;
        }
        if (std::isnan(Py)){
          Py = 0;
        }
        if (std::isnan(Pz)){
          Pz = 0;
        } */
        ret.Push_back(Px);
        ret.Push_back(Py);
        ret.Push_back(Pz);
      }
    return ret; 
    }
  }
