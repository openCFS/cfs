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
      Matrix<Double> mu_FD(dim_,dim_);
      mu_FD.InitValue(0.0);

      LOG_DBG3(eb) << "\n\t HVec = " << HVec.ToString();
      if(timeStep_== 1 && globalIter_ == 1){
/*         mu[0][0] = mu_0*1e6;
        mu[0][1] = 0.0;
        mu[1][0] = 0.0;
        mu[1][1] = mu_0*1e6;
        if(dim_ == 3){
          mu[2][0] = 0.0;
          mu[2][1] = 0.0;
          mu[0][2] = 0.0;
          mu[1][2] = 0.0;
          mu[2][2] = mu_0*1e6;
        } */
        Vector<Double> M;
        M = Evaluate(HVec, false, idx);
        StdVector<Double> B_k(dim_);
        for(UInt i = 0; i < dim_; i++){
            B_k[i] = mu_0 * (HVec[i] + M[i]);
        }
        mu = EvaluateLocalMuFiniteDifferences(HVec, B_k, idx);
        mu_[idx] = mu;
        return mu;
      }

      Vector<Double> M, H, M_dummy;
      M = Evaluate(HVec, false, idx);

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


      //mu = EvaluateLocalMu(delta_H, delta_B, idx);
      //mu = EvaluateLocalMuDFP(delta_H, delta_B, idx);
      mu = EvaluateLocalMuGBM(delta_H, delta_B, idx);
      //mu = EvaluateLocalMuFiniteDifferences(HVec, B_k, idx);
      //mu = EvaluateLocalMuAnhystersisOnly(HVec, B_k, idx);


      M_dummy = Evaluate(HVec, true, idx);
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

      LOG_DBG3(eb) << "\n\t mu = " << mu.ToString();

      mu_[idx] = mu;
      M0_[idx] = M1_[idx];
      H0_[idx] = H1_[idx];

      // mark this element as computed
      hasElemSolution_[idx] = true;
      return mu;
    }

    Matrix<Double> EBHysteresis::EvaluateLocalMuAnhystersisOnly(Vector<Double> HVec, StdVector<Double> B_k, UInt idx){
      Matrix<Double> mu;
      mu.Resize(dim_,dim_);
      Matrix<Double> chi;
      chi.Resize(dim_,dim_);

      Double term1_nominator, term1_denominator, 
             term2_nominator, term2_denominator, 
             term3_nominator, term3_denominator;
      Double mu0;
      mu0 = 4*M_PI*std::pow(10,-7);

      // 2D case (x,y)
      if (dim_ == 2){

        term1_nominator =  2*Ps_*std::atan(std::sqrt((std::pow(HVec[0],2)+std::pow(HVec[1],2)))/A_);
        term1_denominator = M_PI*std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2));
        term2_nominator = 2*M_PI*std::pow(HVec[0],2)*std::atan(std::sqrt((std::pow(HVec[0],2)+std::pow(HVec[1],2)))/A_);
        term2_denominator = M_PI*std::pow(std::pow(HVec[0],2)+std::pow(HVec[1],2),3/2);
        term3_nominator = 2*Ps_*std::pow(HVec[0],2);
        term3_denominator = M_PI*A_*(std::pow(HVec[0],2)+std::pow(HVec[1],2))*(((std::pow(HVec[0],2)+std::pow(HVec[1],2))/(std::pow(A_,2)))+1);
        chi[0][0] = (term1_nominator / term1_denominator) - (term2_nominator / term2_denominator) + (term3_nominator / term3_denominator);

      if (idx == 32){
        printf("std::atan = %f\n",std::atan(std::sqrt((std::pow(HVec[0],2)+std::pow(HVec[1],2)))/A_));
      }  

        term1_nominator =  2*HVec[0]*Ps_*HVec[1];
        term1_denominator = M_PI*A_*(std::pow(HVec[0],2)+std::pow(HVec[1],2))*(((std::pow(HVec[0],2)+std::pow(HVec[1],2))/(std::pow(A_,2)))+1);
        term2_nominator = 2*HVec[0]*M_PI*HVec[1]*std::atan((std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2)))/A_);
        term2_denominator = M_PI*std::pow(std::pow(HVec[0],2)+std::pow(HVec[1],2),3/2);
        chi[1][0] = (term1_nominator / term1_denominator) - (term2_nominator / term2_denominator);
        //chi[1][0] = 0; // for all other methods we set the off diagonals also to zero, so why not here as well?
        chi[0][1] = chi[1][0];


        term1_nominator =  2*Ps_*std::atan((std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2)))/A_);
        term1_denominator = M_PI*std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2));
        term2_nominator = 2*M_PI*std::pow(HVec[1],2)*std::atan(std::sqrt((std::pow(HVec[0],2)+std::pow(HVec[1],2)))/A_);
        term2_denominator = M_PI*std::pow(std::pow(HVec[0],2)+std::pow(HVec[1],2),3/2);
        term3_nominator = 2*Ps_*std::pow(HVec[1],2);
        term3_denominator = M_PI*A_*(std::pow(HVec[0],2)+std::pow(HVec[1],2))*(((std::pow(HVec[0],2)+std::pow(HVec[1],2))/(std::pow(A_,2)))+1);
        chi[1][1] = (term1_nominator / term1_denominator) - (term2_nominator / term2_denominator) + (term3_nominator / term3_denominator);

        mu[0][0] = mu0*(chi[0][0] + 1);
        mu[1][0] = mu0*(chi[1][0] + 0);
        mu[0][1] = mu0*(chi[0][1] + 0);
        mu[1][1] = mu0*(chi[1][1] + 1);
      }else{ // 3D case (x,y,z)
        //EXCEPTION("NO 3D!!!")
        // dMx/dHx
        term1_nominator =  2*Ps_*std::atan2(std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2)),A_);
        term1_denominator = M_PI*std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2));
        term2_nominator = 2*M_PI*std::pow(HVec[0],2)*std::atan2(std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2)),A_);
        term2_denominator = M_PI*std::pow(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2),3/2);
        term3_nominator = 2*Ps_*std::pow(HVec[0],2);
        term3_denominator = M_PI*A_*(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2))*(((std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2))/(std::pow(A_,2)))+1);
        chi[0][0] = (term1_nominator / term1_denominator) - (term2_nominator / term2_denominator) + (term3_nominator / term3_denominator);
        // dMx/dHy
        term1_nominator =  2*HVec[0]*Ps_*HVec[1];
        term1_denominator = M_PI*A_*(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2))*(((std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2))/(std::pow(A_,2)))+1);
        term2_nominator = 2*HVec[0]*M_PI*HVec[1]*std::atan2(std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2))+std::pow(HVec[2],2),A_);
        term2_denominator = M_PI*std::pow(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2),3/2);
        chi[0][1] = (term1_nominator / term1_denominator) - (term2_nominator / term2_denominator);
        chi[0][1] = 0; // for all other methods we set the off diagonals also to zero, so why not here as well?
        // dMx/dHz
        term1_nominator =  2*HVec[0]*Ps_*HVec[2];
        term1_denominator = M_PI*A_*(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2))*(((std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2))/(std::pow(A_,2)))+1);
        term2_nominator = 2*HVec[0]*M_PI*HVec[2]*std::atan2(std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2))+std::pow(HVec[2],2),A_);
        term2_denominator = M_PI*std::pow(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2),3/2);
        chi[0][2] = (term1_nominator / term1_denominator) - (term2_nominator / term2_denominator);
        chi[0][2] = 0; // for all other methods we set the off diagonals also to zero, so why not here as well?

        // dMy/dHy
        term1_nominator =  2*Ps_*std::atan2(std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2)),A_);
        term1_denominator = M_PI*std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2));
        term2_nominator = 2*M_PI*std::pow(HVec[1],2)*std::atan2(std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2)),A_);
        term2_denominator = M_PI*std::pow(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2),3/2);
        term3_nominator = 2*Ps_*std::pow(HVec[1],2);
        term3_denominator = M_PI*A_*(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2))*(((std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2))/(std::pow(A_,2)))+1);
        chi[1][1] = (term1_nominator / term1_denominator) - (term2_nominator / term2_denominator) + (term3_nominator / term3_denominator);

        // dMy/dHz
        term1_nominator =  2*HVec[1]*Ps_*HVec[2];
        term1_denominator = M_PI*A_*(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2))*(((std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2))/(std::pow(A_,2)))+1);
        term2_nominator = 2*HVec[1]*M_PI*HVec[2]*std::atan2(std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2))+std::pow(HVec[2],2),A_);
        term2_denominator = M_PI*std::pow(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2),3/2);
        chi[2][1] = (term1_nominator / term1_denominator) - (term2_nominator / term2_denominator);
        chi[2][1] = 0; // for all other methods we set the off diagonals also to zero, so why not here as well?

        // dMz/dHz
        term1_nominator =  2*Ps_*std::atan2(std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2)),A_);
        term1_denominator = M_PI*std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2));
        term2_nominator = 2*M_PI*std::pow(HVec[2],2)*std::atan2(std::sqrt(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2)),A_);
        term2_denominator = M_PI*std::pow(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2),3/2);
        term3_nominator = 2*Ps_*std::pow(HVec[2],2);
        term3_denominator = M_PI*A_*(std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2))*(((std::pow(HVec[0],2)+std::pow(HVec[1],2)+std::pow(HVec[2],2))/(std::pow(A_,2)))+1);
        chi[2][2] = (term1_nominator / term1_denominator) - (term2_nominator / term2_denominator) + (term3_nominator / term3_denominator);

        // due to symmetric dM/dH tensor:
        chi[1][0] = chi[0][1];
        chi[2][0] = chi[0][2];
        chi[2][1] = chi[1][2];


        mu[0][0] = mu0*(chi[0][0] + 1);
        mu[0][1] = mu0*(chi[0][1] + 0);
        mu[0][2] = mu0*(chi[0][2] + 0);
        mu[1][0] = mu0*(chi[1][0] + 0);
        mu[1][1] = mu0*(chi[1][1] + 1);
        mu[1][2] = mu0*(chi[1][2] + 0);
        mu[2][0] = mu0*(chi[2][0] + 0);
        mu[2][1] = mu0*(chi[2][1] + 0);
        mu[2][2] = mu0*(chi[2][2] + 1);


      }
      if (idx == 32){
        printf("mu_anhyst = %f\n",mu[0][0]);
        printf("Hx = %f\n",HVec[0]);
      }  
      return mu;
    }

    Matrix<Double> EBHysteresis::EvaluateLocalMuFiniteDifferences(Vector<Double> HVec, StdVector<Double> B_k, UInt idx){
      Matrix<Double> mu;
      mu.Resize(dim_,dim_);

      Double h = 1*10^(-13);
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
          // calculate the mu tesnor (dB/dH) for every entry
          mu[0][0] = (B_incx[0] - B_k[0]) / h; mu[0][1] = 0;                        mu[0][2] = 0;
          mu[1][0] = 0;                        mu[1][1] = (B_incy[1] - B_k[1]) / h; mu[1][2] = 0;
          mu[2][0] = 0;                        mu[2][1] = 0;                        mu[2][2] = (B_incz[2] - B_k[2]) / h;

      } 
      if (idx == 32){
        printf("mu_FD = %f\n",mu[0][0]);
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
      return mu;
      //############### delete as soon as it works #######################
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
      

      // #######################################################################################
      // #######################################################################################

      /* just for TEAMProblem32 parameter, since using the mat.xml file only evenly 
      distributed w(k) can be generated and no double Langevin function as anhystersis curve
      can be defined (for normal simulations just comment it out) !! */
      Double M_an,aa, MSa, ab, MSb, J_an,coth_La,coth_Lb;

      // weights from Jacques
      /* chi[0] = 0; chi[1] = 55; chi[2] = 199;
      weight[0] = 0.03; weight[1] = 0.92; weight[2] = 0.05; */ 

      // my weights k = 5
      /* chi[0] = 0; chi[1] = 50.97; chi[2] = 118.73; chi[3] = 222.41; chi[4] = 709.7;
      weight[0] = 0.028703738149827; weight[1] = 0.847068065612296; weight[2] = 0.092764375057132; weight[3] = 0.029978499149262; weight[4] = 0.001485322031482; */ 

      // my weights k = 10
      /* chi[0]=0; chi[1]=35.9778; chi[2]=64.1895; chi[3]=107.8238; 
      chi[4]=151.2311; chi[5]=188.1423; chi[6]=240.4219; chi[7]=160.602; 
      chi[8]=460.8971; chi[9]=709.6528;  
      weight[0]=0.028704; weight[1]=0.39699; weight[2]=0.45008; weight[3]=0.069457; 
      weight[4]=0.023307; weight[5]=0.012757; weight[6]=0.015412; weight[7]=0.00090465; 
      weight[8]=0.00090465; weight[9]=0.0014853; */

      // my weights k = 50
      /* chi[0]=0; chi[1]=7.7505; chi[2]=21.3061; chi[3]=34.1291; 
      chi[4]=44.9371; chi[5]=54.4629; chi[6]=64.2734; chi[7]=73.2565; 
      chi[8]=84.8699; chi[9]=90.5018; chi[10]=103.2363; chi[11]=119.2992; 
      chi[12]=125.3339; chi[13]=129.953; chi[14]=141.0746; chi[15]=153.6859; 
      chi[16]=167.7871; chi[17]=183.3781; chi[18]=198.7852; chi[19]=172.6219; 
      chi[20]=175.3508; chi[21]=180.0113; chi[22]=186.6034; chi[23]=195.1272; 
      chi[24]=205.5826; chi[25]=217.9696; chi[26]=232.2882; chi[27]=248.5385; 
      chi[28]=266.7204; chi[29]=286.8339; chi[30]=308.8791; chi[31]=331.1456; 
      chi[32]=-51.2161; chi[33]=-9.241; chi[34]=32.1366; chi[35]=72.9166; 
      chi[36]=113.0991; chi[37]=152.684; chi[38]=191.6713; chi[39]=230.061; 
      chi[40]=267.8532; chi[41]=305.0478; chi[42]=341.6448; chi[43]=377.6443; 
      chi[44]=413.0462; chi[45]=447.8505; chi[46]=482.0572; chi[47]=515.6664; 
      chi[48]=548.678; chi[49]=581.0921; chi[50]=709.4941; 
      weight[0]=0.028704; weight[1]=0.026403; weight[2]=0.054303; weight[3]=0.12049; 
      weight[4]=0.19363; weight[5]=0.18951; weight[6]=0.13623; weight[7]=0.070348; 
      weight[8]=0.054056; weight[9]=0.020644; weight[10]=0.020644; weight[11]=0.020644; 
      weight[12]=0.0085211; weight[13]=0.0056961; weight[14]=0.0056961; weight[15]=0.0056961; 
      weight[16]=0.0056961; weight[17]=0.0056961; weight[18]=0.0052661; weight[19]=0.0014384; 
      weight[20]=0.0014384; weight[21]=0.0014384; weight[22]=0.0014384; weight[23]=0.0014384; 
      weight[24]=0.0014384; weight[25]=0.0014384; weight[26]=0.0014384; weight[27]=0.0014384; 
      weight[28]=0.0014384; weight[29]=0.0014384; weight[30]=0.0014384; weight[31]=0.0013878; 
      weight[32]=0.00011055; weight[33]=0.00011055; weight[34]=0.00011055; weight[35]=0.00011055; 
      weight[36]=0.00011055; weight[37]=0.00011055; weight[38]=0.00011055; weight[39]=0.00011055; 
      weight[40]=0.00011055; weight[41]=0.00011055; weight[42]=0.00011055; weight[43]=0.00011055; 
      weight[44]=0.00011055; weight[45]=0.00011055; weight[46]=0.00011055; weight[47]=0.00011055; 
      weight[48]=0.00011055; weight[49]=0.00011055; weight[50]=0.0014874; */

      // anhysteresis 
      aa = 9.082; MSa = 0.792; ab = 137.121; MSb = 0.791;

      // #######################################################################################
      // #######################################################################################

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
            /* coth_La = std::cosh(HrS/aa)/std::sinh(HrS/aa);
            coth_Lb = std::cosh(HrS/ab)/std::sinh(HrS/ab);
            J_an = 1*(coth_La - aa/HrS);
            M_an = J_an/(4*M_PI*1e-7); */
            MxS_sol[k] = (2.0 * Ps_/M_PI) * std::atan(HrS/A_) * HrxS_sol[k]/HrS;
            MyS_sol[k] = (2.0 * Ps_/M_PI) * std::atan(HrS/A_) * HryS_sol[k]/HrS;
            MzS_sol[k] = (2.0 * Ps_/M_PI) * std::atan(HrS/A_) * HrzS_sol[k]/HrS;
            /* MxS_sol[k] = M_an * HrxS_sol[k]/HrS;
            MyS_sol[k] = M_an * HryS_sol[k]/HrS;
            MzS_sol[k] = M_an * HrzS_sol[k]/HrS; */
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
        ret.Push_back(Px);
        ret.Push_back(Py);
        ret.Push_back(Pz);
      }
    return ret; 
    }
  }


