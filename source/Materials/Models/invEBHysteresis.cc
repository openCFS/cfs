#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "DataInOut/Logging/LogConfigurator.hh"

#include "invEBHysteresis.hh"
#include "Model.hh"

#include "MatVec/Vector.hh"


#include "Utils/mathParser/mathParser.hh"
#include "Domain/Domain.hh"
#include "Driver/SolveSteps/SolveStepEB.hh"

namespace CoupledField {

DEFINE_LOG(inveb, "invEBHysteresis")


    invEBHysteresis::invEBHysteresis() : Model(),
    numElems_{0}, 
    Js_{0}, A_{0}, numS_{0},chi_factor_{0}, jacobian_method_{0},
    mp_{nullptr}, isFirstTimeFinished_{0},
    timeStep_{0}, globalIter_{0},
    isMH_{false}
    {}

    invEBHysteresis::~invEBHysteresis() {
    }

    void invEBHysteresis::Init(std::map<std::string, double> ParameterMap, shared_ptr<ElemList> entityList, UInt dim) {

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
      Js_ = ParameterMap["Js"];
      A_ = ParameterMap["A"];
      numS_ = ParameterMap["numS"];
      chi_factor_ = ParameterMap["chi_factor"];
      jacobian_method_ = ParameterMap["jacobian_method"];

      isMH_ = ParameterMap["isMH"];
      if(isMH_ == 1.0){
        varHandle_="cacheResult";
      } else {
        varHandle_="step";
      }
      mu_0 = 1.256637061435917e-06;

      Js_ = ParameterMap["Js"];
      A_ = ParameterMap["A"];
      numS_ = ParameterMap["numS"];
      chi_factor_ = ParameterMap["chi_factor"];
      jacobian_method_ = ParameterMap["jacobian_method"];

      B0_.Resize(numElems_, StdVector<Double>(dim_));
      B1_.Resize(numElems_, StdVector<Double>(dim_));
      H0_.Resize(numElems_, StdVector<Double>(dim_));
      H1_.Resize(numElems_, StdVector<Double>(dim_));

      J_x_k_prev_.Resize(numElems_, StdVector<Double>(numS_));
      J_y_k_prev_.Resize(numElems_, StdVector<Double>(numS_));
      J_z_k_prev_.Resize(numElems_, StdVector<Double>(numS_));

      J_x_k_n_.Resize(numElems_, StdVector<Double>(numS_));
      J_y_k_n_.Resize(numElems_, StdVector<Double>(numS_));
      J_z_k_n_.Resize(numElems_, StdVector<Double>(numS_));

      J_x_k_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));
      J_y_k_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));
      J_z_k_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));

      nu_.Resize(numElems_, Matrix<Double>(dim_,dim_));

      isFirstTime_.Resize(numElems_);
      isFirstTime_.Init(1);

      isFirstTimeFinished_ = false;

      timeStep_ = 1;

      mp_ = domain->GetMathParser();
      globalIter_ = 0;
      solve_step_eb_ = dynamic_cast<SolveStepEB*>(domain->GetBasePDE()->GetSolveStep());
      if ( !solve_step_eb_ ) {
        EXCEPTION("Inverse Hysteresis needs solve_step_eb_");
      }

      hasElemSolution_.Resize(numElems_, false);
    }

    Vector<Double> invEBHysteresis::GetFieldIntensity(Vector<Double> BVec, const Integer ElemNum) {
      // define needed variables
      Vector<Double> H(dim_);
      UInt idx = ElemNum2Idx_[ElemNum];

      // compute H
      if (this->solve_step_eb_->update_polarization_eb_ == true) {
        H = Evaluate(BVec, true, idx);
      } else {
        H = Evaluate(BVec, false, idx);
      }
        
      // return
      return H;
    }

    Matrix<Double> invEBHysteresis::ComputeTensorialMaterialParameter(Vector<Double> BVec, const Integer ElemNum) {
      #pragma omp critical 
      {
        if(globalIter_ != mp_->GetExprVars(MathParser::GLOB_HANDLER, "iterationCounter")){
          globalIter_ = mp_->GetExprVars(MathParser::GLOB_HANDLER, "iterationCounter");
          hasElemSolution_.Init(false); //if there is a new iteration, save the values from the previous iteration
        }
        if( (timeStep_ != mp_->GetExprVars(MathParser::GLOB_HANDLER, varHandle_))){
          hasElemSolution_.Init(false);
          J_x_k_n_ = J_x_k_n_tmp_;
          J_y_k_n_ = J_y_k_n_tmp_;
          J_z_k_n_ = J_z_k_n_tmp_;
          timeStep_ = mp_->GetExprVars(MathParser::GLOB_HANDLER, varHandle_);
        }
      }
      // define needed variables
      UInt idx = ElemNum2Idx_[ElemNum];
      Matrix<Double>     nu(dim_,dim_); nu.InitValue(0.0);
      StdVector<Double>& B0 = B0_[idx];
      StdVector<Double>& H0 = H0_[idx];
      Vector<Double>  H(dim_);
      StdVector<Double>  H_old(dim_);
      Vector<Double>  H_dummy(dim_);
      Vector<Double>  delta_H(dim_);
      Vector<Double>  delta_B(dim_);

      // algorithm to determine nu = dH/dB;
      if(hasElemSolution_[idx] == true){ // nu of this element has already been computed
        return nu_[idx];
      }
      if( globalIter_ == 1){ // starting value at nu = diag(1/(mu0*1000));
        for(UInt i = 0; i < dim_; i++) {
          for(UInt j = 0; j < dim_; j++) {
            if (i == j) {
              nu[i][j] = 1/(mu_0*1000);
            } else {
              nu[i][j] = 0;
            }
          }
        }
        nu_[idx] = nu;
        hasElemSolution_[idx] = true;
        return nu;
      }

      H = Evaluate(BVec, false, idx); // get current H
      for(UInt i = 0; i < dim_; i++) { // get the values for H and M from the last timestep and get deltaH and deltaB
        H_old[i] = H0[i]; // get old H
        delta_B[i] = BVec[i] - B0[i];
        delta_H[i] = H[i] - H_old[i];
        B1_[idx][i] = BVec[i];
        H1_[idx][i] = H[i];
      }
      nu = EvaluateLocalNuBFGS(delta_H, delta_B, idx); // determine nu
      // nu = EvaluateLocalMuGBM(delta_B, delta_H, idx);
      
      nu_[idx] = nu; // store all quantities for next iteration step
      H0_[idx] = H1_[idx];
      B0_[idx] = B1_[idx];

      nu_[idx] = nu;
      hasElemSolution_[idx] = true; // mark this element as computed

      // return
      return nu;
    }


    Matrix<Double> invEBHysteresis::EvaluateLocalNuBFGS(Vector<Double> dH, Vector<Double> dB, UInt idx){
      // define needed variables
      Matrix<Double> yyT(dim_, dim_);
      Matrix<Double> BxBxT(dim_, dim_);
      Matrix<Double> B_k1(dim_, dim_);
      Matrix<Double> B = nu_[idx];
      Double yTx;
      Double xTBx;
      Vector<Double> y = dH;
      Vector<Double> x = dB;
      Vector<Double> Bx(dim_);

      // update nu via BFGS formula
      if(dim_ == 2){
        // yyT (outer product)
        yyT[0][0] = y[0]*y[0]; yyT[0][1] = y[0]*y[1]; 
        yyT[1][0] = y[1]*y[0]; yyT[1][1] = y[1]*y[1]; 

        // yTx (inner product)
        yTx = (y[0]*x[0]) + (y[1]*x[1]);

        // xTBx (inner product)
        Bx[0] = B[0][0]*x[0] + B[0][1]*x[1];
        Bx[1] = B[1][0]*x[0] + B[1][1]*x[1];
        xTBx = (x[0]*Bx[0]) + (x[1]*Bx[1]);

        // BxBxT (outer product)
        BxBxT[0][0] = Bx[0]*Bx[0]; BxBxT[0][1] = Bx[0]*Bx[1];
        BxBxT[1][0] = Bx[1]*Bx[0]; BxBxT[1][1] = Bx[1]*Bx[1];
        
        // construct everything
        B_k1[0][0] = B[0][0] + yyT[0][0]/yTx - BxBxT[0][0]/xTBx; 
        B_k1[0][1] = B[0][1] + yyT[0][1]/yTx - BxBxT[0][1]/xTBx;

        B_k1[1][0] = B[1][0] + yyT[1][0]/yTx - BxBxT[1][0]/xTBx; 
        B_k1[1][1] = B[1][1] + yyT[1][1]/yTx - BxBxT[1][1]/xTBx;

        if ( (std::isnan(B_k1[0][0])) || (std::isnan(B_k1[1][1])) || (std::isnan(B_k1[0][1])) || (std::isnan(B_k1[1][0])) ) {
          B_k1[0][0] = 1/mu_0;
          B_k1[1][1] = 0;
          B_k1[0][1] = 0;
          B_k1[1][0] = 1/mu_0;
        }
        if ( (std::isinf(B_k1[0][0])) || (std::isinf(B_k1[1][1])) || (std::isinf(B_k1[0][1])) || (std::isinf(B_k1[1][0])) ) {
          B_k1[0][0] = 1/mu_0;
          B_k1[1][1] = 0;
          B_k1[0][1] = 0;
          B_k1[1][0] = 1/mu_0;
        }
      } else { // 3-D version
        // yyT (outer product)
        yyT[0][0] = y[0]*y[0]; yyT[0][1] = y[0]*y[1]; yyT[0][2] = y[0]*y[2];
        yyT[1][0] = y[1]*y[0]; yyT[1][1] = y[1]*y[1]; yyT[1][2] = y[1]*y[2];
        yyT[2][0] = y[2]*y[0]; yyT[2][1] = y[2]*y[1]; yyT[2][2] = y[2]*y[2];

        // yTx (inner product)
        yTx = (y[0]*x[0]) + (y[1]*x[1]) + (y[2]*x[2]);

        // xTBx (inner product)
        Bx[0] = B[0][0]*x[0] + B[0][1]*x[1] + B[0][2]*x[2];
        Bx[1] = B[1][0]*x[0] + B[1][1]*x[1] + B[1][2]*x[2];
        Bx[2] = B[2][0]*x[0] + B[2][1]*x[1] + B[2][2]*x[2];
        xTBx = (x[0]*Bx[0]) + (x[1]*Bx[1]) + (x[2]*Bx[2]);

        // BxBxT (outer product)
        BxBxT[0][0] = Bx[0]*Bx[0]; BxBxT[0][1] = Bx[0]*Bx[1]; BxBxT[0][2] = Bx[0]*Bx[2];
        BxBxT[1][0] = Bx[1]*Bx[0]; BxBxT[1][1] = Bx[1]*Bx[1]; BxBxT[1][2] = Bx[1]*Bx[2];
        BxBxT[2][0] = Bx[2]*Bx[0]; BxBxT[2][1] = Bx[2]*Bx[1]; BxBxT[2][2] = Bx[2]*Bx[2];

        // construct everything
        B_k1[0][0] = B[0][0] + yyT[0][0]/yTx - BxBxT[0][0]/xTBx; 
        B_k1[0][1] = B[0][1] + yyT[0][1]/yTx - BxBxT[0][1]/xTBx;
        B_k1[0][2] = B[0][2] + yyT[0][2]/yTx - BxBxT[0][2]/xTBx;

        B_k1[1][0] = B[1][0] + yyT[1][0]/yTx - BxBxT[1][0]/xTBx; 
        B_k1[1][1] = B[1][1] + yyT[1][1]/yTx - BxBxT[1][1]/xTBx;
        B_k1[1][2] = B[1][2] + yyT[1][2]/yTx - BxBxT[1][2]/xTBx;

        B_k1[2][0] = B[2][0] + yyT[2][0]/yTx - BxBxT[2][0]/xTBx; 
        B_k1[2][1] = B[2][1] + yyT[2][1]/yTx - BxBxT[2][1]/xTBx;
        B_k1[2][2] = B[2][2] + yyT[2][2]/yTx - BxBxT[2][2]/xTBx;
      }
      return B_k1; // this is the new nu
    }

    Matrix<Double> invEBHysteresis::EvaluateLocalMuGBM(Vector<Double> dH, Vector<Double> dB, UInt idx){
      Matrix<Double> nu;
      nu.Resize(dim_,dim_);

      Vector<Double> dD_minus_epsilon_km1_times_dE(dim_);
      Vector<Double> dHvec(dim_);
      Matrix<Double> rightM(dim_,dim_);
      double offset = 1e-13;

      if(dim_ == 2){
        dH[0] = dH[0] + offset; dH[1] = dH[1] + offset;
        dB[0] = dB[0] + offset; dB[1] = dB[1] + offset;
        dD_minus_epsilon_km1_times_dE[0] = dB[0] - (nu_[idx][0][0]*dH[0] + nu_[idx][0][1]*dH[1]);
        dD_minus_epsilon_km1_times_dE[1] = dB[1] - (nu_[idx][1][0]*dH[0] + nu_[idx][1][1]*dH[1]);

        dHvec[0] = dH[0];
        dHvec[1] = dH[1];
        dD_minus_epsilon_km1_times_dE.ScalarDiv(dHvec.NormL2_squared());
        
        rightM.DyadicMult(dD_minus_epsilon_km1_times_dE, dHvec);
        
        nu[0][0] = nu_[idx][0][0] + rightM[0][0];
        nu[1][0] = nu_[idx][1][0] + rightM[1][0];
        nu[0][1] = nu_[idx][0][1] + rightM[0][1];
        nu[1][1] = nu_[idx][1][1] + rightM[1][1];
                
        LOG_DBG3(eb)<< "\n\t dD_minus_epsilon_km1_times_dE = " << dD_minus_epsilon_km1_times_dE.ToString()
                    << "\n\t dHvec = " << dHvec.ToString()
                    << "\n\t rightM = " << rightM.ToString();

        for (UInt i = 0; i < dim_; i++){
            if (std::isinf(nu[i][i]) || std::isnan(nu[i][i])){
              Matrix<Double> e = nu_[idx];
                nu[i][i] = e[i][i]; //mu_0; //e[i][i]; //mu_0;
            } 
        }
      }else{ // 3D version
        dD_minus_epsilon_km1_times_dE[0] = dB[0] - (nu_[idx][0][0]*dH[0] + nu_[idx][0][1]*dH[1] + nu_[idx][0][2]*dH[2]);
        dD_minus_epsilon_km1_times_dE[1] = dB[1] - (nu_[idx][1][0]*dH[0] + nu_[idx][1][1]*dH[1] + nu_[idx][1][2]*dH[2]);
        dD_minus_epsilon_km1_times_dE[2] = dB[2] - (nu_[idx][2][0]*dH[0] + nu_[idx][2][1]*dH[1] + nu_[idx][2][2]*dH[2]);

        dHvec[0] = dH[0];
        dHvec[1] = dH[1];
        dHvec[2] = dH[2];
        dD_minus_epsilon_km1_times_dE.ScalarDiv(dHvec.NormL2_squared());
        
        rightM.DyadicMult(dD_minus_epsilon_km1_times_dE, dHvec);
        
        nu[1][0] = nu_[idx][1][0] + rightM[1][0];
        nu[0][1] = nu_[idx][0][1] + rightM[0][1];
        nu[0][2] = nu_[idx][0][2] + rightM[0][2];
        nu[2][0] = nu_[idx][2][0] + rightM[2][0];
        nu[1][2] = nu_[idx][1][2] + rightM[1][2];
        nu[2][1] = nu_[idx][2][1] + rightM[2][1];
        nu[0][0] = nu_[idx][0][0] + rightM[0][0];
        nu[1][1] = nu_[idx][1][1] + rightM[1][1];
        nu[2][2] = nu_[idx][2][2] + rightM[2][2];

        LOG_DBG3(eb)<< "\n\t dD_minus_epsilon_km1_times_dE = " << dD_minus_epsilon_km1_times_dE.ToString()
                    << "\n\t dHvec = " << dHvec.ToString()
                    << "\n\t rightM = " << rightM.ToString();

        nu[0][1] = 0.0;
        nu[1][0] = 0.0;
        nu[2][0] = 0.0;
        nu[0][2] = 0.0;
        nu[1][2] = 0.0;
        nu[2][1] = 0.0;
        for (UInt i = 0; i < dim_; i++){
            if (std::isinf(nu[i][i]) || std::isnan(nu[i][i])){
              Matrix<Double> e = nu_[idx];
                nu[i][i] = e[i][i]; //e[i][i]; //mu_0;
            } 
        }
      }
      return nu;  
    }




    Vector<Double> invEBHysteresis::Evaluate(Vector<Double> B_n, bool saveTmpStageVecs, UInt idx) {
      // define needed variables
      Vector<Double>     ret;
      StdVector<Double>  weight(numS_);
      UInt               num_pin;
      StdVector<Double>  chi(numS_);
      Matrix<Double>     J_k(2,numS_); // dim: (2,num_pin)
      Matrix<Double>     J_k_prev(2,numS_); // previous J_k, dim: (2,num_pin)
      Vector<Double>     gradient(2*numS_);
      Matrix<Double>     hessian(2*numS_, 2*numS_), negative_hessian(2*numS_, 2*numS_);
      Matrix<Double>     delta_J_k(2,numS_);
      Matrix<Double>     eta(2,numS_);
      Double             phi_der0, phi_der00, tolerance_newton, eps_newton; // line search / Newton parameter
      UInt               max_newton_iter;  // Newton parameter
      Vector<Double>     solution_linear_system(2*numS_);
      StdVector<Double>  J_x_k_sol, J_y_k_sol;
      Vector<Double>     H_out(2), J_out(2);

      // init. variables
      weight.Init(1.0/numS_);  // hysteresis parameter
      num_pin = numS_;         // hysteresis parameter
      StdVector<Double>& J_x_k_prev = J_x_k_n_[idx]; // get J_k from previous time step
      StdVector<Double>& J_y_k_prev = J_y_k_n_[idx]; 
      for(int kdx = 1; kdx < num_pin; kdx++){ 
        J_k_prev[0][kdx] = J_x_k_prev[kdx];
        J_k_prev[1][kdx] = J_y_k_prev[kdx];
        J_k[0][kdx] = 0.0;
        J_k[1][kdx] = 0.0;
      }
      J_x_k_sol.Resize(numS_, 0.0);
      J_y_k_sol.Resize(numS_, 0.0);
      max_newton_iter = 100; // some Newton cofigurations (should be accessible from the XML-File??)
      tolerance_newton = 1e-12;
      eps_newton = 1e-14;

      // apply Newton method
      for(int ndx = 0; ndx < max_newton_iter; ndx++){ // Newton method iterations
        // calculate the gradient of energy functional
        gradient = CalcGradientMagEnergy(B_n, J_k, J_k_prev);
        // calculate the hessian of energy functional
        hessian = CalcHessianMagEnergy(J_k, J_k_prev);
        // change sign of hessian matrix
        for(int idx = 0; idx < 2*num_pin; idx++){
          for(int jdx = 0; jdx < 2*num_pin; jdx++){
            negative_hessian[idx][jdx] = -hessian[idx][jdx];
          }
        }
        // solve linear system
        solution_linear_system = LUSolve(negative_hessian,gradient);
        for(int j = 1; j <= num_pin; j++){ // arrange solution_linear_system to have a (2,num_pin) matrix
            delta_J_k[0][j-1] = solution_linear_system[2*j-1-1];
            delta_J_k[1][j-1] = solution_linear_system[2*j-1];
        }
        // apply line search
        phi_der0 = 0;
        for(int j = 0; j < num_pin; j++){ 
          phi_der0 = phi_der0 + solution_linear_system[j]*gradient[j];
        }
        eta = LineSearchArmijoInvEB(phi_der0, B_n, J_k, delta_J_k, J_k_prev);
        // update solution J_k
        for(int kdx = 0; kdx < num_pin; kdx++){ 
          J_k[0][kdx] = J_k[0][kdx] + eta[0][kdx]*delta_J_k[0][kdx];
          J_k[1][kdx] = J_k[1][kdx] + eta[1][kdx]*delta_J_k[1][kdx];
        }
        // check stopping criterion
        if (ndx == 0) {
          phi_der00 = phi_der0;
        }
        if ( std::abs(phi_der0) < (tolerance_newton*std::abs(phi_der00) + eps_newton) ) {
          break;
        }
      } // for loop end
      // store J_k for next time step
      if(saveTmpStageVecs){
        for(int kdx = 0; kdx < num_pin; kdx++) {
          J_x_k_sol[kdx] = J_k[0][kdx];
          J_y_k_sol[kdx] = J_k[1][kdx];
        }
        J_x_k_n_tmp_[idx] = J_x_k_sol;
        J_y_k_n_tmp_[idx] = J_y_k_sol;
      }

      // compute total J
      for(int kdx = 0; kdx < num_pin; kdx++) {
        J_out[0] = J_out[0] + weight[kdx]*J_k[0][kdx];
        J_out[1] = J_out[1] + weight[kdx]*J_k[1][kdx];
      }

      // compute H by general constitutive law
      H_out[0] = (1/mu_0)*(B_n[0] - J_out[0]);
      H_out[1] = (1/mu_0)*(B_n[1] - J_out[1]);

      //!! MAKE IT JUST LINEAR (FOR TESTING)!!
      /* H_out[0] = (1/(mu_0*1))*B_n[0];
      H_out[1] = (1/(mu_0*1))*B_n[1]; */

      //!! MAKE IT JUST NONLINEAR (FOR TESTING)!!
      /* Double norm_B_n;
      Double p_0, p_1, p_2; // parameter
      p_0 = 155.9; p_1 = 1.5; p_2 = 9.0;
      norm_B_n = std::sqrt(std::pow(B_n[0],2) + std::pow(B_n[1],2));
      H_out[0] = (p_0 + (p_1*std::pow(norm_B_n,2*p_2)))*B_n[0];
      H_out[1] = (p_0 + (p_1*std::pow(norm_B_n,2*p_2)))*B_n[1]; */

      // return value
      ret.Push_back(H_out[0]);
      ret.Push_back(H_out[1]);
      return ret;

    } 

    Vector<Double> invEBHysteresis::CalcGradientMagEnergy(Vector<Double> B_n, Matrix<Double> J_k, Matrix<Double> J_k_prev) {
    
      // define needed variables
      Vector<Double>  weight(numS_);
      UInt            num_pin;
      Vector<Double>  chi(numS_);
      Vector<Double>  gradient(2*numS_);
      Vector<Double>  gradient_vacuum(2), gradient_anhyst(2), gradient_dissipation(2);
      Vector<Double>  gradient_internal_energy(2);
      Vector<Double>  J_sum(2);
      Vector<Double>  J(2);
      Double delta_J_x, delta_J_y,norm_delta_J, J_x, J_y;
      double epsilon;

      // init. variables
      weight.Init(1.0/numS_);
      num_pin = numS_;
      for(UInt i = 0; i < numS_; i++){ // hysteresis parameter
        chi[i] = chi_factor_ * i / (numS_-1);
      }
      epsilon = 1e-8;

      // pre-computations for the gradient
      J_sum[0] = 0.0; J_sum[1] = 0.0;
      for(int kdx = 0; kdx < num_pin; kdx++){
        J_sum[0] = J_sum[0] + (weight[kdx]*J_k[0][kdx]);
        J_sum[1] = J_sum[1] + (weight[kdx]*J_k[1][kdx]);
      }
      

      for(int kdx = 1; kdx <= num_pin; kdx++){
        // get kth J from J_k
        J_x = J_k[0][kdx-1]; J_y = J_k[1][kdx-1]; J[0] = J_x; J[1] = J_y;
        // calc. differences
        delta_J_x = J_k[0][kdx-1] - J_k_prev[0][kdx-1]; delta_J_y = J_k[1][kdx-1] - J_k_prev[1][kdx-1];
        norm_delta_J = std::sqrt(std::pow(delta_J_x,2) + std::pow(delta_J_y,2));
        // calc. all gradient parts
        // vacuum
        gradient_vacuum[0] = -(weight[kdx-1]/mu_0)*(B_n[0] - J_sum[0]);
        gradient_vacuum[1] = -(weight[kdx-1]/mu_0)*(B_n[1] - J_sum[1]);
        // anhysteresis
        gradient_internal_energy = CalcGradientInternalEnergy(J);
        gradient_anhyst[0] = weight[kdx-1]*gradient_internal_energy[0];
        gradient_anhyst[1] = weight[kdx-1]*gradient_internal_energy[1];
        // dissipation
        gradient_dissipation[0] = weight[kdx-1]*chi[kdx-1]*delta_J_x/std::sqrt(std::pow(norm_delta_J,2) + epsilon);
        gradient_dissipation[1] = weight[kdx-1]*chi[kdx-1]*delta_J_y/std::sqrt(std::pow(norm_delta_J,2) + epsilon);
        gradient[2*kdx-1-1] = gradient_vacuum[0] + gradient_anhyst[0] + gradient_dissipation[0];
        gradient[2*kdx-1] = gradient_vacuum[1] + gradient_anhyst[1] + gradient_dissipation[1];
      }
      return gradient;
    }

    Vector<Double> invEBHysteresis::CalcGradientInternalEnergy(Vector<Double> J) {
    // define needed variables
    Double norm_J, gradient_internal_energy_magnitude;
    Vector<Double>  gradient_internal_energy_direction(2), gradient_internal_energy(2);

    // parameter (anhysteretic)
    double Js = Js_;
    double epsilon = 1e-8;
    double A = A_;

    // calc. magnitude
    norm_J = std::sqrt(std::pow((std::sqrt(std::pow(J[0],2)+std::pow(J[1],2))),2) + epsilon);
    gradient_internal_energy_magnitude = A * std::tan((M_PI/2) * (norm_J / Js));
    
    // calc. direction
    gradient_internal_energy_direction[0] = J[0]/norm_J;
    gradient_internal_energy_direction[1] = J[1]/norm_J;

    // put everything together
    gradient_internal_energy[0] = gradient_internal_energy_magnitude*gradient_internal_energy_direction[0];
    gradient_internal_energy[1] = gradient_internal_energy_magnitude*gradient_internal_energy_direction[1];

    // return
    return gradient_internal_energy;
    }

    Matrix<Double> invEBHysteresis::CalcHessianMagEnergy(Matrix<Double> J_k, Matrix<Double> J_k_prev) {

      // define needed variables
      StdVector<Double>  weight(numS_), chi_k(numS_);
      UInt               num_pin;
      Matrix<Double>     hessian(2*numS_, 2*numS_);
      Matrix<Double>     hessian_k(2, 2);
      Matrix<Double>     hessian_internal_energy(2, 2);
      Matrix<Double>     hessian_dissipation_energy(2, 2);
      Double             J_x, J_y, delta_J_x, delta_J_y, norm_delta_J;
      Vector<Double>     J(2), delta_J(2);
      UInt               chi;

      // init. variables
      weight.Init(1.0/numS_);
      num_pin = numS_;
      for(UInt i = 0; i < numS_; i++){ // hysteresis parameter
        chi_k[i] = chi_factor_ * i / (numS_-1);
      }

      for(int kdx = 1; kdx <= num_pin; kdx++){
        J_x = J_k[0][kdx-1]; J_y = J_k[1][kdx-1]; J[0] = J_x; J[1] = J_y; // get kth J from J_k
        chi = chi_k[kdx-1]; // get kth chi from chi_k
        // calc. differences
        delta_J_x = J_k[0][kdx-1] - J_k_prev[0][kdx-1]; delta_J_y = J_k[1][kdx-1] - J_k_prev[1][kdx-1];
        delta_J[0] = delta_J_x; delta_J[1] = delta_J_y;
        norm_delta_J = std::sqrt(std::pow(delta_J_x,2) + std::pow(delta_J_y,2));

        // vacuum energy
        for(int idx = 1; idx <= num_pin; idx++){
          hessian[2*kdx-1-1][2*idx-1-1] = hessian[2*kdx-1-1][2*idx-1-1] + (weight[idx-1]*weight[kdx-1]/mu_0);
          hessian[2*kdx-1][2*idx-1] = hessian[2*kdx-1][2*idx-1] + (weight[idx-1]*weight[kdx-1]/mu_0);
        }

        // internal energy
        hessian_internal_energy = CalcHessianInternalEnergy(J);
        hessian_k[0][0] = weight[kdx-1]*hessian_internal_energy[0][0];
        hessian_k[1][1] = weight[kdx-1]*hessian_internal_energy[1][1];
        hessian_k[0][1] = weight[kdx-1]*hessian_internal_energy[0][1];
        hessian_k[1][0] = weight[kdx-1]*hessian_internal_energy[1][0];

        // dissipative energy
        hessian_dissipation_energy = CalcHessianDissipationEnergy(delta_J, norm_delta_J, chi);
        hessian_k[0][0] = hessian_k[0][0] + weight[kdx-1]*hessian_dissipation_energy[0][0];
        hessian_k[1][1] = hessian_k[1][1] + weight[kdx-1]*hessian_dissipation_energy[1][1];
        hessian_k[0][1] = hessian_k[0][1] + weight[kdx-1]*hessian_dissipation_energy[0][1];
        hessian_k[1][0] = hessian_k[1][0] + weight[kdx-1]*hessian_dissipation_energy[1][0];

        // put everything together
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                hessian[2*kdx - 2 + i][2*kdx - 2 + j] += hessian_k[i][j];
            }
        }
      }
      // return
      return hessian;
    }

    Matrix<Double> invEBHysteresis::CalcHessianInternalEnergy(Vector<Double> J) {
      // define needed variables
      Matrix<Double> hessian(2, 2);
      Matrix<Double> hessian1(2, 2);
      Matrix<Double> hessian2(2, 2); 
      Double         norm_J, factor1, factor2, sec_factor; 
      
      // parameter (anhysteretic)
      double Js = Js_;
      double epsilon = 1e-8;
      double A = A_;

      // pre-computations for hessian
      norm_J = std::sqrt(std::pow((std::sqrt(std::pow(J[0],2)+std::pow(J[1],2))),2) + epsilon);
      factor1 = A * std::tan((M_PI/2) * (norm_J / Js));
      sec_factor = 1/std::cos(norm_J * M_PI/ 2/ Js);
      factor2 = (A*M_PI/2/Js)*std::pow(sec_factor,2); 

      // part 1
      hessian1[0][0] =  factor1*(1/norm_J - J[0]*J[0]/(norm_J*norm_J*norm_J));
      hessian1[1][1] =  factor1*(1/norm_J - J[1]*J[1]/(norm_J*norm_J*norm_J));
      hessian1[0][1] = -factor1*J[0]*J[1]/(norm_J*norm_J*norm_J);
      hessian1[1][0] = -factor1*J[1]*J[0]/(norm_J*norm_J*norm_J);

      // part 2
      hessian2[0][0] = factor2*J[0]*J[0]/(norm_J*norm_J);
      hessian2[1][1] = factor2*J[1]*J[1]/(norm_J*norm_J);
      hessian2[0][1] = factor2*J[0]*J[1]/(norm_J*norm_J);
      hessian2[1][0] = factor2*J[1]*J[0]/(norm_J*norm_J);

      // sum everything up
      hessian[0][0] = hessian1[0][0] + hessian2[0][0];
      hessian[1][1] = hessian1[1][1] + hessian2[1][1];
      hessian[0][1] = hessian1[0][1] + hessian2[0][1];
      hessian[1][0] = hessian1[1][0] + hessian2[1][0];

      // return
      return hessian;
    }

    Matrix<Double> invEBHysteresis::CalcHessianDissipationEnergy(Vector<Double> delta_J, Double norm_delta_J, UInt chi) {
      // define needed variables
      Matrix<Double> hessian(2, 2);
      Matrix<Double> hessian1(2, 2);
      Matrix<Double> hessian2(2, 2); 
      Double         factor1, factor2; 

      // parameter (dissipation)
      double epsilon = 1e-8;

      // part 1 of dissipation hessian
      factor1 = chi*(1/std::sqrt(std::pow(norm_delta_J,2) + epsilon));
      hessian1[0][0] = factor1; hessian1[1][1] = factor1; 

      // part 2 of hessian
      factor2 = chi/(std::pow((std::pow(norm_delta_J,2) + epsilon),3/2));
      hessian2[0][0] = -factor2*delta_J[0]*delta_J[0];
      hessian2[1][1] = -factor2*delta_J[1]*delta_J[1];
      hessian2[0][1] = -factor2*delta_J[0]*delta_J[1];
      hessian2[1][0] = -factor2*delta_J[1]*delta_J[0];

      // sum everything up
      hessian[0][0] = hessian1[0][0] + hessian2[0][0];
      hessian[1][1] = hessian1[1][1] + hessian2[1][1];
      hessian[0][1] = hessian1[0][1] + hessian2[0][1];
      hessian[1][0] = hessian1[1][0] + hessian2[1][0];

      // return
      return hessian;
    }

    Matrix<Double> invEBHysteresis::LineSearchArmijoInvEB(Double phi_der0, Vector<Double> B_n, Matrix<Double> J_k, Matrix<Double> delta_J_k, Matrix<Double> J_k_prev) {

      // define needed variables
      double         delta, gamma, max_line_search_iter; // line search parameter
      Matrix<Double> eta0(2,numS_);
      Matrix<Double> eta(2,numS_);
      double         phi0, energy; // energies
      Matrix<Double> J_k_update(2,numS_);
      UInt           num_pin;

      // init. variables
      delta = 0.5; 
      gamma = 0.1;
      num_pin = numS_;
      max_line_search_iter = 1000;
      for(int kdx = 0; kdx < num_pin; kdx++){
         eta0[0][kdx] = 0; eta0[1][kdx] = 0;
         eta[0][kdx] = 1; eta[1][kdx] = 1;
      }


      // compute starting energies
      for(int kdx = 0; kdx < num_pin; kdx++){
        J_k_update[0][kdx] = J_k[0][kdx] + eta0[0][kdx]*delta_J_k[0][kdx];
        J_k_update[1][kdx] = J_k[1][kdx] + eta0[1][kdx]*delta_J_k[1][kdx];
      }
      phi0 = CalcMagEnergy(B_n, J_k_update, J_k_prev);
      for(int kdx = 0; kdx < num_pin; kdx++){
        J_k_update[0][kdx] = J_k[0][kdx] + eta[0][kdx]*delta_J_k[0][kdx];
        J_k_update[1][kdx] = J_k[1][kdx] + eta[1][kdx]*delta_J_k[1][kdx];
      }
      energy = CalcMagEnergy(B_n, J_k_update, J_k_prev);

      // do line search
      for(int idx = 0; idx < max_line_search_iter; idx++){
        if (energy <= (phi0 + gamma*eta[0][0]*phi_der0 + 1e-12*phi0)) {
          break;
        } else {
          for(int kdx = 0; kdx < num_pin; kdx++){ // update to new eta and compute J_k_update for that eta
            eta[0][kdx] = eta[0][kdx]*delta;
            eta[1][kdx] = eta[1][kdx]*delta;
            J_k_update[0][kdx] = J_k[0][kdx] + eta[0][kdx]*delta_J_k[0][kdx];
            J_k_update[1][kdx] = J_k[1][kdx] + eta[1][kdx]*delta_J_k[1][kdx];
          }
          energy = CalcMagEnergy(B_n, J_k_update, J_k_prev);
        }
      }

      // return 
      return eta;
    }

    Double invEBHysteresis::CalcMagEnergy(Vector<Double> B, Matrix<Double> J_k_update, Matrix<Double> J_k_prev){
      // define needed variables
      Double epsilon;
      StdVector<Double>  weight(numS_), chi_k(numS_);;
      UInt               num_pin;
      Vector<Double>     J_sum(2);
      Vector<Double>     muH(2);
      Matrix<Double>     diff_J_k(2,numS_);
      Vector<Double>     norm_diff_J_k(numS_);
      Double             energy_vacuum;
      Vector<Double>     energy_internal(numS_), energy_dissipation(numS_);
      Double             magnetic_energy;
      Double             sum_weighted_energy_parts;

      // init. variables
      epsilon = 1e-8;
      weight.Init(1.0/numS_);
      num_pin = numS_;
      J_sum[0] = 0.0; J_sum[1] = 0.0;
      for(UInt i = 0; i < numS_; i++){ // hysteresis parameter
        chi_k[i] = chi_factor_ * i / (numS_-1) + 1e-22;
        chi_k[i] = chi_factor_ * i / (numS_-1) + 0.001;
        chi_k[i] = chi_factor_ * i / (numS_-1) + 1e-22;
      }

      // pre-computations for energies
      for(int kdx = 0; kdx < num_pin; kdx++){
        J_sum[0] = J_sum[0] + (weight[kdx]*J_k_update[0][kdx]);
        J_sum[1] = J_sum[1] + (weight[kdx]*J_k_update[1][kdx]);
        diff_J_k[0][kdx] = J_k_update[0][kdx] - J_k_prev[0][kdx];
        diff_J_k[1][kdx] = J_k_update[1][kdx] - J_k_prev[1][kdx];
        norm_diff_J_k[kdx] = std::sqrt(std::pow(diff_J_k[0][kdx],2) + std::pow(diff_J_k[1][kdx],2));
      }
      muH[0] = B[0] - J_sum[0];
      muH[1] = B[1] - J_sum[1];

      // compute energies
      energy_vacuum = 1/(2*mu_0)*std::pow(std::sqrt(std::pow(muH[0],2) + std::pow(muH[1],2)),2);
      energy_internal = CalcInternalEnergy(J_k_update);
      for(int kdx = 0; kdx < num_pin; kdx++){
        energy_dissipation[kdx] = chi_k[kdx]*std::sqrt(std::pow(norm_diff_J_k[kdx],2) + epsilon);
      }

      // sum everything up
      sum_weighted_energy_parts = 0.0; // consists of weight*(energy_internal + energy_dissipation);
      for(int kdx = 0; kdx < num_pin; kdx++){
        sum_weighted_energy_parts = sum_weighted_energy_parts + weight[kdx]*(energy_internal[kdx] + energy_dissipation[kdx]);
      }
      magnetic_energy = energy_vacuum + sum_weighted_energy_parts;

      // return
      return magnetic_energy;

    }

    Vector<Double> invEBHysteresis::CalcInternalEnergy(Matrix<Double> J_k) {
      // define needed variables
      double Js, A;
      Vector<Double>  norm_J(numS_);
      Vector<Double>  internal_energy(numS_);
      UInt            num_pin;

      // init. variables
      Js = Js_;
      A = A_;
      num_pin = numS_;

      // pre-computations
      for(int kdx = 0; kdx < num_pin; kdx++){
        norm_J[kdx] = std::sqrt(std::pow(J_k[0][kdx],2) + std::pow(J_k[1][kdx],2));
      }
      
      // compute internal energy for all k
      for(int kdx = 0; kdx < num_pin; kdx++){
        if (norm_J[kdx] < Js) {
          internal_energy[kdx] = -(2*A*Js/M_PI)*std::log(std::cos( (M_PI/2)*norm_J[kdx] / Js) );
        } else {
          internal_energy[kdx] = HUGE_VAL;
        }
      }

      // return
      return internal_energy;
    }

    Vector<Double> invEBHysteresis::LUSolve(Matrix<Double> A, Vector<Double> b){
      // define needed variables 
      int                    number_rows_A = 2*numS_;
      Vector<Double>         y, solution_linear_system(number_rows_A);
      Matrix<Double>         L(number_rows_A,number_rows_A), U(number_rows_A,number_rows_A);

      // init. variables 
      // (see definitions)

      // solve the system by LU Decomposition
      LUDecompose(A, L, U); // Step 1: apply LU Decomposition
      y = LUForwardSubstitution(L, b); // Step 2: Solve L * y = b using forward substitution
      solution_linear_system = LUBackwardSubstitution(U, y); // Step 3: Solve U * x = y using backward substitution
      
      // return
      return solution_linear_system;
    }

    bool invEBHysteresis::LUDecompose(Matrix<Double> A, Matrix<Double>& L, Matrix<Double>& U) {
      // define needed variables
      UInt number_rows_A;

      // init. variables
      number_rows_A = 2*numS_;

      // LU Decomposition
      for (int idx = 0; idx < number_rows_A; ++idx) {
        for (int jdx = idx; jdx < number_rows_A; ++jdx) { // Upper triangular matrix U
          U[idx][jdx] = A[idx][jdx];
          for (int kdx = 0; kdx < idx; ++kdx) {
              U[idx][jdx] -= L[idx][kdx] * U[kdx][jdx];
          }
        }
        for (int jdx = idx; jdx < number_rows_A; ++jdx) { // Lower triangular matrix L
          if (idx == jdx) {
              L[idx][idx] = 1; // Diagonal elements of L are set to 1
          } else {
              L[jdx][idx] = A[jdx][idx];
              for (int kdx = 0; kdx < idx; ++kdx) {
                  L[jdx][idx] -= L[jdx][kdx] * U[kdx][idx];
              }
              L[jdx][idx] /= U[idx][idx];
          }
        }
      }
      return true;
    }

    Vector<Double> invEBHysteresis::LUForwardSubstitution(const Matrix<Double>& L, const Vector<Double>& b) {
      // define needed variables
      UInt number_rows_L = 2*numS_;
      Vector<Double> y(number_rows_L);

      // init. variables
      // (see definitions)

      // apply forward substitution
      for (int idx = 0; idx < number_rows_L; ++idx) {
          y[idx] = b[idx];
          for (int jdx = 0; jdx < idx; ++jdx) {
              y[idx] -= L[idx][jdx] * y[jdx];
          }
          y[idx] /= L[idx][idx];  // L[idx][idx] is 1 by definition
      }

      // return
      return y;
    }

    Vector<Double> invEBHysteresis::LUBackwardSubstitution(const Matrix<Double>& U, const Vector<double>& y) {
      // define needed variables
      UInt number_rows_U = 2*numS_;
      Vector<Double> x(number_rows_U);

      // init. variables
      // (see definitions)

      // apply backward substitution
      for (int idx = number_rows_U - 1; idx >= 0; --idx) {
          x[idx] = y[idx];
          for (int jdx = idx + 1; jdx < number_rows_U; ++jdx) {
              x[idx] -= U[idx][jdx] * x[jdx];
          }
          x[idx] /= U[idx][idx];
      }

      // return
      return x;
    }


      







  }
