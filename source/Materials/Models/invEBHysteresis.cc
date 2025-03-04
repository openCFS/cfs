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
      mu_0 = 1.256637061e-06;

      Js_ = ParameterMap["Js"];
      A_ = ParameterMap["A"];
      numS_ = ParameterMap["numS"];
      chi_factor_ = ParameterMap["chi_factor"];
      jacobian_method_ = ParameterMap["jacobian_method"];

      B0_.Resize(numElems_, StdVector<Double>(dim_));
      B1_.Resize(numElems_, StdVector<Double>(dim_));
      J0_.Resize(numElems_, StdVector<Double>(dim_));
      J1_.Resize(numElems_, StdVector<Double>(dim_));

      J_x_k_prev_.Resize(numElems_, StdVector<Double>(numS_));
      J_y_k_prev_.Resize(numElems_, StdVector<Double>(numS_));
      J_z_k_prev_.Resize(numElems_, StdVector<Double>(numS_));

      J_x_k_n_.Resize(numElems_, StdVector<Double>(numS_));
      J_y_k_n_.Resize(numElems_, StdVector<Double>(numS_));
      J_z_k_n_.Resize(numElems_, StdVector<Double>(numS_));

      J_x_k_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));
      J_x_k_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));
      J_x_k_n_tmp_.Resize(numElems_, StdVector<Double>(numS_));

      nu_.Resize(numElems_, Matrix<Double>(dim_,dim_));

      isFirstTime_.Resize(numElems_);
      isFirstTime_.Init(1);

      isFirstTimeFinished_ = false;

      timeStep_ = 1;

      mp_ = domain->GetMathParser();
      globalIter_ = 0;

      hasElemSolution_.Resize(numElems_, false);
    }

    Vector<Double> invEBHysteresis::GetFieldIntensity(Vector<Double> BVec, const Integer ElemNum) {
      // define needed variables
      Vector<Double> H(dim_);
      UInt idx = ElemNum2Idx_[ElemNum];

      // compute H
      H = Evaluate(BVec, false, idx);

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
      /* Vector<Double>     J;
      Vector<Double>     J_dummy; */
      Matrix<Double>     nu(dim_,dim_); nu.InitValue(0.0);
      /* StdVector<Double>& B0 = B0_[idx];
      StdVector<Double>& J0 = J0_[idx];
      StdVector<Double>  H_k(dim_);
      StdVector<Double>  H_k_0(dim_);
      StdVector<Double>  delta_H(dim_);
      StdVector<Double>  delta_B(dim_); */

      
      /* // algorithm to determine nu = dH/dB;
      if(hasElemSolution_[idx] == true){ // nu of this element has already been computed
        return nu_[idx];
      }
      if(timeStep_== 1 && globalIter_ == 1){ // starting value at nu = diag(1/mu0);
        for(UInt i = 0; i < dim_; i++) {
          for(UInt j = 0; j < dim_; j++) {
            if (i == j) {
              nu[i][j] = 1/mu_0;
            } else {
              nu[i][j] = 0;
            }
          }
        }
        nu_[idx] = nu;
        hasElemSolution_[idx] = true;
        return nu;
      }
      J = Evaluate(BVec, false, idx); // get current J
      for(UInt i = 0; i < dim_; i++) { // get the values for H and M from the last timestep and get deltaH and deltaB
        H_k[i] = (1/mu_0) * (BVec[i] - J[i]);
        H_k_0[i] = (1/mu_0) * (B0[i] - J0[i]);
        delta_B[i] = BVec[i] - B0[i];
        delta_H[i] = H_k[i] - H_k_0[i];
        B1_[idx][i] = BVec[i];
        J1_[idx][i] = J[i];
      }
      nu = EvaluateLocalNuBFGS(delta_H, delta_B, idx); // determine nu
      J_dummy = Evaluate(BVec, true, idx); // just to update all J_k values
      nu_[idx] = nu; // store all quantities
      J0_[idx] = J1_[idx];
      B0_[idx] = B1_[idx]; */

      for(UInt i = 0; i < dim_; i++) {
        for(UInt j = 0; j < dim_; j++) {
          if (i == j) {
            nu[i][j] = 1/mu_0;
          } else {
            nu[i][j] = 0;
          }
        }
      }
      nu_[idx] = nu;
      hasElemSolution_[idx] = true; // mark this element as computed

      // return
      return nu;
    }


    Matrix<Double> invEBHysteresis::EvaluateLocalNuBFGS(StdVector<Double> dH, StdVector<Double> dB, UInt idx){
      // define needed variables
      Matrix<Double> yyT(dim_, dim_);
      Matrix<Double> BxBxT(dim_, dim_);
      Matrix<Double> B_k1(dim_, dim_);
      Matrix<Double> B = nu_[idx];
      Double yTx;
      Double xTBx;
      StdVector<Double> y = dH;
      StdVector<Double> x = dB;
      StdVector<Double> Bx;

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
        Bx[2] = B[2][0]*x[0] + B[2][1]*x[1];
        xTBx = (x[0]*Bx[0]) + (x[1]*Bx[1]);

        // BxBxT (outer product)
        BxBxT[0][0] = Bx[0]*Bx[0]; BxBxT[0][1] = Bx[0]*Bx[1];
        BxBxT[1][0] = Bx[1]*Bx[0]; BxBxT[1][1] = Bx[1]*Bx[1];
        
        // construct everything
        B_k1[0][0] = B[0][0] + yyT[0][0]/yTx - BxBxT[0][0]/xTBx; 
        B_k1[0][1] = B[0][1] + yyT[0][1]/yTx - BxBxT[0][1]/xTBx;

        B_k1[1][0] = B[1][0] + yyT[1][0]/yTx - BxBxT[1][0]/xTBx; 
        B_k1[1][1] = B[1][1] + yyT[1][1]/yTx - BxBxT[1][1]/xTBx;
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




    Vector<Double> invEBHysteresis::Evaluate(Vector<Double> B_n, bool saveTmpStageVecs, UInt idx) {
      // define needed variables
      Vector<Double>     ret;
      /* StdVector<Double>  weight(numS_);
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
      StdVector<Double>  J_x_k_sol, J_y_k_sol; */
      Vector<Double>     H_out(2);//, J_out(2);


      /* // init. variables
      weight.Init(1.0/numS_);  // hysteresis parameter
      num_pin = numS_;         // hysteresis parameter
      StdVector<Double>& J_x_k_prev = J_x_k_n_[idx]; // get J_k from previous time step
      StdVector<Double>& J_y_k_prev = J_y_k_n_[idx]; 
      for(int kdx = 1; kdx < num_pin; kdx++){ 
        J_k_prev[0][kdx] = J_x_k_prev[kdx];
        J_k_prev[1][kdx] = J_y_k_prev[kdx];
        J_k[0][kdx] = 0;
        J_k[1][kdx] = 0;
      }
      J_x_k_sol.Resize(numS_, 0.0);
      J_y_k_sol.Resize(numS_, 0.0);

      // some Newton cofigurations (should be accessible from the XML-File??)
      max_newton_iter = 5000;
      tolerance_newton = 10^-12;
      eps_newton = 10^-14;

      
      // apply Newton method
      for(int ndx = 0; ndx < max_newton_iter; ndx++){
        // calculate the gradient of energy functional
        gradient = CalcGradientMagEnergy(B_n, J_k, J_k_prev);
        // calculate the hessian of energy functional
        hessian = CalcHessianMagEnergy(J_k, J_k_prev);
        // change sign of hessian matrix
        for(int idx = 1; idx < 2*num_pin; idx++){
          for(int jdx = 1; jdx < 2*num_pin; jdx++){
            negative_hessian[idx][jdx] = -hessian[idx][jdx];
          }
        }
        // solve linear system
        solution_linear_system = LUSolve(hessian,gradient);
        for(int j = 1; j < num_pin; j++){ // arrange solution_linear_system to have a (2,num_pin) matrix
            delta_J_k[0][j-1] = solution_linear_system[2*j-1-1];
            delta_J_k[1][j-1] = solution_linear_system[2*j-1];
        }
        // line search
        phi_der0 = 0;
        for(int j = 1; j < num_pin; j++){ 
          phi_der0 = phi_der0 + solution_linear_system[j]*gradient[j];
        }
        eta = LineSearchArmijoInvEB(phi_der0, B_n, J_k, delta_J_k, J_k_prev);
        // update solution J_k
        for(int kdx = 1; kdx < num_pin; kdx++){ 
          J_k[0][kdx] = J_k[0][kdx] + eta[0][kdx]*delta_J_k[0][kdx];
          J_k[1][kdx] = J_k[1][kdx] + eta[1][kdx]*delta_J_k[1][kdx];
        }
        // check stopping criterion
        if (ndx == 1) {
          phi_der00 = phi_der0;
        }
        if ( std::abs(phi_der0) < (tolerance_newton*std::abs(phi_der00) + eps_newton) ) {
          break;
        }
      } // for loop end
      // store J_k for next time step
      if(saveTmpStageVecs){
        for(int kdx = 1; kdx < num_pin; kdx++) {
          J_x_k_sol[kdx] = J_k[0][kdx];
          J_y_k_sol[kdx] = J_k[1][kdx];
        }
        J_x_k_n_tmp_[idx] = J_x_k_sol;
        J_y_k_n_tmp_[idx] = J_y_k_sol;
      }

      // compute total J
      for(int kdx = 1; kdx < num_pin; kdx++) {
        J_out[0] = J_out[0] + weight[kdx]*J_k[0][kdx];
        J_out[1] = J_out[1] + weight[kdx]*J_k[1][kdx];
      }

      // compute H by general constitutive law
      H_out[0] = (1/mu_0)*(B_n[0] - J_out[0]);
      H_out[1] = (1/mu_0)*(B_n[1] - J_out[1]); */

      //!! MAKE IT JUST LINEAR (FOR TESTING)!!
      H_out[0] = (1/(mu_0*1))*B_n[0];
      H_out[1] = (1/(mu_0*1))*B_n[1];

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
      int epsilon;

      // init. variables
      weight.Init(1.0/numS_);
      num_pin = numS_;
      for(UInt i = 0; i < numS_; i++){ // hysteresis parameter
        chi[i] = chi_factor_ * i / (numS_-1) + 1e-22;
        chi[i] = chi_factor_ * i / (numS_-1) + 0.001;
        chi[i] = chi_factor_ * i / (numS_-1) + 1e-22;
      }
      J_sum[0] = 0.0; J_sum[1] = 0.0;
      epsilon = 1e-8;

      // pre-computations for the gradient
      for(int kdx = 0; kdx < num_pin; kdx++){
        J_sum[0] = J_sum[0] + (weight[kdx]*J_k[0][kdx]);
        J_sum[1] = J_sum[1] + (weight[kdx]*J_k[1][kdx]);
      }
      

      for(int kdx = 0; kdx < num_pin; kdx++){
        // get kth J from J_k
        J_x = J_k[0][kdx]; J_y = J_k[1][kdx]; J[0] = J_x; J[1] = J_y;
        // calc. differences
        delta_J_x = J_k[0][kdx] - J_k_prev[0][kdx]; delta_J_y = J_k[1][kdx] - J_k_prev[1][kdx];
        norm_delta_J = std::sqrt(std::pow(delta_J_x,2) + std::pow(delta_J_y,2));
        // calc. all gradient parts
        // vacuum
        gradient_vacuum[0] = -(weight[kdx]/mu_0)*(B_n[0] - J_sum[0]);
        gradient_vacuum[1] = -(weight[kdx]/mu_0)*(B_n[1] - J_sum[1]);
        // anhysteresis
        gradient_internal_energy = CalcGradientInternalEnergy(J);
        gradient_anhyst[0] = weight[kdx]*gradient_internal_energy[0];
        gradient_anhyst[1] = weight[kdx]*gradient_internal_energy[1];
        // dissipation
        gradient_dissipation[0] = weight[kdx]*chi[kdx]*delta_J_x/std::sqrt(std::pow(norm_delta_J,2) + epsilon);
        gradient_dissipation[1] = weight[kdx]*chi[kdx]*delta_J_y/std::sqrt(std::pow(norm_delta_J,2) + epsilon);
        gradient[2*kdx-1] = gradient_vacuum[0] + gradient_anhyst[0] + gradient_dissipation[0];
        gradient[2*kdx] = gradient_vacuum[1] + gradient_anhyst[1] + gradient_dissipation[1];
      }
      return gradient;
    }

    Vector<Double> invEBHysteresis::CalcGradientInternalEnergy(Vector<Double> J) {
    // define needed variables
    Double norm_J, gradient_internal_energy_magnitude;
    Vector<Double>  gradient_internal_energy_direction(2), gradient_internal_energy(2);

    // parameter (anhysteretic)
    int Js = Js_;
    int epsilon = 1e-8;
    int A = A_;

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
      int                row_hessian, column_hessian;
      UInt               chi;

      // init. variables
      weight.Init(1.0/numS_);
      num_pin = numS_;
      for(UInt i = 0; i < numS_; i++){ // hysteresis parameter
        chi_k[i] = chi_factor_ * i / (numS_-1) + 1e-22;
        chi_k[i] = chi_factor_ * i / (numS_-1) + 0.001;
        chi_k[i] = chi_factor_ * i / (numS_-1) + 1e-22;
      }

      for(int kdx = 0; kdx < num_pin; kdx++){
        J_x = J_k[0][kdx]; J_y = J_k[1][kdx]; J[0] = J_x; J[1] = J_y; // get kth J from J_k
        chi = chi_k[kdx]; // get kth chi from chi_k
        // calc. differences
        delta_J_x = J_k[0][kdx] - J_k_prev[0][kdx]; delta_J_y = J_k[1][kdx] - J_k_prev[1][kdx];
        delta_J[0] = delta_J_x; delta_J[1] = delta_J_y;
        norm_delta_J = std::sqrt(std::pow(delta_J_x,2) + std::pow(delta_J_y,2));

        // vacuum energy
        for(int idx = 1; idx < num_pin; idx++){
            hessian[2*kdx-1-1][2*idx-1-1] = hessian[2*kdx-1-1][2*idx-1-1] + (weight[idx-1]*weight[kdx-1]/mu_0);
            hessian[2*kdx-1][2*idx-1] = hessian[2*kdx-1][2*idx-1] + (weight[idx-1]*weight[kdx-1]/mu_0);
        }

        // internal energy
        hessian_internal_energy = CalcHessianInternalEnergy(J);
        hessian_k[0][0] = weight[kdx]*hessian_internal_energy[0][0];
        hessian_k[1][1] = weight[kdx]*hessian_internal_energy[1][1];
        hessian_k[0][1] = weight[kdx]*hessian_internal_energy[0][1];
        hessian_k[1][0] = weight[kdx]*hessian_internal_energy[1][0];

        // dissipative energy
        hessian_dissipation_energy = CalcHessianDissipationEnergy(delta_J, norm_delta_J, chi);
        hessian_k[0][0] = hessian_k[0][0] + weight[kdx]*hessian_dissipation_energy[0][0];
        hessian_k[1][1] = hessian_k[1][1] + weight[kdx]*hessian_dissipation_energy[1][1];
        hessian_k[0][1] = hessian_k[0][1] + weight[kdx]*hessian_dissipation_energy[0][1];
        hessian_k[1][0] = hessian_k[1][0] + weight[kdx]*hessian_dissipation_energy[1][0];

        // put everything together
        row_hessian = 2*kdx;
        column_hessian = 2*kdx;
        for (int i = 1; i < 2; ++i) {
          for (int j = 1; j < 2; ++j) {
              hessian[row_hessian + i - 1][column_hessian + j - 1] += hessian_k[i][j];
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
      int Js = Js_;
      int epsilon = 1e-8;
      int A = A_;

      // pre-computations for hessian
      norm_J = std::sqrt(std::pow((std::sqrt(std::pow(J[0],2)+std::pow(J[1],2))),2) + epsilon);
      factor1 = A * std::tan((M_PI/2) * (norm_J / Js));
      sec_factor = 1/std::cos(norm_J * M_PI/ 2/ Js);
      factor2 = (A*M_PI/2/Js)*std::pow(sec_factor,2); 

      // part 1
      hessian1[0][0] =  factor1*(1./norm_J - J[0]*J[0]/(norm_J*norm_J*norm_J));
      hessian1[1][1] =  factor1*(1./norm_J - J[1]*J[1]/(norm_J*norm_J*norm_J));
      hessian1[0][1] = -factor1*(1./norm_J - J[0]*J[1]/(norm_J*norm_J*norm_J));
      hessian1[1][0] = -factor1*(1./norm_J - J[1]*J[0]/(norm_J*norm_J*norm_J));

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
      int epsilon = 1e-8;

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
      Matrix<Double> J_k_update;
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
      int Js, A;
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
      Matrix<Double>         L, U;

      // init. variables 
      // (see definitions)

      // solve the system by LU Decomposition
      LUDecompose(A, L, U); // apply LU Decomposition
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
