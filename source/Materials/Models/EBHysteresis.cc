#include <iostream>
#include <fstream>
#include <cmath>

#include "DataInOut/Logging/LogConfigurator.hh"

#include "EBHysteresis.hh"
#include "Model.hh"

#include "MatVec/Vector.hh"


#include "Utils/mathParser/mathParser.hh"
#include "Domain/Domain.hh"
#include "Domain/CoefFunction/CoefFunction.hh"


namespace CoupledField {

DEFINE_LOG(eb, "EBHysteresis")


  EBHysteresis::EBHysteresis() : Model(),
                                 numElems_{0},
                                 Ps_{0}, A_{0}, mu0_{0}, numS_{0}, chi_factor_{0}, jacobian_method_{0},
                                 mp_{nullptr}, iterTracker4Mu_{0},
                                 isMH_{false}
  {
  }

  EBHysteresis::~EBHysteresis() 
  {
  }


  void EBHysteresis::Init(std::map<std::string, double> ParameterMap, std::map<std::string, string> StringParameterMap, shared_ptr<ElemList> entityList, UInt dim) 
  {
    
    dim_ = dim;
    mu0_ = 4 * M_PI * 1e-7;

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
    // numS_ = ParameterMap["numS"]; 
    // chi_factor_ = ParameterMap["chi_factor"];
    jacobian_method_ = ParameterMap["jacobian_method"];
    anhyst_type_ = StringParameterMap["anhyst_type"];
    anhyst_formula_ = StringParameterMap["anhyst_formula"];
    //anhyst_type_=1; // hardcoded to test 
    if (anhyst_type_ == "analytic_anhysteresis")
    {
      if(anhyst_formula_ == "atan"){
        Ps_ = ParameterMap["Ps"];
        A_ = ParameterMap["A"];
      }
      if(anhyst_formula_ == "pacejka"){
        m_sat_pacejka_ = ParameterMap["m_sat"];
        a_pacejka_ = ParameterMap["a"];
        b_pacejka_ = ParameterMap["b"];
        c_pacejka_ = ParameterMap["c"];
      }

    }
    else if (anhyst_type_ == "multiscale_anhysteresis")
    {
      // multiscale anhysteresis model
      SMSM_model_ = std::make_unique<SMSM>(ParameterMap["Ps"],
                                           ParameterMap["AS"],
                                           ParameterMap["K1"],
                                           ParameterMap["K2"],
                                           ParameterMap["lambda100"],
                                           ParameterMap["lambda111"],
                                           dim_);
    }
    
    pinning_forces_weight_ = StringParameterMap["weights_file_path"]+"/"+StringParameterMap["pinning_forces_weights_file"];
    std::ifstream file(pinning_forces_weight_);
    double x, y;
    // Check if the file was successfully opened
    if (!file) {
        EXCEPTION("Error: Unable to open file: " << pinning_forces_weight_ << std::endl);
    }

    kappa_file_.Clear();
    omega_file_.Clear();
    while (file >> x >> y) {
        kappa_file_.push_back(x);
        omega_file_.push_back(y);
    }
    // Close the file
    file.close();
    numS_ = kappa_file_.GetSize(); // at this point the numS is overwritten - change later
    isMH_ = ParameterMap["isMH"];
    if(isMH_ == 1.0)
    {
      varHandle_="cacheResult";
    } else {
      varHandle_="step";
    }

    switch (UInt(ParameterMap["approx_type"]))
    {
    case 1:
      approx_type_ = "fullEB";
      break;
    case 2:
      approx_type_ = "approxVPM";
      break;
    }


    Htotal_prev_.Resize(numElems_, StdVector<Double>(dim_));
    Mprev_iter_.Resize(numElems_, StdVector<Double>(dim_));

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

    mu_.Resize(numElems_, Matrix<Double>(dim_, dim_));

    mp_ = domain->GetMathParser();
    iterTracker4Mu_ = 0;
    saveTmpStageVecs_ = false;

    alreadyHasMu_.Resize(numElems_, false);
  }

  Double EBHysteresis::ComputeMaterialParameter(Vector<Double> HVec, const Integer ElemNum)
  {
    // This method gives the rho_art for the h-form.
    Double beta = 1; // default value for n
    // SVD approach
    Matrix<Double> dBdH = ComputeTensorialMaterialParameter(HVec, ElemNum);
    double frobeniusNormSq;
    double determinant;
    double determinantSq; 
    double discriminant;
    double maxEigenvalueATA;
    if (dim_ == 2) {
      frobeniusNormSq = dBdH[0][0]*dBdH[0][0] + dBdH[0][1]*dBdH[0][1] + dBdH[1][0]*dBdH[1][0] + dBdH[1][1]*dBdH[1][1];
      determinant = dBdH[0][0]*dBdH[1][1]*dBdH[1][1] - dBdH[0][1]*dBdH[0][1]*dBdH[1][0]*dBdH[1][0];
      determinantSq = determinant * determinant;
      discriminant = std::sqrt(std::pow(frobeniusNormSq, 2) - 4 * determinantSq);
      maxEigenvalueATA = (frobeniusNormSq + discriminant) / 2.0;
      return std::sqrt(maxEigenvalueATA);
    } else {
      std::vector<std::vector<double>> ATA(3, std::vector<double>(3));
      for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
              ATA[i][j] = dBdH[0][i] * dBdH[0][j] + 
                          dBdH[1][i] * dBdH[1][j] + 
                          dBdH[2][i] * dBdH[2][j];
          }
      }
      // Power Iteration to find the largest eigenvalue of ATA
      std::vector<double> v = {1.0, 1.0, 1.0}; // Initial non-zero vector
      double maxEigenvalue = 0.0;
      int num_iterations = 20;
      for (int i = 0; i < num_iterations; ++i) {
          // Multiply: new_v = ATA * v
          std::vector<double> new_v(3, 0.0);
          for(int row = 0; row < 3; ++row) {
              for(int col = 0; col < 3; ++col) {
                  new_v[row] += ATA[row][col] * v[col];
              }
          }

          // Find the norm (magnitude) of the new vector. This value converges to the max eigenvalue.
          maxEigenvalue = std::sqrt(new_v[0] * new_v[0] + new_v[1] * new_v[1] + new_v[2] * new_v[2]);

          // Normalize the vector for the next iteration
          v[0] = new_v[0] / maxEigenvalue;
          v[1] = new_v[1] / maxEigenvalue;
          v[2] = new_v[2] / maxEigenvalue;
      }
      return std::sqrt(maxEigenvalue);
    }
  }

    

  Vector<Double> EBHysteresis::GetFluxDensity(Vector<Double> HVec, const Integer ElemNum) 
  {
    Vector<Double> B(dim_);
    UInt idx = ElemNum2Idx_[ElemNum];

    LOG_DBG3(eb) << "\n\t HVec = " << HVec.ToString();
    Vector<Double> M;
    M = Evaluate(HVec, idx);

    LOG_DBG3(eb) << "\n\t M = " << M.ToString();
    for (UInt i = 0; i < dim_; ++i)
    {
      B[i] = mu0_ * (HVec[i] + M[i]);
    }
    return B;
  }




  Vector<Double> EBHysteresis::GetFluxDensity(Vector<Double> HVec, const Integer ElemNum,
                                              LocPointMapped lpm, PtrCoefFct stressCoef)
  {
    Vector<Double> B(dim_);
    UInt idx = ElemNum2Idx_[ElemNum];
    Vector<Double> sigma;
    stressCoef->GetVector(sigma, lpm);

    // // dirty hack for 2d setups: add zero z-component to the stress tensor
    // if(dim_ == 2){
    //   sigma.Push_back(0.0);
    //   sigma.Push_back(0.0);
    //   sigma.Push_back(0.0);
    //   sigma[5] = sigma[2];
    //   sigma[2] = 0.0;
    // }
    SMSM_model_->Register_stress(sigma);

    LOG_DBG3(eb) << "\n\t sigma = " << sigma.ToString();
    LOG_DBG3(eb) << "\n\t HVec = " << HVec.ToString();
    Vector<Double> M;
    M = Evaluate(HVec, idx);

    LOG_DBG3(eb) << "\n\t M = " << M.ToString();

    for (UInt i = 0; i < dim_; ++i)
    {
      B[i] = mu0_ * (HVec[i] + M[i]);
    }
    return B;
  }




  Matrix<Double> EBHysteresis::ComputeTensorialMaterialParameter(Vector<Double> HVec, const Integer ElemNum)
  {
    UInt idx = ElemNum2Idx_[ElemNum];
    // the idea of this if is that the mu_ should only get updated for a new iteration (we might have several integration points
    // in one element and since we only use ONE hysteresis operator for one element, this is correct).
    if(iterTracker4Mu_ != mp_->GetExprVars(MathParser::GLOB_HANDLER, "iterationCounter")){
      iterTracker4Mu_ = mp_->GetExprVars(MathParser::GLOB_HANDLER, "iterationCounter");
      LOG_DBG3(eb) << "\n\t Trigger new iteration"<< std::endl;
      alreadyHasMu_.Init(false);
    }
   
    if(alreadyHasMu_[idx] == true){
      if(idx == 1){
        LOG_DBG2(eb) << "\t mu from alreadyHasMu_ = " << mu_[idx].ToString();
      }
      return mu_[idx];
    }
    Vector<Double> test_dmu_domega;
    test_dmu_domega = ComputeMaterialParamaterDerivative(HVec, idx);
    
    Vector<Double> M;
    Matrix<Double> mu(dim_, dim_);
    mu.InitValue(0.0);

    // To obtain a good starting point for the quasi-Newton method in the first time step
    // and first Newton iteration the Jacobian is approximated by forward finite differences
    //if(timeStep_== 1 && iterTracker4Mu_ == 1){
    if(iterTracker4Mu_ == 1)
    {
      Vector<Double> M;
      M = Evaluate(HVec, idx);
      StdVector<Double> B_k(dim_);
      for(UInt i = 0; i < dim_; i++)
      {
          B_k[i] = mu0_ * (HVec[i] + M[i]);
      }
      mu = EvaluateLocalMuFiniteDifferences(HVec, B_k, idx);

      mu_[idx] = mu;  
      alreadyHasMu_[idx] = true;
      return mu;
    }

    M = Evaluate(HVec, idx);

    // get the values for H and M from the last timestep and get deltaH and deltaB
    StdVector<Double> B_k(dim_);
    StdVector<Double> B_k_0(dim_);
    StdVector<Double> delta_H(dim_);
    StdVector<Double> delta_B(dim_);
    for (UInt i = 0; i < dim_; i++)
    {
      B_k[i] = mu0_ * (HVec[i] + M[i]);
      B_k_0[i] = mu0_ * (Htotal_prev_[idx][i] + Mprev_iter_[idx][i]);
      delta_H[i] = HVec[i] - Htotal_prev_[idx][i];
      delta_B[i] = B_k[i] - B_k_0[i];
    }

    if ((numS_ > 1) || (anhyst_type_ == "multiscale_anhysteresis") || (anhyst_formula_ == "pacejka"))
    { // hysteretic case
      switch (jacobian_method_)
      {
      case 1:
        // use finite differences
            if(idx == 1){
              LOG_DBG2(eb) << "\t use EvaluateLocalMuFiniteDifferences" ;
            }
        mu = EvaluateLocalMuFiniteDifferences(HVec, B_k, idx);
        break;
      case 2:
        // use Broyden method
            if(idx == 1){
              LOG_DBG2(eb) << "\t use EvaluateLocalMuGBM" ;
            }
        mu = EvaluateLocalMuGBM(delta_H, delta_B, idx);
        break;
      case 3:
        // use simple finite differences
           if(idx == 1){
              LOG_DBG2(eb) << "\t use EvaluateLocalMu" ;
            }
        mu = EvaluateLocalMu(delta_H, delta_B, idx);
        break;
      case 4:
        // use simple finite differences
            if(idx == 1){
              LOG_DBG2(eb) << "\t use EvaluateLocalMuBFGS" ;
            }
        mu = EvaluateLocalMuBFGS(delta_H, delta_B, idx);
        break;
      }
    }
    else
    { // nonlinear case (only anhysteresis)
      if(idx == 1){
        LOG_DBG2(eb) << "\t mu from real EvaluateLocalMuAnhystersisOnly";
      }
      mu = EvaluateLocalMuAnhystersisOnly(HVec, idx);
    }

    for (UInt i = 0; i < dim_; ++i) {
        Htotal_prev_[idx][i] = HVec[i];
        Mprev_iter_[idx][i] = M[i];
    }
    mu_[idx] = mu;    
    alreadyHasMu_[idx] = true;
    return mu;
  }


  void EBHysteresis::AllowUpdates(bool allow)
  {
    saveTmpStageVecs_ = allow;
  }

  void EBHysteresis::UpdateStates()
  {

    #pragma omp critical
        {
          HxS_n_ = HxS_n_tmp_;
          HyS_n_ = HyS_n_tmp_;
          HzS_n_ = HzS_n_tmp_;
          MxS_n_ = MxS_n_tmp_;
          MyS_n_ = MyS_n_tmp_;
          MzS_n_ = MzS_n_tmp_;
        }
  }

  /*
  =========================================================================================
      Methods regarding the inverse scheme
  =========================================================================================
  */
  Vector<Double> EBHysteresis::ComputeMaterialParamaterDerivative(Vector<Double> HVec, UInt idx){
    Vector<Double> dmu_domega(numS_);
    Vector<Double> M;
    Double numerator, denominator;
    Double H_norm, M_norm;
    if (dim_ == 2){ // 2-D case
      // DEFINE AND INIT. ALL NEEDED VARIABLES (2-D Case)
      StdVector<Double> Mx_j, My_j;
      Mx_j.Resize(numS_, 0.0);
      My_j.Resize(numS_, 0.0);
      Mx_j = MxS_n_tmp_[idx];
      My_j = MyS_n_tmp_[idx];
      Double Hx = HVec[0];
      Double Hy = HVec[1];
      Double H_norm = std::sqrt(std::pow(Hx,2) + std::pow(Hy,2));
      M = Evaluate(HVec, idx);
      Double Mx = M[0];
      Double My = M[1];
      Double M_norm = std::sqrt(std::pow(Mx,2) + std::pow(My,2));

      // ACTUAL COMPUTATION ACCORDING TO EQUATION (2-D Case)
      for (UInt jdx = 0; jdx < numS_; jdx++){
        numerator = (Hx*Mx_j[jdx]+Hy*My_j[jdx]) + (Mx*Mx_j[jdx]+My*My_j[jdx]);
        denominator = std::sqrt(std::pow(H_norm,2) + 2*(Hx*Mx+Hy*My) + std::pow(M_norm,2));
        dmu_domega[jdx] = (mu0_/H_norm)*(numerator/denominator);
      }
    } else { // 3-D case
      // DEFINE AND INIT. ALL NEEDED VARIABLES (3-D Case)
      StdVector<Double> Mx_j, My_j, Mz_j;
      Mx_j.Resize(numS_, 0.0);
      My_j.Resize(numS_, 0.0);
      Mz_j.Resize(numS_, 0.0);
      Mx_j = MxS_n_tmp_[idx];
      My_j = MyS_n_tmp_[idx];
      Mz_j = MzS_n_tmp_[idx];
      Double Hx = HVec[0];
      Double Hy = HVec[1];
      Double Hz = HVec[2];
      M = Evaluate(HVec, idx);
      Double Mx = M[0];
      Double My = M[1];
      Double Mz = M[2];
      // ACTUAL COMPUTATION ACCORDING TO EQUATION (3-D Case)
      H_norm = std::sqrt(std::pow(Hx,2) + std::pow(Hy,2) + std::pow(Hz,2));
      M_norm = std::sqrt(std::pow(Mx,2) + std::pow(My,2) + std::pow(Mz,2));
      for (UInt jdx = 0; jdx < numS_; jdx++){
        numerator = (Hx*Mx_j[jdx]+Hy*My_j[jdx]+Hz*Mz_j[jdx]) + (Mx*Mx_j[jdx]+My*My_j[jdx]+Mz*Mz_j[jdx]);
        denominator = std::sqrt(std::pow(H_norm,2) + 2*(Hx*Mx+Hy*My+Hz*Mz) + std::pow(M_norm,2));
        dmu_domega[jdx] = (mu0_/H_norm)*(numerator/denominator);
      }      
    }
    
    // RETURN
    return dmu_domega;
  }



  /*
  =========================================================================================
      Different versions for evaluating the local Jacobian (mu)
  =========================================================================================
  */
  Matrix<Double> EBHysteresis::EvaluateLocalMuAnhystersisOnly(Vector<Double> HVec, UInt idx)
  {
    Matrix<Double> mu; 
    mu.Resize(dim_,dim_);
    Matrix<Double> chi; 
    chi.Resize(dim_,dim_);
    Matrix<Double> T2; 
    T2.Resize(dim_,dim_);
    Matrix<Double> identity; 
    identity.Resize(dim_,dim_); 

    Double t0, t1, t3, factor1,factor2,factor3;
    Double mu0 = 1.256637061e-06;

    // 2D case (x,y)
    if (dim_ == 2){
      t0 = std::sqrt(HVec[0]*HVec[0] + HVec[1]*HVec[1]);
      t1 = t0/A_;
      t3 = std::atan(t1);
      factor1 = (2*Ps_/(t1*t1+1))/(M_PI*t0*t0*A_);
      factor2 = (2*Ps_*t3)/(M_PI*t0*t0*t0);
      factor3 = (2*Ps_*t3)/(M_PI*t0);

      for (UInt i = 0; i < dim_; ++i){
        for (UInt j = 0; j < dim_; ++j){
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

      for (UInt i = 0; i < dim_; ++i){
        for (UInt j = 0; j < dim_; ++j){
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

      eh_x[0] = h; 
      eh_x[1] = 0;
      eh_y[0] = 0; 
      eh_y[1] = h;
      for(UInt i = 0; i < dim_; i++){
        H_incx[i] = HVec[i] + eh_x[i];
        H_incy[i] = HVec[i] + eh_y[i];
      }
      // Evaluate EB Hysteresis Model for the increment AND without updating the stage values!!
      M_incx = Evaluate(H_incx, idx);
      M_incy = Evaluate(H_incy, idx);
      for (UInt i = 0; i < dim_; i++)
      {
        B_incx[i] = mu0_ * (H_incx[i] + M_incx[i]);
        B_incy[i] = mu0_ * (H_incy[i] + M_incy[i]);
      }
      // calculate the mu tesnor (dB/dH) for every entry
      mu[0][0] = (B_incx[0] - B_k[0]) / h; 
      mu[0][1] = 0;
      mu[1][0] = 0;
      mu[1][1] = (B_incy[1] - B_k[1]) / h;

    }else{ // 3D case (x,y,z)

      eh_x[0] = h; 
      eh_x[1] = 0; 
      eh_x[2] = 0;
      eh_y[0] = 0; 
      eh_y[1] = h; 
      eh_y[2] = 0;
      eh_z[0] = 0; 
      eh_z[1] = 0; 
      eh_z[2] = h;
      for(UInt i = 0; i < dim_; i++){
        H_incx[i] = HVec[i] + eh_x[i];
        H_incy[i] = HVec[i] + eh_y[i];
        H_incz[i] = HVec[i] + eh_z[i];
      }
      // Evaluate EB Hysteresis Model for the increment AND without updating the stage values!!
      M_incx = Evaluate(H_incx, idx);
      M_incy = Evaluate(H_incy, idx);
      M_incz = Evaluate(H_incz, idx);
      for (UInt i = 0; i < dim_; i++)
      {
        B_incx[i] = mu0_ * (H_incx[i] + M_incx[i]);
        B_incy[i] = mu0_ * (H_incy[i] + M_incy[i]);
        B_incz[i] = mu0_ * (H_incz[i] + M_incz[i]);
      }
      // calculate the mu tensor (dB/dH) for every entry
      mu[0][0] = (B_incx[0] - B_k[0]) / h; 
      mu[0][1] = 0;
      mu[0][2] = 0;
      mu[1][0] = 0;
      mu[1][1] = (B_incy[1] - B_k[1]) / h; 
      mu[1][2] = 0;
      mu[2][0] = 0;
      mu[2][1] = 0;
      mu[2][2] = (B_incz[2] - B_k[2]) / h;

    } 
    /* if (idx == 15){
      std::cout << "FD" << std::endl;
      std::cout << "mu[0][0]: " << mu[0][0] << ", mu[0][1]: " << mu[0][1] << ", mu[0][2]: " << mu[0][2] << std::endl;
      std::cout << "mu[1][0]: " << mu[1][0] << ", mu[1][1]: " << mu[1][1] << ", mu[1][2]: " << mu[1][2] << std::endl;
      std::cout << "mu[2][0]: " << mu[2][0] << ", mu[2][1]: " << mu[2][1] << ", mu[2][2]: " << mu[2][2] << std::endl;
    } */
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
              mu[i][i] = e[i][i];
          } 
      }
    if (idx == 32){
      printf("mu_simpleFD = %f\n",mu[0][0]);
    }
    return mu;
  }

Matrix<Double> EBHysteresis::EvaluateLocalMuBFGS(StdVector<Double> dH, StdVector<Double> dB, UInt idx){

      // define needed variables
      Matrix<Double> yyT(dim_, dim_);
      Matrix<Double> BxBxT(dim_, dim_);
      Matrix<Double> B_k1(dim_, dim_);
      Matrix<Double> B = mu_[idx];
      Double yTx;
      Double xTBx;
      StdVector<Double> y = dB;
      StdVector<Double> x = dH;
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
          B_k1[0][0] = mu0_;
          B_k1[1][1] = 0;
          B_k1[0][1] = 0;
          B_k1[1][0] = mu0_;
        }
        if ( (std::isinf(B_k1[0][0])) || (std::isinf(B_k1[1][1])) || (std::isinf(B_k1[0][1])) || (std::isinf(B_k1[1][0])) ) {
          B_k1[0][0] = mu0_;
          B_k1[1][1] = 0;
          B_k1[0][1] = 0;
          B_k1[1][0] = mu0_;
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


  Matrix<Double> EBHysteresis::EvaluateLocalMuDFP(StdVector<Double> dH, StdVector<Double> dB, UInt idx)
  {
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
              mu[i][i] = e[i][i]; 
          } 
      }
    }else{ // 3D version
      
      
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
      for (UInt i = 0; i < dim_; i++)
      {
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


  Matrix<Double> EBHysteresis::EvaluateLocalMuGBM(StdVector<Double> dH, StdVector<Double> dB, UInt idx)
  {
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
    }else{ // 3D version
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



  /*
  =========================================================================================
      Different versions of Evaluate function, depending on which models are used
  =========================================================================================
  */
  Vector<Double> EBHysteresis::Evaluate(Vector<Double> Hn, UInt idx)
  {
    Vector<Double> ret;
    StdVector<Double> weight(numS_);
    weight.Init(1.0 / numS_);
    StdVector<Double> chi(numS_);
    for (UInt i = 0; i < numS_; i++)
    {
      // chi[i] = chi_factor_ * i / (numS_ - 1) + 1e-22;
      chi[i] = kappa_file_[i];
      weight[i] = omega_file_[i];
    }

    if(anhyst_type_=="analytic_anhysteresis"){
      if(anhyst_formula_=="pacejka"){
        //check pacejka parameters
        if(m_sat_pacejka_<=0){
          EXCEPTION("m_sat cannot be negative or 0");
        }
        if(a_pacejka_<0){
          EXCEPTION("a cannot be negative");
        }
        if((b_pacejka_*c_pacejka_)<0){
          EXCEPTION("b times c cannot be negative");
        }     
      }
    }

    double sum = 0.0;
    for (size_t i = 0; i < weight.size(); ++i) {
      sum += weight[i];
    }
    //check weight distribution
    if((sum-1)>0.01 ||(sum-1)<-0.01){
      EXCEPTION("The sum of the weights must be equal to 1!");
    }


    if (dim_ == 2)
    {// make deciscion not on anhyst type but on analytic 
      if ((approx_type_ == "fullEB") && (anhyst_type_ == "multiscale_anhysteresis"))
      {
        // fullEB + Multiscale mode
        #pragma omp critical
        ret = Eval_2D_EBM_MSM(Hn, saveTmpStageVecs_, idx, weight, chi);
      }
      else if ((approx_type_ == "approxVPM") && (anhyst_type_ == "multiscale_anhysteresis"))
      {
        // VPM + Multiscale model
        #pragma omp critical
        ret = Eval_2D_VPM_MSM(Hn, saveTmpStageVecs_, idx, weight, chi);
      }
      else if ((approx_type_ == "approxVPM") && (anhyst_type_ == "analytic_anhysteresis"))
      {
        // VPM + atan anhysteresis
        ret = Eval_2D_VPM(Hn, saveTmpStageVecs_, idx, weight, chi);
      }
      else if ((approx_type_ == "fullEB") && (anhyst_type_ == "analytic_anhysteresis"))
      {
        // fullEB + atan anhysteresis
        ret = Eval_2D_EBM_ATAN(Hn, saveTmpStageVecs_, idx, weight, chi);
      }
    }
    else if (dim_ == 3)
    {
      if ((approx_type_ == "fullEB") && (anhyst_type_ == "multiscale_anhysteresis"))
      {
        // fullEB + Multiscale model
        #pragma omp critical
        ret = Eval_3D_EBM_MSM(Hn, saveTmpStageVecs_, idx, weight, chi);
      }
      else if ((approx_type_ == "approxVPM") && (anhyst_type_ == "multiscale_anhysteresis"))
      {
        // VPM + Multiscale model
        #pragma omp critical
        ret = Eval_3D_VPM_MSM(Hn, saveTmpStageVecs_, idx, weight, chi);
      }
      else if ((approx_type_ == "approxVPM") && (anhyst_type_ == "analytic_anhysteresis"))
      {
        // VPM + atan anhysteresis
        ret = Eval_3D_VPM(Hn, saveTmpStageVecs_, idx, weight, chi);
      }
      else if ((approx_type_ == "fullEB") && (anhyst_type_ == "analytic_anhysteresis"))
      {
        // fullEB + atan anhysteresis
        EXCEPTION("EBM+ATAN currently only implemented for 2D!");
      }
    }
    return ret;
  }


  Vector<Double> EBHysteresis::Eval_3D_VPM(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi)
  {
    Vector<Double> ret;
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

    Double HrxS_prev, HryS_prev, HrzS_prev, HrS, Px, Py, Pz, condition1, Hirr_x, Hirr_y, Hirr_z;

    StdVector<Double> &HxS_prev = HxS_n_[idx];
    StdVector<Double> &HyS_prev = HyS_n_[idx];
    StdVector<Double> &HzS_prev = HzS_n_[idx];

    for (UInt k = 0; k < numS_; k++)
    {
      HrxS_prev = HxS_prev[k];
      HryS_prev = HyS_prev[k];
      HrzS_prev = HzS_prev[k];
      condition1 = (std::pow((Hex_x - HrxS_prev) / chi[k], 2) + std::pow((Hex_y - HryS_prev) / chi[k], 2) + std::pow((Hex_z - HrzS_prev) / chi[k], 2));
      if (condition1 <= 1.0)
      {
        HrxS_sol[k] = HrxS_prev;
        HryS_sol[k] = HryS_prev;
        HrzS_sol[k] = HrzS_prev;
      }
      else
      {
        // dissipation now reached (just arrived at the border of the sphere)
        if (k == 0)
        {
          HrxS_sol[k] = Hex_x;
          HryS_sol[k] = Hex_y;
          HrzS_sol[k] = Hex_z;
        }
        else
        {
          //
          Hirr_x = chi[k] * (Hex_x - HrxS_prev) / std::sqrt((std::pow((Hex_x - HrxS_prev), 2) + std::pow((Hex_y - HryS_prev), 2)) + std::pow((Hex_z - HrzS_prev), 2));
          Hirr_y = chi[k] * (Hex_y - HryS_prev) / std::sqrt((std::pow((Hex_x - HrxS_prev), 2) + std::pow((Hex_y - HryS_prev), 2)) + std::pow((Hex_z - HrzS_prev), 2));
          Hirr_z = chi[k] * (Hex_z - HrzS_prev) / std::sqrt((std::pow((Hex_x - HrxS_prev), 2) + std::pow((Hex_y - HryS_prev), 2)) + std::pow((Hex_z - HrzS_prev), 2));
          HrxS_sol[k] = Hex_x - Hirr_x;
          HryS_sol[k] = Hex_y - Hirr_y;
          HrzS_sol[k] = Hex_z - Hirr_z;
        }
      }
      HrS = std::sqrt(std::pow(HrxS_sol[k], 2) + std::pow(HryS_sol[k], 2) + std::pow(HrzS_sol[k], 2));
      if (std::sqrt(std::pow(HrS, 2)) > 1.0e-12)
      { 
        if(anhyst_formula_=="atan"){
          MxS_sol[k] = (2.0 * Ps_ / M_PI) * std::atan(HrS / A_) * HrxS_sol[k] / HrS;
          MyS_sol[k] = (2.0 * Ps_ / M_PI) * std::atan(HrS / A_) * HryS_sol[k] / HrS;
          MzS_sol[k] = (2.0 * Ps_ / M_PI) * std::atan(HrS / A_) * HrzS_sol[k] / HrS;
        }
        else if(anhyst_formula_=="pacejka"){
          MxS_sol[k] = m_sat_pacejka_ * std::sin(std::atan(a_pacejka_*HrS+b_pacejka_*std::atan(c_pacejka_*HrS))) * HrxS_sol[k] / HrS;
          MyS_sol[k] = m_sat_pacejka_ * std::sin(std::atan(a_pacejka_*HrS+b_pacejka_*std::atan(c_pacejka_*HrS))) * HryS_sol[k] / HrS;
          MzS_sol[k] = m_sat_pacejka_ * std::sin(std::atan(a_pacejka_*HrS+b_pacejka_*std::atan(c_pacejka_*HrS))) * HrzS_sol[k] / HrS;
        }
      }
      else
      {
        MxS_sol[k] = 0.0;
        MyS_sol[k] = 0.0;
        MzS_sol[k] = 0.0;
      }
    }
    Px = 0.0;
    Py = 0.0;
    Pz = 0.0;
    for (UInt k = 0; k < numS_; k++)
    {
      Px += weight[k] * MxS_sol[k];
      Py += weight[k] * MyS_sol[k];
      Pz += weight[k] * MzS_sol[k];
    }
    if (saveTmpStageVecs)
    {
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
    //std::cout << "Name: " << pinning_forces_weight_ << std::endl;
    return ret;
  }


  Vector<Double> EBHysteresis::Eval_2D_VPM(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi)
  {
    Vector<Double> ret;
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

    Double HrxS_prev, HryS_prev, HrS, Px, Py, condition1, Hirr_x, Hirr_y;

    StdVector<Double> &HxS_prev = HxS_n_[idx];
    StdVector<Double> &HyS_prev = HyS_n_[idx];

    for (UInt k = 0; k < numS_; k++)
    {
      HrxS_prev = HxS_prev[k];
      HryS_prev = HyS_prev[k];
      condition1 = (std::pow((Hex_x - HrxS_prev) / chi[k], 2) + std::pow((Hex_y - HryS_prev) / chi[k], 2));
      if (condition1 <= 1.0)
      {
        HrxS_sol[k] = HrxS_prev;
        HryS_sol[k] = HryS_prev;
      }
      else
      {
        // dissipation now reached (just arrived at the border of the sphere)
        if (k == 0)
        {
          HrxS_sol[k] = Hex_x;
          HryS_sol[k] = Hex_y;
        }
        else
        {
          //
          Hirr_x = chi[k] * (Hex_x - HrxS_prev) / std::sqrt((std::pow((Hex_x - HrxS_prev), 2) + std::pow((Hex_y - HryS_prev), 2)));
          Hirr_y = chi[k] * (Hex_y - HryS_prev) / std::sqrt((std::pow((Hex_x - HrxS_prev), 2) + std::pow((Hex_y - HryS_prev), 2)));
          HrxS_sol[k] = Hex_x - Hirr_x;
          HryS_sol[k] = Hex_y - Hirr_y;
        }
      }
      HrS = std::sqrt(std::pow(HrxS_sol[k], 2) + std::pow(HryS_sol[k], 2));
      if (std::sqrt(std::pow(HrS, 2)) > 1.0e-12)
      { 
        if(anhyst_formula_=="atan"){
          MxS_sol[k] = (2.0 * Ps_ / M_PI) * std::atan(HrS / A_) * HrxS_sol[k] / HrS;
          MyS_sol[k] = (2.0 * Ps_ / M_PI) * std::atan(HrS / A_) * HryS_sol[k] / HrS;
        }
        if(anhyst_formula_=="pacejka"){
          MxS_sol[k] = m_sat_pacejka_ * std::sin(std::atan(a_pacejka_*HrS+b_pacejka_*std::atan(c_pacejka_*HrS))) * HrxS_sol[k] / HrS;
          MyS_sol[k] = m_sat_pacejka_ * std::sin(std::atan(a_pacejka_*HrS+b_pacejka_*std::atan(c_pacejka_*HrS))) * HryS_sol[k] / HrS;
        }
      }
      else
      {
        MxS_sol[k] = 0.0;
        MyS_sol[k] = 0.0;
      }
    }
    Px = 0.0;
    Py = 0.0;
    for (UInt k = 0; k < numS_; k++)
    {
      Px += weight[k] * MxS_sol[k];
      Py += weight[k] * MyS_sol[k];
    }
    if (saveTmpStageVecs)
    {
      HxS_n_tmp_[idx] = HrxS_sol;
      HyS_n_tmp_[idx] = HryS_sol;
      MxS_n_tmp_[idx] = MxS_sol;
      MyS_n_tmp_[idx] = MyS_sol;
    }
    ret.Push_back(Px);
    ret.Push_back(Py);
    return ret;
  }


  Vector<Double> EBHysteresis::Eval_2D_VPM_MSM(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi)
  {
    Vector<Double> Hn3d(3);
    Hn3d[0] = Hn[0];
    Hn3d[1] = Hn[1];
    Hn3d[2] = 0.0;
    Vector<Double> M3d = Eval_3D_VPM_MSM(Hn3d, saveTmpStageVecs, idx, weight, chi);
    Vector<Double> M(2);
    M[0] = M3d[0];
    M[1] = M3d[1];
    return M;
  }

  Vector<Double> EBHysteresis::Eval_3D_VPM_MSM(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi)
  {
    Vector<Double> ret;
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

    Double HrxS_prev, HryS_prev, HrzS_prev, HrS, Px, Py, Pz, condition1, Hirr_x, Hirr_y, Hirr_z;

    StdVector<Double> &HxS_prev = HxS_n_[idx];
    StdVector<Double> &HyS_prev = HyS_n_[idx];
    StdVector<Double> &HzS_prev = HzS_n_[idx];

    for (UInt k = 0; k < numS_; k++)
    {
      HrxS_prev = HxS_prev[k];
      HryS_prev = HyS_prev[k];
      HrzS_prev = HzS_prev[k];
      condition1 = (std::pow((Hex_x - HrxS_prev) / chi[k], 2) + std::pow((Hex_y - HryS_prev) / chi[k], 2) + std::pow((Hex_z - HrzS_prev) / chi[k], 2));
      if (condition1 <= 1.0)
      {
        HrxS_sol[k] = HrxS_prev;
        HryS_sol[k] = HryS_prev;
        HrzS_sol[k] = HrzS_prev;
      }
      else
      {
        // dissipation now reached (just arrived at the border of the sphere)
        if (k == 0)
        {
          HrxS_sol[k] = Hex_x;
          HryS_sol[k] = Hex_y;
          HrzS_sol[k] = Hex_z;
        }
        else
        {
          //
          Hirr_x = chi[k] * (Hex_x - HrxS_prev) / std::sqrt((std::pow((Hex_x - HrxS_prev), 2) + std::pow((Hex_y - HryS_prev), 2)) + std::pow((Hex_z - HrzS_prev), 2));
          Hirr_y = chi[k] * (Hex_y - HryS_prev) / std::sqrt((std::pow((Hex_x - HrxS_prev), 2) + std::pow((Hex_y - HryS_prev), 2)) + std::pow((Hex_z - HrzS_prev), 2));
          Hirr_z = chi[k] * (Hex_z - HrzS_prev) / std::sqrt((std::pow((Hex_x - HrxS_prev), 2) + std::pow((Hex_y - HryS_prev), 2)) + std::pow((Hex_z - HrzS_prev), 2));
          HrxS_sol[k] = Hex_x - Hirr_x;
          HryS_sol[k] = Hex_y - Hirr_y;
          HrzS_sol[k] = Hex_z - Hirr_z;
        }
      }
      HrS = std::sqrt(std::pow(HrxS_sol[k], 2) + std::pow(HryS_sol[k], 2) + std::pow(HrzS_sol[k], 2));
      if (std::sqrt(std::pow(HrS, 2)) > 1.0e-12)
      {
        // Use the MSMS for calculation of the new stage magnetization vector
        StdVector<Double> dirH(3);
        dirH[0] = HrxS_sol[k] / HrS;
        dirH[1] = HryS_sol[k] / HrS;
        dirH[2] = HrzS_sol[k] / HrS;
        LOG_DBG3(eb) << "\n\t INPUT (H) OF SMSM: [" << HrxS_sol[k] << "," << HryS_sol[k] << ", " << HrzS_sol[k] << "]";
        SMSM_model_->Eval(HrS, dirH);
        Vector<Double> M = SMSM_model_->GetM();
        LOG_DBG3(eb) << "\n\t OUTPUT (M) OF SMSM: " << M.ToString();

        // MSM version
        MxS_sol[k] = M[0];
        MyS_sol[k] = M[1];
        MzS_sol[k] = M[2];
      }
      else
      {
        MxS_sol[k] = 0.0;
        MyS_sol[k] = 0.0;
        MzS_sol[k] = 0.0;
      }
    }
    Px = 0.0;
    Py = 0.0;
    Pz = 0.0;
    for (UInt k = 0; k < numS_; k++)
    {
      Px += weight[k] * MxS_sol[k];
      Py += weight[k] * MyS_sol[k];
      Pz += weight[k] * MzS_sol[k];
    }
    if (saveTmpStageVecs)
    {
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
    return ret;
  }

  

  Vector<Double> EBHysteresis::Eval_2D_EBM_ATAN(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi)
  {
    Vector<Double> ret;
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
    StdVector<Double> &HxS_prev = HxS_n_[idx];
    StdVector<Double> &HyS_prev = HyS_n_[idx];
    StdVector<Double> &MxS_prev = MxS_n_[idx];
    StdVector<Double> &MyS_prev = MyS_n_[idx];

    for (UInt k = 0; k < numS_; k++)
    {
      HrxS_prev = HxS_prev[k];
      HryS_prev = HyS_prev[k];
      MxSprev = MxS_prev[k];
      MySprev = MyS_prev[k];
      // if(std::sqrt( std::pow(Hex_x - HrxS_prev, 2) + std::pow(Hex_y - HryS_prev, 2) ) <= chi[k]){
      if ((std::pow((Hex_x - HrxS_prev) / chi[k], 2) + std::pow((Hex_y - HryS_prev) / chi[k], 2)) <= 1.0)
      {
        HrxS_sol[k] = HrxS_prev;
        HryS_sol[k] = HryS_prev;
      }
      else
      {
        StdVector<Double> i_phi_initial(2);
        i_phi_initial.Init(0.0);
        // dissipation now reached (just arrived at the border of the sphere)
        if (k == 0)
        {
          i_phi_initial[0] = 0.0;
          i_phi_initial[1] = 0.0;
        }
        else
        {
          // use the direction of the vector play model as initial direction
          i_phi_initial[0] = (1.0 / chi[k] * (HrxS_prev - Hex_x)) / (std::pow((Hex_x - HrxS_prev) / chi[k], 2) + std::pow((Hex_y - HryS_prev) / chi[k], 2));
          i_phi_initial[1] = (1.0 / chi[k] * (HryS_prev - Hex_y)) / (std::pow((Hex_x - HrxS_prev) / chi[k], 2) + std::pow((Hex_y - HryS_prev) / chi[k], 2));
        }
        phi = std::atan2(i_phi_initial[1], i_phi_initial[0]); // already in rad
        err = 1.0;
        iter = 0;
        while (err > 1.0e-3 && iter < 10)
        {
          ux = Hex_x + chi[k] * std::cos(phi);
          uy = Hex_y + chi[k] * std::sin(phi);
          uabs = std::sqrt(std::pow(ux, 2) + std::pow(uy, 2));
          ux1 = -chi[k] * std::sin(phi);
          uy1 = chi[k] * std::cos(phi);
          ux2 = -chi[k] * std::cos(phi);
          uy2 = -chi[k] * std::sin(phi);
          Man = 2 * Ps_ / M_PI * std::atan(uabs / A_);
          Man1 = 2 * Ps_ / (A_ * M_PI) * (1.0 / (1.0 + std::pow(uabs / A_, 2)));
          g1 = (Man / uabs) * (ux * ux1 + uy * uy1) - MxSprev * ux1 - MySprev * uy1;
          g2 = (uabs * Man1 - Man) / (std::pow(uabs, 3)) * std::pow(ux * ux1 + uy * uy1, 2) + Man / uabs * (std::pow(ux1, 2) + std::pow(uy1, 2)) + Man / uabs * (ux * ux2 + uy * uy2) - MxSprev * ux2 - MySprev * uy2;
          if (std::sqrt(std::pow(g2, 2)) < 1e-8)
          {
            phiNew = 0;
          }
          else
          {
            phiNew = phi - g1 / g2;
          }
          err = std::sqrt(std::pow(phiNew - phi, 2));
          phi = phiNew;
          iter++;
        }
        error[k] = err;
        dir[k] = phi;
        numIter[k] = iter;
        HrxS_sol[k] = Hex_x + chi[k] * std::cos(phi);
        HryS_sol[k] = Hex_y + chi[k] * std::sin(phi);
      }
      HrS = std::sqrt(std::pow(HrxS_sol[k], 2) + std::pow(HryS_sol[k], 2));
      if (std::sqrt(std::pow(HrS, 2)) > 1.0e-12)
      {
        MxS_sol[k] = (2.0 * Ps_ / M_PI) * std::atan(HrS / A_) * HrxS_sol[k] / HrS;
        MyS_sol[k] = (2.0 * Ps_ / M_PI) * std::atan(HrS / A_) * HryS_sol[k] / HrS;
      }
      else
      {
        MxS_sol[k] = 0.0;
        MyS_sol[k] = 0.0;
      }
    }
    Px = 0.0;
    Py = 0.0;
    for (UInt k = 0; k < numS_; k++)
    {
      Px += weight[k] * MxS_sol[k];
      Py += weight[k] * MyS_sol[k];
    }
    if (saveTmpStageVecs)
    {
      HxS_n_tmp_[idx] = HrxS_sol;
      HyS_n_tmp_[idx] = HryS_sol;
      MxS_n_tmp_[idx] = MxS_sol;
      MyS_n_tmp_[idx] = MyS_sol;
    }
    ret.Push_back(Px);
    ret.Push_back(Py);
    return ret;
  }



  Vector<Double> EBHysteresis::Eval_3D_EBM_MSM(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi)
  {
    Vector<Double> Hn2d(2);
    Hn2d[0] = Hn[0];
    Hn2d[1] = Hn[1];
    Vector<Double> M2d = Eval_2D_EBM_MSM(Hn2d, saveTmpStageVecs, idx, weight, chi);
    Vector<Double> M(3);
    M[0] = M2d[0];
    M[1] = M2d[1];
    M[2] = 0.0;
    return M;
  }

  Vector<Double> EBHysteresis::Eval_2D_EBM_MSM(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi)
  {
    Vector<Double> ret;
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

    Double HrxS_prev, HryS_prev, MxSprev, MySprev, phi, err, F_prime, F_prime_prime,
           ux, uy, uabs, ux1, uy1, ux2, uy2, Man, Man1, g1, g2, phiNew, HrS, Px, Py;
    Matrix<Double> dM_dHrev(2,2);
    
    
    UInt iter;
    StdVector<Double> &HxS_prev = HxS_n_[idx];
    StdVector<Double> &HyS_prev = HyS_n_[idx];
    StdVector<Double> &MxS_prev = MxS_n_[idx];
    StdVector<Double> &MyS_prev = MyS_n_[idx];
    LOG_DBG3(eb) << "\n\t##################### STARTING EVAL OF EB2D MODEL #####################";
    LOG_DBG3(eb) << "\n\t Stress tensor = "<<SMSM_model_->GetSigma().ToString();    
    LOG_DBG3(eb) << "\n\t chi = "<< chi.ToString();

    for (UInt k = 0; k < numS_; k++)
    {
      LOG_DBG3(eb) << "\n\t Starting with k = "<<k;
      HrxS_prev = HxS_prev[k];
      HryS_prev = HyS_prev[k];
      MxSprev = MxS_prev[k];
      MySprev = MyS_prev[k];
      // if(std::sqrt( std::pow(Hex_x - HrxS_prev, 2) + std::pow(Hex_y - HryS_prev, 2) ) <= chi[k]){
      if ((k > 0) && ((std::pow((Hex_x - HrxS_prev) / chi[k], 2) + std::pow((Hex_y - HryS_prev) / chi[k], 2)) <= 1.0))
      {
        LOG_DBG3(eb) << "\n\t\t Inside of sphere, keep Hrev";
        HrxS_sol[k] = HrxS_prev;
        HryS_sol[k] = HryS_prev;
      }
      else if ((k > 0) && ((std::pow((Hex_x - HrxS_prev) / chi[k], 2) + std::pow((Hex_y - HryS_prev) / chi[k], 2)) > 1.0))
      {
        LOG_DBG3(eb) << "\n\t\t Touching sphere -> perform optimization";
        StdVector<Double> i_phi_initial(2);
        i_phi_initial.Init(0.0);
        // dissipation now reached (just arrived at the border of the sphere)
        if (k == 0)
        {
          i_phi_initial[0] = 0.0;
          i_phi_initial[1] = 0.0;
        }
        else
        {
          // use the direction of the vector play model as initial direction
          i_phi_initial[0] = -(1.0 / chi[k] * (HrxS_prev - Hex_x)) / (std::pow((Hex_x - HrxS_prev) / chi[k], 2) + std::pow((Hex_y - HryS_prev) / chi[k], 2));
          i_phi_initial[1] = -(1.0 / chi[k] * (HryS_prev - Hex_y)) / (std::pow((Hex_x - HrxS_prev) / chi[k], 2) + std::pow((Hex_y - HryS_prev) / chi[k], 2));
          LOG_DBG3(eb) << "\n\t\t Initial VPM solution: i_phi = "<< i_phi_initial.ToString();
        }
        phi = std::atan2(i_phi_initial[1], i_phi_initial[0]); // already in rad
        LOG_DBG3(eb) << "\n\t\t Initial VPM solution: phi = "<< phi;

        err = 1.0;
        iter = 0;
        F_prime = 1.0;

        LOG_DBG3(eb) << "\n\t\t Start optimization";
        while ( (err > 1.0e-5) && (std::abs(F_prime) > 1.0e-5) && (iter < 10))
        {
          Calc_derivs(F_prime, F_prime_prime, dM_dHrev, MxSprev, MySprev,
                      Hex_x, Hex_y, HrxS_prev, HryS_prev, phi, chi[k]);
          LOG_DBG3(eb) << "\n\t\t\t dM_dHrev = "<< dM_dHrev.ToString();
          LOG_DBG3(eb) << "\n\t\t\t F_prime = "<< F_prime;
          LOG_DBG3(eb) << "\n\t\t\t F_prime_prime = "<< F_prime_prime;
          if (std::abs(F_prime_prime) < 1e-6)
          {
            phiNew = phi;
          }
          else
          {
            Double eta = Energy_linesearch(Hex_x, Hex_y,
                                       HrxS_prev, HryS_prev,
                                       MxSprev, MySprev, phi, chi[k],
                                       F_prime, F_prime_prime);
            LOG_DBG3(eb) << "\n\t\t\t Linesearch eta = "<< eta;
            phiNew = phi - (F_prime / F_prime_prime)*eta;
          }
          err = std::sqrt(std::pow(phiNew - phi, 2));
          phi = phiNew;
          iter++;
          LOG_DBG3(eb) << "\n\t\t\t New phi = "<< phi;
          LOG_DBG3(eb) << "\n\t\t\t error = "<< err;
        }


        error[k] = err;
        dir[k] = phi;
        numIter[k] = iter;
        HrxS_sol[k] = Hex_x - chi[k] * std::cos(phi);
        HryS_sol[k] = Hex_y - chi[k] * std::sin(phi);
      }
      else
      {
        HrxS_sol[k] = Hex_x;
        HryS_sol[k] = Hex_y;
      }

      HrS = std::sqrt(std::pow(HrxS_sol[k], 2) + std::pow(HryS_sol[k], 2));
      if (std::sqrt(std::pow(HrS, 2)) > 1.0e-12)
      {


        // Use the MSMS for calculation of the new stage magnetization vector
        StdVector<Double> dirH(2);
        dirH[0] = HrxS_sol[k] / HrS;
        dirH[1] = HryS_sol[k] / HrS;
        LOG_DBG3(eb) << "\n\t INPUT (H) OF SMSM: [" << HrxS_sol[k] << "," << HryS_sol[k] << ", " << 0.0 << "]";
        SMSM_model_->Eval(HrS, dirH);
        Vector<Double> M = SMSM_model_->GetM();
        LOG_DBG3(eb) << "\n\t OUTPUT (M) OF SMSM: " << M.ToString();

        // MSM version
        MxS_sol[k] = M[0];
        MyS_sol[k] = M[1];
      }
      else
      {
        MxS_sol[k] = 0.0;
        MyS_sol[k] = 0.0;
      }
    }
    Px = 0.0;
    Py = 0.0;
    for (UInt k = 0; k < numS_; k++)
    {
      Px += weight[k] * MxS_sol[k];
      Py += weight[k] * MyS_sol[k];
    }
    if (saveTmpStageVecs)
    {
      HxS_n_tmp_[idx] = HrxS_sol;
      HyS_n_tmp_[idx] = HryS_sol;
      MxS_n_tmp_[idx] = MxS_sol;
      MyS_n_tmp_[idx] = MyS_sol;
    }
    ret.Push_back(Px);
    ret.Push_back(Py);
    return ret;
  }

  void EBHysteresis::Calc_derivs(Double &F_prime, Double &F_prime_prime, Matrix<Double> &dM_dHrev,
                                Double MxSprev, Double MySprev,
                                Double Hex_x, Double Hex_y,
                                Double HrxS_prev, Double HryS_prev,
                                Double phi, Double chi)
  {
    Double ux = Hex_x - chi * std::cos(phi);
    Double uy = Hex_y - chi * std::sin(phi);
    Double uabs = std::sqrt(ux*ux + uy*uy);

    StdVector<Double> dirH(2);
    dirH[0] = ux / uabs;
    dirH[1] = uy / uabs;

    StdVector<Double> de_dphi(2);
    StdVector<Double> d2e_dphi2(2);
    de_dphi[0] = -std::sin(phi);   de_dphi[1] = std::cos(phi);
    d2e_dphi2[0] = -std::cos(phi); d2e_dphi2[1] = -std::sin(phi);
    SMSM_model_->Eval(uabs, dirH);
    dM_dHrev = SMSM_model_->GetdMdH();
    Vector<Double> M = SMSM_model_->GetM();
    
    Double M_Mprev_x = M[0] - MxSprev;
    Double M_Mprev_y = M[1] - MySprev;

    LOG_DBG3(eb) << "\n\t phi: " << phi;
    LOG_DBG3(eb) << "\n\t M[0] - MxSprev: " << M_Mprev_x;
    LOG_DBG3(eb) << "\n\t M[1] - MySprev: " << M_Mprev_y;
    LOG_DBG3(eb) << "\n\t de_dphi: " << de_dphi.ToString();
    LOG_DBG3(eb) << "\n\t d2e_dphi2: " << d2e_dphi2.ToString();
    F_prime =  chi*(M_Mprev_x*de_dphi[0] + M_Mprev_y*de_dphi[1]);

    Double dMdHrev_dedphi_x = dM_dHrev[0][0]*de_dphi[0] + dM_dHrev[0][1]*de_dphi[1];
    Double dMdHrev_dedphi_y = dM_dHrev[1][0]*de_dphi[0] + dM_dHrev[1][1]*de_dphi[1];
    Double dedephi_dMdHrev_dedphi = dMdHrev_dedphi_x*de_dphi[0] + dMdHrev_dedphi_y*de_dphi[1];
    Double d2edphi2_M = d2e_dphi2[0] * M_Mprev_x + d2e_dphi2[1] * M_Mprev_y;
    F_prime_prime = -chi*chi * dedephi_dMdHrev_dedphi + chi*d2edphi2_M;
  }


  Double EBHysteresis::Energy_linesearch(Double Hx, Double Hy,
                                       Double Hprev_x, Double Hprev_y,
                                       Double Mprev_x, Double Mprev_y, Double phi_orig, Double chi,
                                       Double& F_prime_orig, Double& F_prime_prime_orig)
  {

    Double G_prev = 1.0e12;
    Double s = 1.0e-5;
    Double e = 1.0;
    UInt num = 3;
    Double etaopt = 0.0;
    Double F_prime, F_prime_prime, phi_k, G, eta;
    Matrix<Double> dM_dHrev;


    for (UInt eta_iter = 0 ; eta_iter < num+1; eta_iter++)
    {
      eta = s+eta_iter*(e-s)/num;
      phi_k = phi_orig - eta * F_prime_orig/F_prime_prime_orig;
      
      
      LOG_DBG3(eb) << "\n\t phi_orig: " << phi_orig;
      LOG_DBG3(eb) << "\n\t F_prime_orig: " << F_prime_orig;
      LOG_DBG3(eb) << "\n\t F_prime_prime_orig: " << F_prime_prime_orig;
      LOG_DBG3(eb) << "\n\t eta: " << eta;
      LOG_DBG3(eb) << "\n\t phi_k: " << phi_k;


      Calc_derivs(F_prime, F_prime_prime, dM_dHrev, Mprev_x, Mprev_y,
                      Hx, Hy, Hprev_x, Hprev_y, phi_k, chi);
      LOG_DBG3(eb) << "\n\t F_prime: " << F_prime;
                      
      G = F_prime * F_prime_orig/F_prime_prime_orig;
      LOG_DBG3(eb) << "\n\t G: " << G;

      if(std::abs(G) <= G_prev)
      {
         G_prev = std::abs(G);
         etaopt = eta;
      }
    }
    return etaopt;
  }
}

