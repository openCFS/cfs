#include <torch/script.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <array>
#include <sstream>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "SMSM.hh"


namespace CoupledField
{

DEFINE_LOG(smsm, "SMSM")

  SMSM::SMSM(Double Ms, Double AS, Double K1, Double K2, Double lambda100, Double lambda111, UInt dim)
  {
    l100_ = lambda100;
    l111_ = lambda111;
    AS_ = AS;
    K1_ = K1;
    K2_ = K2;
    Ms_ = Ms;
    mu0_ = 4*M_PI*1e-7;
    dim_ = dim;

    if(dim_ == 3){
      // Define the file path
      // std::string filePath = "/home/klaus/Devel/CFS_SRC/cfs/source/Materials/Models/TABSPHEREI4S4_2562.txt";
      std::string filePath = "/Users/kroppert/Devel/CFS_SRC/cfs/source/Materials/Models/TABSPHEREI4S4_2562.txt";
      std::ifstream file(filePath);

      // Check if the file was successfully opened
      if (!file.is_open()) {
          EXCEPTION("Error: Unable to open file " << filePath);
      }

      numRows_ = 2562;
      const UInt numCols = 3;

      // Pre-allocate vector to store the data
      TABgamma_.Resize(numRows_, Vector<Double>(numCols));

      //data.reserve(numRows);

      // Read the file line by line
      std::string line;
      size_t currentRow = 0;
      while (std::getline(file, line) && currentRow < numRows_) {
          std::stringstream ss(line);
          for (unsigned int col = 0; col < numCols; ++col) {
              if (!(ss >> TABgamma_[currentRow][col])) {
                  EXCEPTION("Error: Failed to parse line " << currentRow + 1 << ": " << line);
                  break;
              }
          }
          currentRow++;
      }
      // Close the file
      file.close();

      SIGMAloc_.Resize(3, 3);
      SIGMAloc_.Init();

      epsmu11_.Resize(numRows_);
      epsmu12_.Resize(numRows_);
      epsmu13_.Resize(numRows_);
      epsmu22_.Resize(numRows_);
      epsmu23_.Resize(numRows_);
      epsmu33_.Resize(numRows_);
      Wan_.Resize(numRows_);

      for (UInt i = 0; i < numRows_; ++i) {
          epsmu11_[i] = 1.5 * l100_ * (std::pow(TABgamma_[i][0], 2) - 1.0 / 3.0);
          epsmu22_[i] = 1.5 * l100_ * (std::pow(TABgamma_[i][1], 2) - 1.0 / 3.0);
          epsmu33_[i] = 1.5 * l100_ * (std::pow(TABgamma_[i][2], 2) - 1.0 / 3.0);

          epsmu12_[i] = 1.5 * l111_ * (TABgamma_[i][0]*TABgamma_[i][1]);
          epsmu13_[i] = 1.5 * l111_ * (TABgamma_[i][0]*TABgamma_[i][2]);
          epsmu23_[i] = 1.5 * l111_ * (TABgamma_[i][1]*TABgamma_[i][2]);

          Wan_[i] = K2_ * (std::pow(TABgamma_[i][0], 2) *
                          std::pow(TABgamma_[i][1], 2) *
                          std::pow(TABgamma_[i][2], 2)) +
                    K1_ * (std::pow(TABgamma_[i][0], 2) * std::pow(TABgamma_[i][1], 2) +
                          std::pow(TABgamma_[i][1], 2) * std::pow(TABgamma_[i][2], 2) +
                          std::pow(TABgamma_[i][2], 2) * std::pow(TABgamma_[i][0], 2));
      }
    }else{
      // 2D CASE
      // here we can create a uniform meshin in 2d (circle)

      numRows_ = 360*3;
      const UInt numCols = 2;
      Double delta_alpha = 2*M_PI/numRows_;

      // Pre-allocate vector to store the data
      TABgamma_.Resize(numRows_, Vector<Double>(numCols));
      for(UInt i = 0; i < numRows_; ++i){
        TABgamma_[i][0] = cos(i*delta_alpha);
        TABgamma_[i][1] = sin(i*delta_alpha);
      }


      // 2D case, here we assume a plane stress state!!!!!
      SIGMAloc_.Resize(2, 2);
      SIGMAloc_.Init();

      epsmu11_.Resize(numRows_);
      epsmu12_.Resize(numRows_);
//      epsmu13_.Resize(numRows_, 0.0);
      epsmu22_.Resize(numRows_);
//      epsmu23_.Resize(numRows_, 0.0);
      epsmu33_.Resize(numRows_, -1.5 * l100_ * 1.0 / 3.0);
      Wan_.Resize(numRows_);

      for (UInt i = 0; i < numRows_; ++i) {
          epsmu11_[i] = 1.5 * l100_ * (std::pow(TABgamma_[i][0], 2) - 1.0 / 3.0);
          epsmu22_[i] = 1.5 * l100_ * (std::pow(TABgamma_[i][1], 2) - 1.0 / 3.0);
          epsmu12_[i] = 1.5 * l111_ * (TABgamma_[i][0]*TABgamma_[i][1]);
          Wan_[i] = K2_ * (std::pow(TABgamma_[i][0], 2) *
                          std::pow(TABgamma_[i][1], 2)) +
                    K1_ * (std::pow(TABgamma_[i][0], 2) * std::pow(TABgamma_[i][1], 2));
      }
    }
    LOG_DBG3(smsm) << "\n\t max(Wan_) after precompute step = " << Wan_.Max() 
                   << "\n\t min(Wan_) after precompute step = " << Wan_.Min();
  }

  SMSM::~SMSM()
  {
  }

  void SMSM::Eval(Double valH, StdVector<Double> dirHloc)
  {
    if(dim_ == 2){
      this->Eval2D(valH, dirHloc);
    }else{
      this->Eval3D(valH, dirHloc);
    }
  }


  void SMSM::Eval3D(Double valH, StdVector<Double> dirHloc)
  {

  Vector<Double> Wtot(numRows_);
  Double As_times_Wtot_exp = 0.0; 
  Double gamma0, gamma1, gamma2, e11, e12, e13, e23, e33, e22, f;
  for (UInt i = 0; i < numRows_; ++i) {
    gamma0 = TABgamma_[i][0];
    gamma1 = TABgamma_[i][1];
    gamma2 = TABgamma_[i][2];
    e11 = epsmu11_[i];
    e12 = epsmu12_[i];
    e13 = epsmu13_[i];
    e23 = epsmu23_[i];
    e33 = epsmu33_[i];
    e22 = epsmu22_[i];

    Wtot[i] =   // mechanic part
                -(e11 * SIGMAloc_[0][0]
                + e22 * SIGMAloc_[1][1]
                + e33 * SIGMAloc_[2][2]
                + 2.0 * e12 * SIGMAloc_[0][1]
                + 2.0 * e13 * SIGMAloc_[0][2]
                + 2.0 * e23 * SIGMAloc_[1][2])
                + // magnetic part
                -mu0_*Ms_*valH *  (
                        gamma0 *dirHloc[0]
                       +gamma1 *dirHloc[1]
                       +gamma2 *dirHloc[2])
                + // anisotropic part
                Wan_[i];
  }
  // the problem is that the argument of the exp function is too large (also occurs sometimes in machine learning,
  // where a similar distribution called softmax is employed). The trick is to add a constant, in our case a negative value,
  // to the argument in both the nominator and denominator, which cancells out but prevents overflows
  // Get the maximum value of Wtot and then compute the denominator:
  Double Wtot_max = Wtot.Max();
  #pragma omp parallel for reduction(+:As_times_Wtot_exp)
  for (UInt i = 0; i < numRows_; ++i) {
    As_times_Wtot_exp += std::exp(-AS_*Wtot[i] - Wtot_max*AS_);
  } 

  LOG_DBG3(smsm) << "\n\t SIGMAloc_ " << SIGMAloc_.ToString()
                 << "\n\t max(Wtot) after precompute step = " << Wtot.Max() 
                 << "\n\t min(Wtot) after precompute step = " << Wtot.Min();
  LOG_DBG3(smsm) << "\n\t Sum(exp(As*Wtot) = " << As_times_Wtot_exp;
  if(std::isinf(As_times_Wtot_exp)){
    EXCEPTION("Sum(exp(As*Wtot) = "<<As_times_Wtot_exp);
  }



  Double eps_11 = 0.0;
  Double eps_12 = 0.0;
  Double eps_13 = 0.0;
  Double eps_22 = 0.0;
  Double eps_23 = 0.0;
  Double eps_33 = 0.0;
  Double dMdH_00 = 0.0;
  Double dMdH_11 = 0.0;
  Double dMdH_22 = 0.0;
  Double dMdH_01 = 0.0;
  Double dMdH_02 = 0.0;
  Double dMdH_10 = 0.0;
  Double dMdH_20 = 0.0;
  Double dMdH_21 = 0.0;
  Double dMdH_12 = 0.0;
  Double MMoy0 = 0.0;
  Double MMoy1 = 0.0;
  Double MMoy2 = 0.0;
  Double falpha_Malpha_0, falpha_Malpha_1, falpha_Malpha_2;
  Vector<Double> Malpha(3);
  Vector<Double> falpha_Malpha(3);
  for (UInt i = 0; i < numRows_; ++i) {
    gamma0 = TABgamma_[i][0];
    gamma1 = TABgamma_[i][1];
    gamma2 = TABgamma_[i][2];
    e11 = epsmu11_[i];
    e12 = epsmu12_[i];
    e13 = epsmu13_[i];
    e23 = epsmu23_[i];
    e33 = epsmu33_[i];
    e22 = epsmu22_[i];


    // here we must also use the scaling-trick described above for the nominator
    f = std::exp(-AS_*Wtot[i] - Wtot_max*AS_) / As_times_Wtot_exp;
    
    MMoy0 += Ms_*f*gamma0;
    MMoy1 += Ms_*f*gamma1;
    MMoy2 += Ms_*f*gamma2;

    eps_11 += Ms_*f * e11;
    eps_12 += Ms_*f * e12;
    eps_13 += Ms_*f * e13;
    eps_22 += Ms_*f * e22;
    eps_23 += Ms_*f * e23;
    eps_33 += Ms_*f * e33;

    falpha_Malpha_0 = gamma0*Ms_ * f;
    falpha_Malpha_1 = gamma1*Ms_ * f;
    falpha_Malpha_2 = gamma2*Ms_ * f;
      
    dMdH_00 += falpha_Malpha_0*gamma0*Ms_;
    dMdH_11 += falpha_Malpha_1*gamma1*Ms_;
    dMdH_22 += falpha_Malpha_2*gamma2*Ms_;
    dMdH_01 += falpha_Malpha_0*gamma1*Ms_;
    dMdH_02 += falpha_Malpha_0*gamma2*Ms_;
    dMdH_10 += falpha_Malpha_1*gamma0*Ms_;
    dMdH_20 += falpha_Malpha_2*gamma0*Ms_;
    dMdH_21 += falpha_Malpha_2*gamma1*Ms_;
    dMdH_12 += falpha_Malpha_1*gamma2*Ms_;
  }

 


  dMdH_00 = mu0_*AS_*(dMdH_00 - MMoy0*MMoy0);
  dMdH_11 = mu0_*AS_*(dMdH_11 - MMoy1*MMoy1);
  dMdH_22 = mu0_*AS_*(dMdH_22 - MMoy2*MMoy2);
  dMdH_01 = mu0_*AS_*(dMdH_01 - MMoy0*MMoy1);
  dMdH_02 = mu0_*AS_*(dMdH_02 - MMoy0*MMoy2);
  dMdH_10 = mu0_*AS_*(dMdH_10 - MMoy1*MMoy0);
  dMdH_20 = mu0_*AS_*(dMdH_20 - MMoy2*MMoy0);
  dMdH_21 = mu0_*AS_*(dMdH_21 - MMoy2*MMoy1);
  dMdH_12 = mu0_*AS_*(dMdH_12 - MMoy1*MMoy2);


  dMdH_.Resize(3,3);
  dMdH_.Init();
  dMdH_[0][0] = dMdH_00;
  dMdH_[1][1] = dMdH_11;
  dMdH_[2][2] = dMdH_22;
  dMdH_[0][1] = dMdH_01;
  dMdH_[0][2] = dMdH_02;
  dMdH_[1][0] = dMdH_10;
  dMdH_[2][0] = dMdH_20;
  dMdH_[2][1] = dMdH_21;
  dMdH_[1][2] = dMdH_12;

  MMoy_.Resize(3);
  MMoy_.Init(0.0);
  MMoy_[0] = MMoy0;
  MMoy_[1] = MMoy1;
  MMoy_[2] = MMoy2;

  LOG_DBG3(smsm) << "\n\t dMdH_ = " << dMdH_.ToString();

 
  epsmumoy_.Resize(3,3);
  epsmumoy_[0][0] = eps_11;
  epsmumoy_[0][1] = eps_12;
  epsmumoy_[0][2] = eps_13;
  epsmumoy_[1][0] = eps_12;
  epsmumoy_[1][1] = eps_22;
  epsmumoy_[1][2] = eps_23;
  epsmumoy_[2][0] = eps_13;
  epsmumoy_[2][1] = eps_23;
  epsmumoy_[2][2] = eps_33;

  }


  void SMSM::Eval2D(Double valH, StdVector<Double> dirHloc)
  {

  Vector<Double> Wtot(numRows_);
  Double As_times_Wtot_exp = 0.0; 
  Double gamma0, gamma1, gamma2, e11, e12, e13, e23, e33, e22, f;
  for (UInt i = 0; i < numRows_; ++i) {
    gamma0 = TABgamma_[i][0];
    gamma1 = TABgamma_[i][1];
    e11 = epsmu11_[i];
    e12 = epsmu12_[i];
    e22 = epsmu22_[i];
    e33 = epsmu33_[i];

    Wtot[i] =   // mechanic part
                -(e11 * SIGMAloc_[0][0]
                + e22 * SIGMAloc_[1][1]
                + 2.0 * e12 * SIGMAloc_[0][1])
                + // magnetic part
                -mu0_*Ms_*valH *  (
                        gamma0 *dirHloc[0]
                       +gamma1 *dirHloc[1])
                + // anisotropic part
                Wan_[i];
  }
  // the problem is that the argument of the exp function is too large (also occurs sometimes in machine learning,
  // where a similar distribution called softmax is employed). The trick is to add a constant, in our case a negative value,
  // to the argument in both the nominator and denominator, which cancells out but prevents overflows
  // Get the maximum value of Wtot and then compute the denominator:
  Double Wtot_max = Wtot.Max();
  #pragma omp parallel for reduction(+:As_times_Wtot_exp)
  for (UInt i = 0; i < numRows_; ++i) {
    As_times_Wtot_exp += std::exp(-AS_*Wtot[i] - Wtot_max*AS_);
  } 

  LOG_DBG3(smsm) << "\n\t SIGMAloc_ " << SIGMAloc_.ToString()
                 << "\n\t max(Wtot) after precompute step = " << Wtot.Max() 
                 << "\n\t min(Wtot) after precompute step = " << Wtot.Min();
  LOG_DBG3(smsm) << "\n\t Sum(exp(As*Wtot) = " << As_times_Wtot_exp;
  if(std::isinf(As_times_Wtot_exp)){
    EXCEPTION("Sum(exp(As*Wtot) = "<<As_times_Wtot_exp);
  }



  Double eps_11 = 0.0;
  Double eps_12 = 0.0;
  Double eps_22 = 0.0;
  Double eps_33 = 0.0;
  Double dMdH_00 = 0.0;
  Double dMdH_11 = 0.0;
  Double dMdH_01 = 0.0;
  Double dMdH_10 = 0.0;
  Double MMoy0 = 0.0;
  Double MMoy1 = 0.0;
  Double falpha_Malpha_0, falpha_Malpha_1;
  Vector<Double> Malpha(3);
  Vector<Double> falpha_Malpha(3);
  for (UInt i = 0; i < numRows_; ++i) {
    gamma0 = TABgamma_[i][0];
    gamma1 = TABgamma_[i][1];
    e11 = epsmu11_[i];
    e12 = epsmu12_[i];
    e33 = epsmu33_[i];
    e22 = epsmu22_[i];

    // here we must also use the scaling-trick described above for the nominator
    f = std::exp(-AS_*Wtot[i] - Wtot_max*AS_) / As_times_Wtot_exp;
    
    MMoy0 += Ms_*f*gamma0;
    MMoy1 += Ms_*f*gamma1;

    eps_11 += Ms_*f * e11;
    eps_12 += Ms_*f * e12;
    eps_22 += Ms_*f * e22;
    eps_33 += Ms_*f * e33;

    falpha_Malpha_0 = gamma0*Ms_ * f;
    falpha_Malpha_1 = gamma1*Ms_ * f;
      
    dMdH_00 += falpha_Malpha_0*gamma0*Ms_;
    dMdH_11 += falpha_Malpha_1*gamma1*Ms_;
    dMdH_01 += falpha_Malpha_0*gamma1*Ms_;
    dMdH_10 += falpha_Malpha_1*gamma0*Ms_;
  }

  dMdH_00 = mu0_*AS_*(dMdH_00 - MMoy0*MMoy0);
  dMdH_11 = mu0_*AS_*(dMdH_11 - MMoy1*MMoy1);
  dMdH_01 = mu0_*AS_*(dMdH_01 - MMoy0*MMoy1);
  dMdH_10 = mu0_*AS_*(dMdH_10 - MMoy1*MMoy0);

  dMdH_.Resize(2,2);
  dMdH_.Init();
  dMdH_[0][0] = dMdH_00;
  dMdH_[1][1] = dMdH_11;
  dMdH_[0][1] = dMdH_01;
  dMdH_[1][0] = dMdH_10;
  
  MMoy_.Resize(2);
  MMoy_.Init(0.0);
  MMoy_[0] = MMoy0;
  MMoy_[1] = MMoy1;
  
  LOG_DBG3(smsm) << "\n\t dMdH_ = " << dMdH_.ToString();

 
  epsmumoy_.Resize(3,3);
  epsmumoy_[0][0] = eps_11;
  epsmumoy_[0][1] = eps_12;
  epsmumoy_[0][2] = 0.0;
  epsmumoy_[1][0] = eps_12;
  epsmumoy_[1][1] = eps_22;
  epsmumoy_[1][2] = 0.0;
  epsmumoy_[2][0] = 0.0;
  epsmumoy_[2][1] = 0.0;
  epsmumoy_[2][2] = 0.0;
  }


  void SMSM::Register_stress(Vector<Double> sigma)
  {
    if(dim_ == 2){
      SIGMAloc_[0][0] = sigma[0]; SIGMAloc_[0][1] = sigma[2]; 
      SIGMAloc_[1][0] = sigma[2]; SIGMAloc_[1][1] = sigma[1]; 
    }else{
      SIGMAloc_[0][0] = sigma[0]; SIGMAloc_[0][1] = sigma[5]; SIGMAloc_[0][2] = sigma[4];
      SIGMAloc_[1][0] = sigma[5]; SIGMAloc_[1][1] = sigma[1]; SIGMAloc_[1][2] = sigma[3]; 
      SIGMAloc_[2][0] = sigma[4]; SIGMAloc_[2][1] = sigma[3]; SIGMAloc_[2][2] = sigma[2];
    }

    //SIGMAloc_.Mult(1.0e6); // in Pascals
    //Double TABvalsig = 0.0;
    //SIGMAloc_.Mult(TABvalsig);

  }












  // we only need to access these variable inside of this cc file, so internal linkage is ok
  static torch::jit::script::Module module_;
  static torch::Device device_ = torch::Device(torch::kCPU);
  // In class as member variables
  static torch::Tensor dirHloc_x_;
  static torch::Tensor dirHloc_y_;
  static torch::Tensor input_;


  SMSM_PYTORCH::SMSM_PYTORCH(Double Ms, Double AS, Double K1, Double K2, Double lambda100, Double lambda111, UInt dim)
  {
    l100_ = lambda100;
    l111_ = lambda111;
    AS_ = AS;
    K1_ = K1;
    K2_ = K2;
    Ms_ = Ms;
    mu0_ = 4*M_PI*1e-7;
    dim_ = dim;
    


    // std::string model_path = "/Users/kroppert/Devel/CFS_SRC/cfs/source/Materials/Models/nn_tracedmodel_nostress.pt";
    std::string model_path = "/Users/kroppert/Devel/CFS_SRC/cfs/source/Materials/Models/nn_tracedmodel_nostress.pt";
    try {
        // Deserialize the ScriptModule from a file using torch::jit::load().
        module_ = torch::jit::load(model_path);
        module_ = torch::jit::optimize_for_inference(module_);
    } catch (const c10::Error& e) {
        std::cerr << "Error loading the model\n";
    }

    // Determine if MPS (Metal Performance Shaders) or CPU is available ... take care, moving stuff around is expensive!!!
    // if (torch::hasMPS()) {
    //     device_ = torch::Device(torch::kMPS);
    // }

    dirHloc_x_ = torch::zeros({1});
    dirHloc_y_ = torch::zeros({1});
    // ADAPT NUMBER OF INPUT PARAMETERS FOR THE SPECIFIC MODEL!!!!!!
    input_ = torch::zeros({1, 5}, torch::TensorOptions().dtype(torch::kFloat32).device(device_));
    // Move the model to the correct device
    module_.to(device_);
  }

  SMSM_PYTORCH::~SMSM_PYTORCH()
  {
  }

  void SMSM_PYTORCH::Eval(Double valH, StdVector<Double> dirHloc)
  {
    if(dim_ == 2){
      this->Eval2D(valH, dirHloc);
    }else{
      this->Eval3D(valH, dirHloc);
    }
  }


  void SMSM_PYTORCH::Eval3D(Double valH, StdVector<Double> dirHloc)
  {


  }


  void SMSM_PYTORCH::Eval2D(Double valH, StdVector<Double> dirHloc)
  {
    
  // Disable gradient calculations
  torch::NoGradGuard no_grad;

  dirHloc_x_.fill_(dirHloc[0]);
  dirHloc_y_.fill_(dirHloc[1]);

  // Calculate phi_H tensors
  torch::Tensor phi_H = torch::atan2(dirHloc_y_, dirHloc_x_);

  float sigma_value = 0.0; // hardcoded for now!!!!

  /* for the 
  nn_tracedmodel_nostress
  model 
  */
  // Fill the input tensor
  input_.index_put_({torch::indexing::Slice(), 0}, float(this->AS_));  // AS
  input_.index_put_({torch::indexing::Slice(), 1}, float(this->Ms_));  // Ms
  // input_.index_put_({torch::indexing::Slice(), 2}, sigma_value);  // sigma_value
  input_.index_put_({torch::indexing::Slice(), 2}, torch::cos(phi_H));  // cos(phi_H)
  input_.index_put_({torch::indexing::Slice(), 3}, torch::sin(phi_H));  // sin(phi_H)
  input_.index_put_({torch::indexing::Slice(), 4}, float(valH));  // Hm


  /* for the 
  gpytorch_tracedmodel_nostress
  model
  */ 
  // Fill the input tensor
  // input_.index_put_({torch::indexing::Slice(), 0}, torch::cos(phi_H));  // cos(phi_H)
  // input_.index_put_({torch::indexing::Slice(), 1}, torch::sin(phi_H));  // sin(phi_H)
  // input_.index_put_({torch::indexing::Slice(), 2}, float(valH));  // Hm




  // StandardScaler mean and scale values (extracted from Python)
  // for the cos(phi) and sin(phi) enter 0.0!
  std::vector<float> means_input = {2.96177837e-03, 1.40515493e+06, 0.0, 0.0, 1.57281174e+01};  
  std::vector<float> scales_input = {1.12018447e-03, 1.18703615e+05, 0.0, 0.0, 8.69155511e+03}; 
  std::vector<float> means_output = {389.51755523, -22.13245129};  
  std::vector<float> scales_output = {936219.22701107, 936385.22703299}; 
  // Standardize the input and output except for the columns in exception_array
  std::vector<int> exception_array = {2, 3};  // Don't standardize cos(phi_H) and sin(phi_H)


  // // Print the input tensor values before standardization
  // auto input_data_before = input_.accessor<float, 2>();  // Assuming it's a 2D tensor
  // std::cout << "Input tensor values (before standardization):" << std::endl;
  // for (int i = 0; i < input_.size(0); ++i) {
  //     for (int j = 0; j < input_.size(1); ++j) {
  //         std::cout << "input_[" << i << "][" << j << "] = " << input_data_before[i][j] << std::endl;
  //     }
  // }



  for (int col = 0; col < input_.size(1); ++col) {
      if (std::find(exception_array.begin(), exception_array.end(), col) != exception_array.end() || scales_input[col] == 0.0) {
          continue;  // Skip this column if it's in the exception array or if the scale is zero
      }
      // Standardize: (x - mean) / scale
      input_.index({torch::indexing::Slice(), col}) = 
          (input_.index({torch::indexing::Slice(), col}) - means_input[col]) / scales_input[col];
  }



  // // Print the input tensor values after standardization
  // auto input_data_after = input_.accessor<float, 2>();  // Assuming it's a 2D tensor
  // std::cout << "Input tensor values (after standardization):" << std::endl;
  // for (int i = 0; i < input_.size(0); ++i) {
  //     for (int j = 0; j < input_.size(1); ++j) {
  //         std::cout << "input_[" << i << "][" << j << "] = " << input_data_after[i][j] << std::endl;
  //     }
  // }


  // Prepare input as a vector of IValues (inputs for the model)
  std::vector<torch::jit::IValue> inputs;
  inputs.push_back(input_);

  // Execute the model and turn its output into a tensor (for the NN model)
  at::Tensor mean_output = module_.forward(inputs).toTensor();

  // Execute the model and turn its output into a tensor (for the GPR model the output is a tuple)
  // auto output_tuple = module_.forward(inputs).toTuple();  // Get the tuple output
  // at::Tensor mean_output = output_tuple->elements()[0].toTensor();  // Extract the mean
  // at::Tensor variance_output = output_tuple->elements()[1].toTensor();  // Extract the variance

  // Move the output back to CPU to print
  mean_output = mean_output.to(torch::kCPU);


  //Inverse transform the output
  for (int col = 0; col < mean_output.size(1); ++col) {
      if (scales_output[col] != 0.0) {
          // Inverse of the standardization: x_original = x_scaled * scale + mean
          mean_output.index({torch::indexing::Slice(), col}) = 
              mean_output.index({torch::indexing::Slice(), col}) * scales_output[col] + means_output[col];
      }
  }

  // Print ypred_test[:,0] and ypred_test[:,1] (first two columns of output)
  auto mean_output_data = mean_output.accessor<float, 2>();  // Assuming output is a 2D tensor

  // std::cout << "\ndirHloc:" << dirHloc.ToString() << std::endl;
  // std::cout << "\nvalH:" << valH << std::endl;
  // for (int i = 0; i < 2; ++i) {
  //     std::cout << "Mx: " << mean_output_data[i][0] << ", My: " << mean_output_data[i][1] << std::endl;
  // }


  MMoy_.Resize(2);
  MMoy_.Init(0.0);
  MMoy_[0] = mean_output_data[0][0];
  MMoy_[1] = mean_output_data[0][1];

  }


  void SMSM_PYTORCH::Register_stress(Vector<Double> sigma)
  {
    if(dim_ == 2){
      SIGMAloc_[0][0] = sigma[0]; SIGMAloc_[0][1] = sigma[2]; 
      SIGMAloc_[1][0] = sigma[2]; SIGMAloc_[1][1] = sigma[1]; 
    }else{
      SIGMAloc_[0][0] = sigma[0]; SIGMAloc_[0][1] = sigma[5]; SIGMAloc_[0][2] = sigma[4];
      SIGMAloc_[1][0] = sigma[5]; SIGMAloc_[1][1] = sigma[1]; SIGMAloc_[1][2] = sigma[3]; 
      SIGMAloc_[2][0] = sigma[4]; SIGMAloc_[2][1] = sigma[3]; SIGMAloc_[2][2] = sigma[2];
    }

    //SIGMAloc_.Mult(1.0e6); // in Pascals
    //Double TABvalsig = 0.0;
    //SIGMAloc_.Mult(TABvalsig);

  }
} // Namespace end
