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

  SMSM::SMSM()
  {
    // Define the file path
    std::string filePath = "/home/kroppert/Devel/CFS_SRC/latest_trunkGit/cfs/source/Materials/Models/TABSPHEREI4S4_2562.txt";
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




    l100_ = 0.0;
    l111_ = 0.0;
    AS_ = 1.0e-3;
    K1_ = 0.0;
    K2_ = 0.0;
    Ms_ = 1.23e6;
    mu0_ = 1.256637061e-06;
    
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
    LOG_DBG3(smsm) << "\n\t max(Wan_) after precompute step = " << Wan_.Max() 
                   << "\n\t min(Wan_) after precompute step = " << Wan_.Min();
  }

  SMSM::~SMSM()
  {
  }

  void SMSM::Eval(Double valH, StdVector<Double> dirHloc)
  {

  Vector<Double> Wtot;
  Wtot.Resize(numRows_);
  Double As_times_Wtot_exp = 0.0; 
  for (UInt i = 0; i < numRows_; ++i) {
    Wtot[i] =   // mechanic part
                -(epsmu11_[i] * SIGMAloc_[0][0]
                + epsmu22_[i] * SIGMAloc_[1][1]
                + epsmu33_[i] * SIGMAloc_[2][2]
                + 2.0 * epsmu12_[i] * SIGMAloc_[0][1]
                + 2.0 * epsmu13_[i] * SIGMAloc_[0][2]
                + 2.0 * epsmu23_[i] * SIGMAloc_[1][2])
                + // magnetic part
                -mu0_*Ms_*valH *  (
                        TABgamma_[i][0] *dirHloc[0]
                       +TABgamma_[i][1] *dirHloc[1]
                       +TABgamma_[i][2] *dirHloc[2])
                + // anisotropic part
                Wan_[i];
    //As_times_Wtot_exp += std::exp(-AS_*Wtot[i]);
  }
  // the problem is that the argument of the exp function is too large (also occurs sometimes in machine learning,
  // where a similar distribution called softmax is employed). The trick is to add a constant, in our case a negative value,
  // to the argument in both the nominator and denominator, which cancells out but prevents overflows
  // Get the maximum value of Wtot and then compute the denominator:
  Double Wtot_max = Wtot.Max();
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

  Vector<Double> falpha;
  falpha.Resize(numRows_);
  
  MMoy_.Resize(3);
  MMoy_.Init(0.0);
  Double eps_11 = 0.0;
  Double eps_12 = 0.0;
  Double eps_13 = 0.0;
  Double eps_22 = 0.0;
  Double eps_23 = 0.0;
  Double eps_33 = 0.0;
  for (UInt i = 0; i < numRows_; ++i) {
    // here we must also use the scaling-trick described above for the nominator
    falpha[i] = std::exp(-AS_*Wtot[i] - Wtot_max*AS_) / As_times_Wtot_exp;
    MMoy_[0] += Ms_*falpha[i]*TABgamma_[i][0];
    MMoy_[1] += Ms_*falpha[i]*TABgamma_[i][1];
    MMoy_[2] += Ms_*falpha[i]*TABgamma_[i][2];

    eps_11 += Ms_*falpha[i] * epsmu11_[i];
    eps_12 += Ms_*falpha[i] * epsmu12_[i];
    eps_13 += Ms_*falpha[i] * epsmu13_[i];
    eps_22 += Ms_*falpha[i] * epsmu22_[i];
    eps_23 += Ms_*falpha[i] * epsmu23_[i];
    eps_33 += Ms_*falpha[i] * epsmu33_[i];
  }
  LOG_DBG3(smsm) << "\n\t max(falpha) after precompute step = " << falpha.Max() 
                 << "\n\t min(falpha) after precompute step = " << falpha.Min();

  // Synthesis of results
  
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

  void SMSM::Register_stress(Vector<Double> sigma)
  {
    SIGMAloc_[0][0] = sigma[0]; SIGMAloc_[0][1] = sigma[5]; SIGMAloc_[0][2] = sigma[4];
    SIGMAloc_[1][0] = sigma[5]; SIGMAloc_[1][1] = sigma[1]; SIGMAloc_[1][2] = sigma[3]; 
    SIGMAloc_[2][0] = sigma[4]; SIGMAloc_[2][1] = sigma[3]; SIGMAloc_[2][2] = sigma[2];
    SIGMAloc_.Mult(1.0e6); // in Pascals
    //Double TABvalsig = 0.0;
    //SIGMAloc_.Mult(TABvalsig);

  }

} // Namespace end
