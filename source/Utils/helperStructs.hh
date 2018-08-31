#ifndef HELPER_STRUCTS_2018
#define HELPER_STRUCTS_2018

#include "General/Environment.hh"

namespace CoupledField {
  
  // structs used for Hysteresis and coefFunctionHyst
  
  struct InitialInput {
    bool useInitialInput;
    bool prescribeOutput;
    bool scaleBySaturation;
    Vector<Double> inputVector;
  };
  
  struct ParameterInversion {
    // common parameter for inversion methods
    UInt inversionMethod;
    UInt maxNumIts;
    Double tolH;
    Double tolB;
    // trust region regularization > for LM only
    UInt maxNumRegIts;
    Double alphaRegStart;
    Double alphaRegMin;
    Double alphaRegMax;
    Double trustLow;
    Double trustMid;
    Double trustHigh;
    // linesearch parameter > for Newton and KrylovNewton
    UInt maxNumLSIts;    
    Double alphaLSMin;
    Double alphaLSMax;
    
    Double jacRes;
    int jacImplementation;
    bool stopLineSearchAtLocalMin;
    
    // parameter for projected Levenberg-Marquardt
    // > "Levenberg-Marquardt methods for constrained nonlinear equations with strong local convergence properties"
    //     -Kanzow,Yamashita,Fukushima
    Double projLM_mu;
    Double projLM_rho;
    Double projLM_beta;
    Double projLM_sigma; 
    Double projLM_gamma;  
    Double projLM_tau; 
    Double projLM_c; 
    Double projLM_p;
  };
   
  struct ParameterPreisachWeights {
    UInt numRows_;
    std::string weightType_;
    // weightType = 0 > use constant weight
    Double constWeight_;
    // weightType = 1 > muDat
    Double muDat_A_;
    Double muDat_h1_;
    Double muDat_sigma1_;
    Double muDat_eta_;
    // weightType = 2 > extended muDat (use muDat parameter in addition)
    Double muDat_h2_;
    Double muDat_sigma2_;
        
    // weightType = 3 > tensor given
    Matrix<Double> weightTensor_;
    
    // anhysteretic parameter
    Double anhysteretic_a_;
    Double anhysteretic_b_;
    Double anhysteretic_c_;
    Double anhystAtSat_normalized_;
    int anhysteretic_cInAtan_;
    bool anhystOnly_;
    bool anhystCountingToOutputSat_;
    // 0 (default): integrate over all weights, scale down to 1, use outputSat for scaling of Preisach operator
    // 1 (default for muDat): integrate over all weights, scale down to 1, use preisachSat for scaling of Preisach operator
    //                        with preisachSat = outputSat*(1.0 - anhystAtSat_normalized_)
    Double intOverWeights_;
    
    
  };
  
//   ParameterPreisachWeights POL_weightParams_;
//  ParameterPreisachWeights STRAIN_weightParams_;
//    
//  ParameterPreisachOperators POL_operatorParams_;
//  ParameterPreisachOperators STRAIN_operatorParams_;
//  
//  ParameterInversion LM_inversion_;
//  InitialInput POL_initial_;
  
  struct ParameterIrrStrainsAndCoupling {
    
    bool couplingDefined_inMatFile_;
    int strainForm_;
    bool useStrainForm_;
    
    int ci_size_;
    Matrix<Double> ci_;
    Double c1_;
    Double c2_;
    Double c3_;
    Double d0_;
    Double d1_;
    int irrStrainForm_;
    Double sSat_;
    int scaleTosSat_;
    int paramsForHalfRange_;
    
    bool ownHystOperator_;
    
  };
  
  
  struct ParameterPreisachOperators {
    
    int methodType_; // SCALAR or VECTOR (> CoefFunction::CoefDimType)
    std::string methodName_; // scalarPreisach, vectorPreisach_Sutor, vectorPreisach_Mayergoyz
    /*
     * common for all implemented operators
     * (scalar, vector sutor, vector mayergoyz isotropic)
     */
    Double inputSat_;
    Double outputSat_; // saturation of Polarization or Strains but not of flux!
    Double preisachSat_; // depends on anhysteretic parameter! 
    // if anhystCountingToOutputSat_ == false > preisachSat = outputSat = value of preisach operator alone at inputSat
    // if anhystCountingToOutputSat_ == true > preisachSat = outputSat*(1-anhystAtSat_normalized_)
    bool fieldsAlignedAboveSat_;
    bool hystOutputRestrictedToSat_; // meant is value of preisachSat_ here
    bool hasInverseModel_;
    
    /*
     * scalar model
     */
    Vector<Double> fixDirection_;
    
    /*
     * vector model sutor
     */
    UInt evalVersion_;
    Double rotResistance_;
    Double angularDistance_;
    Double angularResolution_;
    Double amplitudeResolution_;
    Double angularClipping_;
    UInt printOut_;
    UInt bmpResolution_;
    bool isClassical_;
    bool scaleUpToSaturation_;
    
    /*
     * vector model mayergoyz (isotropic!)
     */
    bool isIsotropic_;
    int numDirections_;
    int outputClipping_;
  };
  
} // end of CoupledField

#endif
