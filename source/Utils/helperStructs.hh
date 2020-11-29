#ifndef HELPER_STRUCTS_2018
#define HELPER_STRUCTS_2018

#include "General/Environment.hh"

namespace CoupledField {

  // structs used for Hysteresis and coefFunctionHyst

  //! Data from tracing hyst operator
  struct TracedData {
    // for identification
//    std::string nameTag_;
//    UInt numberTag_;
//
//    bool operator<(const TracedData& a) const
//    {
//        return numberTag_ < a.numberTag_;
//    }

//    UInt tracedSteps_;

    Double maxSlope_;
    Double minSlope_;
    Double maxPolarization_;
    Double negCoercivity_;
    
//    bool forStrain_;
  };

  //! Initial input for hysteresis operator
  struct InitialInput {
    bool useInitialInput;
    bool prescribeOutput;
    bool scaleBySaturation;
    Vector<Double> inputVector;
  };
  
  //! Parameter for local inversion of hysteresis operator (i.e. inversion on element level)
  struct ParameterInversion {
    // common parameter for inversion methods
    localInversionFlag inversionMethod;
    UInt maxNumIts;
    Double tolH;
    Double tolB;
    bool tolH_useAsRelativeNorm;
    bool tolB_useAsRelativeNorm;
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

    // iteration via FP > value obtained from tracing hyst operator > see coefFunctionHyst
    Double minSlopeHystOperator;
    Double maxSlopeHystOperator;
    Double safetyFactor_C;

    bool printWarnings;
  };

  struct LinesearchParameter {
    // checkingOrder_: as name suggests, order in which multiple criteria shall be checked
    UInt checkingOrder_;

    bool operator<(const LinesearchParameter& a) const
    {
        return checkingOrder_ < a.checkingOrder_;
    }

    // nameTag_: for human readable output later on
    std::string nameTag_;

//    // lsType_: to distinguish differnet linesearches > defined in Environment.hh
////    typedef enum { INVALID_LINESEARCH = -2, NO_LINESEARCH = -1, LS_TRIAL_AND_ERROR = 0, LS_GOLDENSECTION = 1, LS_QUADINTERPOL = 2, LS_VECBASEDPOLY = 3,
////    LS_THREEPOINTPOLY = 4, LS_BACKTRACKING_ARMIJO = 5, LS_BACKTRACKING_ARMIJO_NONMONOTONIC = 6, LS_BACKTRACKING_GRADFREE = 7 } linesearchFlag;
//    linesearchFlag lsType_;

    /*
     * Common parameter for most linesearch types
     */
    // maxIter_: maximal number of iterations
    UInt maxIter_;

    // minEta_: minimal (absolute value) of stepping length
    Double minEta_;

    // stoppingTol_: tolerance up to which a given criterion has to be satisfied (e.g., value of residual)
    Double stoppingTol_;

    // measurePerformance_: shall runtime be measured?
    bool measurePerformance_;

    // allowNegativeEta_: if true, eta might be negative as well, i.e., search direction is reverted
    // should not be possible in backtracking case; only for exact linesearches (GoldenSection, QuadraticPolynomial)
    bool allowNegativeEta_;

    /*
     * Specific parameter
     */
    // LS_TRIAL_AND_ERROR (backtracking linesearch; stopping criterion = residual)
    // etaStart_: starting value for eta (default = 1.0)
    Double etaStart_;

    // decreaseFactor_: factor by which eta shall be decreased (default = 0.5)
    Double decreaseFactor_;

    // residualDecreaseForSuccess_: additional stopping criterion; if residual decreased by this factor compared
    //  to initial residual, we stop
    Double residualDecreaseForSuccess_;

    // If iteration stopped due to (b; residual increased twice in a row) and the currently best
    // found eta is etaStart, we can check if larger eta (i.e. overrelaxation) would lead to a
    // better value for eta; in this case 1.0/decreaseFactor will be utilized as increase factor;
    bool testOverrelaxation_;

    // LS_BACKTRACKING_ARMIJO, LS_BACKTRACKING_ARMIJO_NONMONOTONIC, LS_BACKTRACKING_GRADFREE (backtracking linesearch; stopping criterion = Armijo-type rule)
    // (non-monotonic) Armijo-rule: f(x + eta*dx) <= f_max_M(x) + rho*eta*gradF(x)^T*dx
    // etaStart and decreaseFactor as above;
    // rho_: parameter as in Armijo-rule; default 0.9
    Double rho_;

    // historyLengthArmijo_: number of entries from which f_max_M(x) is determined
    //  value = 1: just compare to previous value
    //  value = M: compare to maximum of previous M iterations
    UInt historyLengthArmijo_;

    // LS_THREEPOINTPOLY
    // stop if currentResidual < (1-alphaCheck*etaCurrent)*initialResidual
    // alphaCheck = 0.001
    Double alphaCheck_;

    // safeguarding parameter > limit stepping length
    //    if(eta < sigma0*etaCurrent){
    //      eta = sigma0*etaCurrent;
    //    } else if(eta > sigma1*etaCurrent){
    //      eta = sigma1*etaCurrent;
    //    }
    Double sigma0_; //sigma0 = 0.1;
    Double sigma1_; //sigma1 = 0.5;

    // LS_VECBASEDPOLY
    /* determine new eta based on h(x)
     *  1             if x >= 0.5 + epsilon => stopping criterion
     *  x             if epsilon <= x < 0.5 + epsilon
     *  epsilon       if 0 <= x < epsilon
     * -epsilon       if -epsilon <= x < 0
     * x              if -0.5-epsilon <= x < -epsilon
     * -0.5-epsilon   if x < -0.5-epsilon
     */
    // Default: epsilon = 0.16; // slightly smaller than 1/6
    Double epsilon_;

    // LS_GOLDENSECTION, LS_QUADINTERPOL
    // no specific parameter required; starting value for eta is determined by inspecting initial search interval
  };


  struct StoppingCriterion {
    // checkingOrder_: as name suggests, order in which multiple criteria shall be checked
    UInt checkingOrder_;

    bool operator<(const StoppingCriterion& a) const
    {
        return checkingOrder_ < a.checkingOrder_;
    }

    // nameTag_: for human readable output later on
    std::string nameTag_;

    // isResidual_: true > check residual; false > check update/increment of unknown
    bool isResidual_;

    // isRelative_: use relative value if true; else use absolute value
    bool isRelative_;

    // scalingFactor_: may be used to scale residual e.g., to the size of the residual
    bool scaleToSystemSize_;
    Double scalingFactor_;

    // stoppingTolerance_: actual stopping tolerance
    Double stoppingTolerance_;

    // firstCheck_: first iteration at which criterion shall be checked (for failback criteria which shall be checked
    //  after e.g., the 10th iteration for the first time
    UInt firstCheck_;

    // checkingFrequency_: number of iterations that has to be passed before criterion is checked again;
    // e.g. 1 = check every iteration
    //      2 = check every second iteration
    UInt checkingFrequency_;

    // lastCheckedValue_: optional storage for last checked value;
    mutable Double lastCheckedValue_;
  };

  //! all kind or parameter related to Preisach weight function
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
    
    // weightType = 4 and 5 > Lorentzian function
    Double muLorentz_A_;
    Double muLorentz_h1_;
    Double muLorentz_sigma1_;
    Double muLorentz_h2_;
    Double muLorentz_sigma2_;

    // anhysteretic parameter
    Double anhysteretic_a_;
    Double anhysteretic_b_;
    Double anhysteretic_c_;
    Double anhysteretic_d_;
    Double anhystAtSat_normalized_;
//    int anhysteretic_cInAtan_;
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

  //! Parameter for magnetostrictive and piezoelectric coupling in context of hysteresis
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

  //! Parameter describing the different Preisach hysteresis models
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
    
    /*
     * Added specifically for revised version
     * here, we get a guaranteed alignment of input and output if abs(input) > inputSat_/kappa_rev;
     * we can further have alignment if we enforce it by hand via mat.xml by setting flag enforceSatOutputAtSatInput; 
     * note: for kappa_rev >= 1, alignment is not necessarily achieved by inputSat_/kappa_rev as the
     * angular lag deltaPhi first vanishes at inputSat_! so inputSat_ is the lower limit!
     * 
     * for classic version, we have guaranteed alignment of input and output if abs(input) > inputSat_
     * for Mayergoyz, there is no guaranteed alignment due to remanent contributions (those can be zero in
     * specific cases, i.e., if whole excitation history is uniaxial, but ingeneral no alignment can be ensured)
     * 
     * Consequence: 
     *  Classic: inputForAlignment_ = inputSat_
     *  Revised: inputForAlignment_ = inputSat_/kappa_rev for kappa_rev < 1; inputSat for kappa_rev >= 1 or if enforceSatOutputAtSatInput=true
     *  Mayergoyz: inputForAlignment_ = numeric_limits<Double>::max()
     * 
     * inputForAlignment_ can be used in Hysteresis.cc to check whether we are in the region where
     * input and output align ABOVE inputSat_ so that we can utilize simple bisection for anhyst part
     */
    Double inputForAlignment_;
    
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
    Vector<Double> startingAxisMG_;
    // additional parameter for fitting rotational loss
    // source: Dlala - "Improving Loss Properties of the Mayergoyz Vector Hysteresis Model"
    Double lossParam_a;
    Double lossParam_b;
    int restrictionOfPsi_;
    bool useAbsoluteValueOfdPhi_;
    bool normalizeXInExponentOfG_;
    Double scalingFactorXInExponent_;

    // debugging option
    bool checkInversionResult_;
    bool printWarnings_;
  };

  //! Parameter for testing various aspects of hysteresis like evaluation, inversion, approximation of Jacobian etc.
  //! will be read from input.xml in SolveStepHyst
  struct TestOptions {
    bool testForwardEvaluation;
    bool testInversion;

    bool testJacobianApproximations;
    bool traceHysteresisOperator;
  };

} // end of CoupledField

#endif
