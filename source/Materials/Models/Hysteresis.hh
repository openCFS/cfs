#ifndef FILE_HYSTERESIS
#define FILE_HYSTERESIS

#include <string>
#include "General/defs.hh"
#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include "General/Exception.hh"
#include "Utils/tools.hh"
#include "Utils/helperStructs.hh"

namespace CoupledField {

  class Hysteresis
  {
  public:
    Hysteresis(Integer numElem);
    
    Hysteresis(Integer numElem, ParameterPreisachOperators operatorParams, ParameterPreisachWeights weightParams);
    
    Hysteresis(Integer numElem, Double XSaturated, Double YSaturated, Double hystSaturated, 
    Double anhystA, Double anhystB, Double anhystC, bool anhystOnly, bool anhyst_cInAtan);

    ///
    virtual ~Hysteresis();

    // for extension of scalar model
    virtual void UpdateRotationState(Vector<Double> flux_in, Matrix<Double> eps_mu, UInt idx){
      EXCEPTION( "UpdateRotationState not implemented in base-Class" );
    }
    
    virtual void UpdateRotationStateWithFluxDensity(Vector<Double> flux_in, Matrix<Double> eps_mu, UInt idx){
      EXCEPTION( "UpdateRotationStateWithFluxDensity not implemented in base-Class" );
    }
    virtual void UpdateRotationStateWithFieldIntensity(Vector<Double> field_in, UInt idx){
      EXCEPTION( "UpdateRotationStateWithFieldIntensity not implemented in base-Class" );
    }
    
    virtual void EvaluateRotationState(UInt idx){
      EXCEPTION( "EvaluateRotationState not implemented in base-Class" );
    }
        
    virtual Vector<Double> getRotationDirection(UInt idx){
      EXCEPTION( "getRotationDirection not implemented in base-Class" );
    }
    
    //!
    virtual Double computeValue(Double& xVal, Integer idxElem, bool overwrite, int& successFlag) {
      EXCEPTION( "computeValue not implemented in base-Class" );
      return 0.0;
    };

    virtual Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idxElem, bool overwrite,
      bool debugOut, int& successFlag) {
      EXCEPTION( "computeValue_vec not implemented in base-Class" );
      Vector<Double> Yout;
      return Yout;
    };
    
    virtual Vector<Double> computeValue_vecMeasure(Vector<Double>& xVal, Integer idxElem, bool overwrite,
      bool debugOut, int& successFlag, Double& time) {
      EXCEPTION( "computeValue_vec not implemented in base-Class" );
      Vector<Double> Yout;
      return Yout;
    };

    virtual Vector<Double> computeInput_vec(Vector<Double> yVal, Integer operatorIndex, 
      Matrix<Double> mu, bool fieldsAlignedAboveSat, bool hystOutputRestrictedToSat, 
      int& successFlag){
      EXCEPTION("computeInput_vec not implemented in base-class");
    }
    
//    virtual Vector<Double> computeInput_vec_withStatistics(Vector<Double> yVal, Vector<Double> prevYval,
//      Vector<Double> prevXval, Vector<Double> prevHystval, Integer operatorIndex, 
//      Matrix<Double> mu, bool overwriteDirection, 
//      UInt& totalNumberOfLMIterations, UInt& totalNumberOfLinesearchIterations, 
//      UInt& maximalNumberOfLinesearchIterations, UInt& successCode, 
//      Double& minAlpha, Double& maxAlpha, Double& avgAlpha ){
//      EXCEPTION("computeInput_vec_withStatistics not implemented in base-class");
//    }
        
    // do not overwrite memory per default
    virtual Double computeInputAndUpdate(Double Yin, Double eps_mu, Integer operatorIndex, 
      bool overwrite, int& successFlag){
      EXCEPTION( "computeInputAndUpdate not implemented in base-Class" );
      return 0.0;
    };
    
    //!
    virtual Double computeValueAndUpdate(Double xVal, Integer idxElem, 
      bool overwrite, int& successFlag ) {
      EXCEPTION( "computeValueAndUpdate not implemented in base-Class" );
      return 0.0;
    };
    
    virtual Double computeValueAndUpdateMeasure( Double Xin, Integer idx,
          bool overwrite, int& successFlag, Double& time ){
      EXCEPTION( "computeValueAndUpdate not implemented in base-Class" );
      return 0.0;
    };

    virtual void SetPreviousYval( Double yval, Integer idxElem )  {
      EXCEPTION( "SetPreviousYval not implemented in base-Class" );
    };

    //! 
    virtual Double getValue(  Integer idxElem) {
      EXCEPTION( "getValue not implemented in base-Class" );
      return 0.0;
    }

    virtual Double getActXval ( UInt idxElem ) {
      EXCEPTION( "getActXval not implemented in base-Class" );
      return 0.0;
    }

    //! 
    virtual UInt getStringLength( Integer idxElem ) {
      EXCEPTION( "getStringLength not implemented in base-Class" );
      return 0;
    }

    //!
    virtual Double GetIncX() {
      EXCEPTION( "GetIncX not implemented in base-Class" );
      return 1.0;
    };

    //!
    virtual Double updateMinMaxList(Double newX, Integer idxElem, 
      bool overwrite, int& successFlag) {
      EXCEPTION( "updateMinMaxList not implemented in base-Class" );
      return 0;
    };

    //! 
    virtual void SetTimeStepVal(Double dt) {
      EXCEPTION( "SetTimeStepVal not implemented in base-Class" );
    };

    //!
    virtual Double EvalEverett(Double x1, Double x2, Integer idx ) {
      EXCEPTION( " EvalEverett not implemented in base-Class" );
      return 0.0;
    };

    virtual void switchingStateToBmp(UInt numPixel, std::string filename, UInt idElem,bool overLayWithRotState = false){
      EXCEPTION( " switchingStateToBmp not implemented in base-Class" );
    };

    virtual void rotationStateToBmp(UInt numPixel, std::string filename, UInt idElem){
      EXCEPTION( " switchingStateToBmp not implemented in base-Class" );
    };

    virtual std::string runtimeToString(){
      EXCEPTION( "runtimeToString not implemented in base-Class" );
    };

    virtual void setFlags(UInt performanceFlag){
      EXCEPTION( "setFlag not inplemented in base-Class")
    };
    
//    virtual void SetParamsForInversion(UInt maxIter, Double resTolH, Double resTolB, Double jacobiResolution,
//          bool useTikhonov, Double alphaLSStart, Double alphaLSMin, Double alphaLSMax, Double angClipping){
//      EXCEPTION( "Only implemented and required for VectorPreisach model");
//    };

    inline Double evalAnhystPart_normalized(Double xNormalizedUnclipped){
      // returns normalized anhysteretic part
      // note: at this point we assume that a,b,c are meant for p,e in range [-1,1]
      //       the script for determining these parameters assuems p,e in range [-0.5,0.5] 
      //      > parameters have to be transfereed to correct range
      //         prior to calling this function; this can either be done by user or by CFS if flag in mat.xml is set
      if(anhyst_cInAtan_ == true){
        return anhyst_A_*std::atan(anhyst_B_*(xNormalizedUnclipped + anhyst_C_));
      } else {
        return anhyst_A_*std::atan(anhyst_B_*xNormalizedUnclipped) + anhyst_C_*xNormalizedUnclipped;
      }
    }
    
    static Double evalAnhystPart_normalized_atSaturation(Double anhystA, Double anhystB, Double anhystC, bool cInAtan){
      // anhyst parameter are assumed to be defined for full range i.e. from -1 to +1
//      std::cout << "anhystA: " << anhystA << std::endl;
//      std::cout << "anhystB: " << anhystB << std::endl;
//      std::cout << "anhystC: " << anhystC << std::endl;
      
      if(cInAtan){
        return anhystA*std::atan(anhystB + anhystC);
      } else {
        return anhystA*std::atan(anhystB) + anhystC;
      }
    }
    
    
    Double bisectForAnhyst_normalized(Double Ytarget_normalized, Double Xdown_normalized, Double Xup_normalized, 
      Double Poffset_normalized, Double eps_mu_normalized, Double tol, Vector<Double> dir, UInt idx);
    
    Double bisectForAnhyst(Double Ytarget, Double Xdown, Double Xup, Double Poffset, Double eps_mu, Double tol){
      Vector<Double> zeroVec = Vector<Double>(dim_);
      zeroVec.Init();
      return bisectForAnhyst(Ytarget, Xdown, Xup, Poffset, eps_mu, tol, zeroVec, 0);
    }
    
    Double bisectForAnhyst(Double Ytarget, Double Xdown, Double Xup, Double Poffset, Double eps_mu, Double tol, Vector<Double> dir, UInt idx);
    
    // from VecPreisachv10 > put into baseclass to make it available for Mayergoyz model, too
    void SetParamsForInversion(ParameterInversion inversionParams){
          
      INV_params_ = inversionParams;
//      std::cout << "INV_params_.inversionMethod " << INV_params_.inversionMethod << std::endl;
//      std::cout << "INV_params_.maxNumIts " << INV_params_.maxNumIts << std::endl;
//      std::cout << "INV_params_.tolH " << INV_params_.tolH << std::endl;
//      std::cout << "INV_params_.tolB " << INV_params_.tolB << std::endl;
//      std::cout << "INV_params_.maxNumRegIts " << INV_params_.maxNumRegIts << std::endl;
//      std::cout << "INV_params_.alphaRegStart " << INV_params_.alphaRegStart << std::endl;
//      std::cout << "INV_params_.alphaRegMin " << INV_params_.alphaRegMin << std::endl;
//      std::cout << "INV_params_.alphaRegMax " << INV_params_.alphaRegMax << std::endl;
//      std::cout << "INV_params_.trustLow " << INV_params_.trustLow << std::endl;
//      std::cout << "INV_params_.trustMid " << INV_params_.trustMid << std::endl;
//      std::cout << "INV_params_.trustHigh " << INV_params_.trustHigh << std::endl;
//      std::cout << "INV_params_.maxNumLSIts " << INV_params_.maxNumLSIts << std::endl;
//      std::cout << "INV_params_.alphaLSMin " << INV_params_.alphaLSMin << std::endl;
//      std::cout << "INV_params_.alphaLSMax " << INV_params_.alphaLSMax << std::endl;
//      std::cout << "INV_params_.jacRes " << INV_params_.jacRes << std::endl;
//      std::cout << "INV_params_.jacImplementation " << INV_params_.jacImplementation << std::endl;
//      std::cout << "INV_params_.stopLineSearchAtLocalMin " << INV_params_.stopLineSearchAtLocalMin << std::endl;    
//      
//      std::cout << "INV_params_.projLM_mu " << INV_params_.projLM_mu << std::endl;
//      std::cout << "INV_params_.projLM_rho " << INV_params_.projLM_rho << std::endl;
//      std::cout << "INV_params_.projLM_beta " << INV_params_.projLM_beta << std::endl;
//      std::cout << "INV_params_.projLM_sigma " << INV_params_.projLM_sigma << std::endl;
//      std::cout << "INV_params_.projLM_gamma " << INV_params_.projLM_gamma << std::endl;
//      std::cout << "INV_params_.projLM_tau " << INV_params_.projLM_tau << std::endl;
//      std::cout << "INV_params_.projLM_c " << INV_params_.projLM_c << std::endl;  
//      std::cout << "INV_params_.projLM_p " << INV_params_.projLM_p << std::endl;  
      
    }

    bool checkInversionOutput(Vector<Double>& xComputed, Vector<Double>& yTarget, 
			Matrix<Double>& mu, Double tol, Double& resYNorm, Integer operatorIdx, bool overwriteMemory, bool output = false);
		
    bool checkConvergence(Vector<Double>& jacTres, Double& errorNorm, Double tol);
     
    Integer checkIncrementTrustRegion(Vector<Double>& x_new, 
          Vector<Double>& res_cur, Vector<Double>& res_new,
          Vector<Double>& jac_dx, Double& alpha, Double alphaStepUp, Double alphaStepDown, int stayBelowSat);
    
    Vector<Double> computeResidual(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& hystVal, Matrix<Double> mu_inv);
        
    Matrix<Double> computeJacobian(Vector<Double>& xVal, Vector<Double>& hystVal, 
          Matrix<Double> mu_inv, Integer operatorIdx, Double sign, 
					bool overwriteMemory, int stayBelowSat, Double scalingForJacDiagonal);
    
    Vector<Double> computeJacobianTimesVector(Vector<Double>& x, Vector<Double>& v, 
    Vector<Double>& y, Vector<Double>& hyst_x, Matrix<Double> mu_inv, Integer operatorIdx);
        
    Vector<Double> computeUpdate_Newton(Vector<Double>& x, Vector<Double>& y, 
      Vector<Double>& hyst_x, Matrix<Double> mu_inv, Integer operatorIdx, Double scalingForJacDiagonal);
    
    Vector<Double> computeUpdate_Krylov(Vector<Double>& x, Vector<Double>& y, 
      Vector<Double>& hyst_x, Matrix<Double> mu_inv, Integer operatorIdx);
    
    Vector<Double> computeUpdate_LM_withFixedAlpha(Vector<Double> jacTres_neg, 
          Matrix<Double>& jacTjac, Vector<Double>& res, Double alphaFix);
    
    Vector<Double> computeUpdate_LM(Vector<Double> jacTres_neg, 
          Matrix<Double>& jacTjac, Double& alphaIn, Double& alphaAcc, Double& alphaMinReq, Double alphaMax);
    
    Vector<Double> projectToSolutionSpace(Vector<Double> x, int stayBelowSat);
    
    bool computeUpdateLM_Projected(Vector<Double>& xk, Vector<Double>& resk, Vector<Double>& y,
      Matrix<Double>& mu_inv, Matrix<Double>& jac, Vector<Double>& jacTres,
      Vector<Double>& dx, UInt operatorIdx, int stayBelowSat);

    bool computeUpdateLM_TrustRegion(Vector<Double>& xStart, Vector<Double>& xCurrent, Vector<Double>& xUpdate, 
      Vector<Double>& hystCurrent, Vector<Double>& resCurrent, Vector<Double>& yTarget, 
      Matrix<Double>& mu_inv, Matrix<Double>& jacCurrent, Vector<Double>& jacTresCurrent, 
      int operatorIdx, int stayBelowSat, 
      Double& alpha, Double alphaMin, Double alphaMax, UInt& numberOfIterations);

    bool computeUpdateNewton_Linesearch(Vector<Double>& xStart, Vector<Double>& xCurrent, Vector<Double>& xUpdate, 
      Vector<Double>& hystCurrent, Vector<Double>& resCurrent, Vector<Double>& yTarget, 
      Matrix<Double>& mu_inv, Matrix<Double>& jacCurrent, Vector<Double>& jacTresCurrent, 
      int operatorIdx, int stayBelowSat,
      Double& alpha, Double& alphaMin, Double& alphaMax, bool stopLineSearchAtLocalMin,
      Double scalingForJacDiagonal);
  
    bool computeUpdate(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& res, 
      Vector<Double>& xUpdate, Matrix<Double>& jac, Matrix<Double>& jacT, Matrix<Double> mu, Matrix<Double> mu_inv, 
      Integer operatorIdx, bool overwriteMemory,
      Double& alpha, Double alphaMin, Double alphaMax,  
      UInt& numberOfIterations, Vector<Double>& xStart,int stayBelowSat,Vector<Double> sol);

    Vector<Double> computeInput_vec_withPrevStates(Vector<Double> yVal, Vector<Double> prevYval,
      Vector<Double> prevXval, Vector<Double> prevHystval, Integer operatorIndex, 
      Matrix<Double> mu, bool fieldsAlignedAboveSat, bool hystOutputRestrictedToSat, int& successFlag);
    
    Vector<Double> computeInput_vec_withStatistics(Vector<Double> yVal, Vector<Double> prevYval,
      Vector<Double> prevXval, Vector<Double> prevHystval, Integer operatorIndex, 
      Matrix<Double> mu, bool fieldsAlignedAboveSat, bool hystOutputRestrictedToSat,
      UInt& totalNumberOfLMIterations, UInt& totalNumberOfLinesearchIterations, 
      UInt& maximalNumberOfLinesearchIterations, int& successFlag, Double& minAlpha, 
      Double& maxAlpha, Double& avgAlpha,Vector<Double> sol );
    
    Vector<Double> obtainStartVector(Vector<Double> previousXVector, Vector<Double> previousYVector,
      Vector<Double> currentYVector, Double YdiffToRemancence, Double YdiffToSaturation, int stayBelowSat);
    
  protected:
    Double anhyst_A_;
    Double anhyst_B_;
    Double anhyst_C_;
    bool anhyst_cInAtan_;
    bool anhystOnly_;
    Double XSaturated_;
    Double PSaturated_;
    Double hystSaturated_; // saturation value of hyst operator alone (without anhyst part); can be the same as PSaturated
    // depending on parameter anhystCountingToOutputSat_ (part of weight params)
    UInt dim_;
    
    // additional for inversion
    Vector<Double>* prevXVal_;
    Vector<Double>* prevHVal_;
    
    // For vector inversion
    ParameterInversion INV_params_;
    
  private:

    Integer numElements_;

  };


} //end of namespace


#endif

