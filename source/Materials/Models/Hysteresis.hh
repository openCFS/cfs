#ifndef FILE_HYSTERESIS
#define FILE_HYSTERESIS

#include <string>
#include "General/defs.hh"
#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include "General/Exception.hh"
#include "Utils/tools.hh"

namespace CoupledField {

  class Hysteresis
  {
  public:
    Hysteresis(Integer numElem);
    
    Hysteresis(Integer numElem, Double XSaturated, Double YSaturated, Double anhystA, Double anhystB, Double anhystC, bool anhystOnly);

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
      return anhyst_A_*std::atan(anhyst_B_*xNormalizedUnclipped) + anhyst_C_*xNormalizedUnclipped;
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
    void SetParamsForInversion(UInt maxIter, Double resTolH, Double resTolB, Double jacobiResolution,
         bool useTikhonov, Double alphaLSStart, Double alphaLSMin, Double alphaLSMax, Double angClipping){
      INV_maxIter_ = maxIter;
      INV_resTolH_ = resTolH;
      INV_resTolB_ = resTolB;
      INV_jacobiResolution_ = jacobiResolution;
      INV_useTikhonov_ = useTikhonov;
      INV_alphaLSStart_ = alphaLSStart;
      INV_alphaLSMin_ = alphaLSMin;
      INV_alphaLSMax_ = alphaLSMax;
      INV_angClipping_ = angClipping;
    }

    bool checkInversionOutput(Vector<Double>& xComputed, Vector<Double>& yTarget, 
			Matrix<Double>& mu, Double tol, Double& resYNorm, Integer operatorIdx, bool overwriteMemory, bool output = false);
		
    bool checkConvergence(Vector<Double>& res, Matrix<Double>& jacT, Double& errorNorm, Double tol);
    
//    Double computeRho(Vector<Double>& xNew, Vector<Double>& xUpdate, 
//				Vector<Double>& res, Vector<Double>& resShifted, Matrix<Double>& jac);
//    
//    Integer checkIncrement(Vector<Double>& xNew, Vector<Double>& xUpdate, 
//		Vector<Double>& res, Vector<Double>& resShifted, Matrix<Double>& jac, Double& alpha, int stayBelowSat);
//    
    Integer checkIncrementTrustRegion(Vector<Double>& x_new, 
          Vector<Double>& res_cur, Vector<Double>& res_new,
          Vector<Double>& jac_dx, Double& alpha, int stayBelowSat);
    
//    Integer checkIncrementOLD(Vector<Double>& xNew, Vector<Double>& xUpdate, 
//		Vector<Double>& res, Vector<Double>& resShifted, Matrix<Double>& jac, Double& alpha, int stayBelowSat);
//    
    Vector<Double> computeResidual(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& hystVal, Matrix<Double> mu_inv);
        
    Matrix<Double> computeJacobian(Vector<Double>& xVal, Vector<Double>& hystVal, 
          Matrix<Double> mu_inv, Integer operatorIdx, Double sign, UInt implementation, 
					bool overwriteMemory, int stayBelowSat);
    
    Vector<Double> computeUpdate_Newton(Vector<Double>& x, Vector<Double>& y, 
      Vector<Double>& hyst_x, Matrix<Double> mu_inv, Integer operatorIdx);
    
    Vector<Double> computeUpdate_Krylov(Vector<Double>& x, Vector<Double>& y, 
      Vector<Double>& hyst_x, Matrix<Double> mu_inv, Integer operatorIdx);
    
    Vector<Double> computeUpdate_LM(Vector<Double> jacTres_neg, 
          Matrix<Double>& jacTjac, Double& alphaIn, Double& alphaAcc, Double& alphaMinReq, Double alphaMax);
    
    Vector<Double> computeJacobianTimesVector(Vector<Double>& x, Vector<Double>& v, 
      Vector<Double>& y, Vector<Double>& hyst_x, Matrix<Double> mu_inv, Integer operatorIdx);
      
    bool computeUpdate_LM_full(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& res, 
      Vector<Double>& xUpdate, Matrix<Double>& jac, Matrix<Double>& jacT, Matrix<Double> mu, Matrix<Double> mu_inv, 
      Integer operatorIdx, bool overwriteMemory,
      Double& alpha, Double alphaMin, Double alphaMax,  
      UInt& numberOfIterations, Vector<Double>& xStart, Double factorToSat,int stayBelowSat,Vector<Double> sol);

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
    bool anhystOnly_;
    Double XSaturated_;
    Double PSaturated_;
    UInt dim_;
    
    // additional for inversion
    Vector<Double>* prevXVal_;
    Vector<Double>* prevHVal_;
    
    // For inversion via Levenberg Marquardt
    UInt INV_maxIter_;
    Double INV_resTolH_;
    Double INV_resTolB_;
    Double INV_jacobiResolution_;
    bool INV_useTikhonov_;
    Double INV_alphaLSStart_;
    Double INV_alphaLSMin_;
    Double INV_alphaLSMax_;
    Double INV_angClipping_;
  private:

    Integer numElements_;

  };


} //end of namespace


#endif

