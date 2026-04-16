#ifndef FILE_HYSTERESIS
#define FILE_HYSTERESIS

#include <string>
#include "General/defs.hh"
#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include "General/Exception.hh"
#include "Utils/ToolsFull.hh"
#include "Utils/helperStructs.hh"
#include "Utils/Timer.hh"

namespace CoupledField {

  class Hysteresis
  {
  public:
    /*
     * Base routines I - Constructors, Destructor
     */
    Hysteresis(Integer numElem);

    Hysteresis(Integer numElem, ParameterPreisachOperators operatorParams, ParameterPreisachWeights weightParams);

    Hysteresis(Integer numElem, Double XSaturated, Double YSaturated, Double hystSaturated,
    Double anhystA, Double anhystB, Double anhystC, Double anhystD, bool anhystOnly);

    virtual ~Hysteresis();

    /*
     * Base routines II - Forward evaluation of operator
     */
    // Scalar Model
    virtual Double computeValue(Double& xVal, Integer idxElem, bool overwrite, int& successFlag) {
      EXCEPTION( "computeValue not implemented in base-Class" );
      return 0.0;
    };

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

    // Vector Model
    virtual Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idxElem, bool overwrite,
      bool debugOut, int& successFlag, bool skipAnhystPart = false) {
      EXCEPTION( "computeValue_vec not implemented in base-Class" );
      Vector<Double> Yout;
      return Yout;
    };
    
    Vector<Double> computeValue_vecMeasure(Vector<Double>& xVal, Integer idx,
    bool overwrite,bool debugOutput,int& successFlag, Double& time){

      Timer* timer = new Timer();
      Double startTime = timer->GetCPUTime();
      timer->Start();

      Vector<Double> Yvec = computeValue_vec(xVal, idx, overwrite, debugOutput, successFlag);

      timer->Stop();
      Double endTime = timer->GetCPUTime();
      time = endTime-startTime;

      return Yvec;
    };

    /*
     * Base routines III - Backward evaluation of operator / Inversion
     */
    // Scalar Model
    virtual Double computeInputAndUpdate(Double Yin, Double eps_mu, Integer operatorIndex,
      bool overwrite, int& successFlag){
      EXCEPTION( "computeInputAndUpdate not implemented in base-Class" );
      return 0.0;
    };
    
    // Vector Model
    virtual Vector<Double> computeInput_vec(Vector<Double> yVal, Integer operatorIndex,
      Matrix<Double> mu, bool fieldsAlignedAboveSat, bool hystOutputRestrictedToSat,
      int& successFlag, bool useEverett=false, bool overwrite=false){
      EXCEPTION("computeInput_vec not implemented in base-class");
    }

    virtual Vector<Double> computeInput_vec_Everett(Vector<Double> yVal, Integer operatorIndex,
          Matrix<Double> mu, bool overwrite, int& successFlag){
      EXCEPTION("computeInput_vec_Everett only implemented for Mayergoyz vector model");
    }

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

    virtual void collectParallelProjections(bool switchOnOff, UInt cnt){   
      EXCEPTION( " collectParallelProjections not implemented in base-Class" );
    };
    
    virtual void rotationListToTxt(std::string filename, UInt idElem, bool append, std::string optionalHeader){
      EXCEPTION( " rotationListToTxt not implemented in base-Class" );
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
//      if(anhyst_cInAtan_ == true){
//        return anhyst_A_*std::atan(anhyst_B_*(xNormalizedUnclipped + anhyst_C_));
//      } else {
        return anhyst_A_*std::atan(anhyst_B_*xNormalizedUnclipped + anhyst_D_) + anhyst_C_*xNormalizedUnclipped;
//      }
    }

    static Double evalAnhystPart_normalized_atSaturation(Double anhystA, Double anhystB, Double anhystC, Double anhystD){
      // anhyst parameter are assumed to be defined for full range i.e. from -1 to +1
//      std::cout << "anhystA: " << anhystA << std::endl;
//      std::cout << "anhystB: " << anhystB << std::endl;
//      std::cout << "anhystC: " << anhystC << std::endl;

//      if(cInAtan){
//        return anhystA*std::atan(anhystB + anhystC);
//      } else {
        return anhystA*std::atan(anhystB + anhystD) + anhystC;
//      }
    }

    Double bisectForAnhyst_normalized(Double Ytarget_normalized, Double Xdown_normalized, Double Xup_normalized,
      Double Poffset_normalized, Double eps_mu_normalized, Double tol, Vector<Double> dir, UInt idx, int& successFlag);

    Double bisectForAnhyst(Double Ytarget, Double Xdown, Double Xup, Double Poffset, Double eps_mu, Double tol, int& successFlag){
      Vector<Double> zeroVec = Vector<Double>(dim_);
      zeroVec.Init();
      return bisectForAnhyst(Ytarget, Xdown, Xup, Poffset, eps_mu, tol, zeroVec, 0, successFlag);
    }

    Double bisectForAnhyst(Double Ytarget, Double Xdown, Double Xup, Double Poffset, Double eps_mu, Double tol, Vector<Double> dir, UInt idx, int& successFlag);

    // from VecPreisachv10 > put into baseclass to make it available for Mayergoyz model, too
    void SetParamsForInversion(ParameterInversion inversionParams){
      INV_params_ = inversionParams;
      invParamsSet_ = true; 
    }

    bool checkInversionOutput(Vector<Double>& xComputed, Vector<Double>& yTarget, Vector<Double>& yObtained,
			Matrix<Double>& mu, Double tol, Double& resYNorm, Integer operatorIdx, bool overwriteMemory, bool output = false);

    bool checkConvergence(Vector<Double>& jacTres, Double& errorNorm, Double tol);

    Integer checkIncrementTrustRegion(Vector<Double>& x_new,
          Vector<Double>& res_cur, Vector<Double>& res_new,
          Vector<Double>& jac_dx, Double& alpha, Double alphaStepUp, Double alphaStepDown, int stayBelowSat);

    Vector<Double> computeResidual(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& hystVal, Matrix<Double> mu_inv);

    Matrix<Double> computeJacobian(Vector<Double>& xVal, Vector<Double>& hystVal,
          Matrix<Double> mu_inv, Integer operatorIdx, Double sign,
					bool overwriteMemory, int stayBelowSat, Double scalingForJacDiagonal);

    // basically new version of computeJacobianTimesVector but
    // based on Kelly 1995
    Vector<Double> ApproximateDirectionalDerivative(Vector<Double>& x, Vector<Double>& v,
          Vector<Double>& y, Vector<Double>& hyst_x, Matrix<Double> mu_inv, Integer operatorIdx);
    Matrix<Double> ApproximateFDJacobian(Vector<Double>& x, Vector<Double>& y, Vector<Double>& hyst_x,
          Matrix<Double> mu_inv, Integer operatorIdx, int version);
    bool CompareJacobianApproximations(Vector<Double>& xVal, Vector<Double>& y,
          Matrix<Double> mu_inv, Integer operatorIdx);


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
      Double scalingForJacDiagonal, UInt& numberOfIterations);

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

    // for scalar preisach
    virtual void setFixDirection(Vector<Double> newDirection){
      EXCEPTION("Only for scalar model");
    }

    virtual Vector<Double> getFixDirection(){
      EXCEPTION("Only for scalar model");
    }
    
    // Helper functions
    static bool checkVectorEquality(Vector<Double>& e1, Vector<Double>& e2, UInt criterion, Double tol){

      if(e1.GetSize() != e2.GetSize()){
        return false;
      }

      if(criterion == 1){
        // check angle between vectors; tol will be the allowed angle in rad
        Double angle = 0.0;
        e1.Inner(e2,angle);
        angle = std::acos(angle);
        if(angle <= tol){
          return true;
        } else {
          return false;
        }
      }

      else if(criterion == 2){
        // check single entries; only if all entries are closer than tol to each other return true
        for(UInt i = 0; i < e1.GetSize(); i++){
          if(abs(e1[i]-e2[i]) > tol){
            return false;
          }
        }
        return true;
      }

      else {
        EXCEPTION("Unknown criterion for vector comparison.");
      }

    }


  protected:
    Double anhyst_A_;
    Double anhyst_B_;
    Double anhyst_C_;
    Double anhyst_D_;
//    bool anhyst_cInAtan_;
    bool anhystOnly_;
    Double XSaturated_;
    Double XForAlignment_;
    Double PSaturated_;
    /*
     * saturation value of hyst operator alone (without anhyst part); can be the same as PSaturated
     * depending on parameter anhystCountingToOutputSat_ (part of weight params)
     */
    Double hystSaturated_;
    UInt dim_;

    // Rrequired for inversion
    Vector<Double>* prevXVal_;
    Vector<Double>* prevHVal_;
    ParameterInversion INV_params_;
    bool invParamsSet_;

    bool checkInversionResult_;
    bool printWarnings_;

  private:
    Integer numElements_;

  };

} //end of namespace

#endif

