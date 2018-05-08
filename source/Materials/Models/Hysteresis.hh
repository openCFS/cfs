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
    
    virtual void EvaluateRotationState(UInt idx){
      EXCEPTION( "EvaluateRotationState not implemented in base-Class" );
    }
    
    virtual Vector<Double> getRotationDirectionAndUpdate(Vector<Double> flux_in, Matrix<Double> eps_mu, UInt idx){
      EXCEPTION( "getRotationDirectionAndUpdate not implemented in base-Class" );
    }
    
    //!
    virtual Double computeValue(Double& xVal, Integer idxElem, bool overwrite = true) {
      EXCEPTION( "computeValue not implemented in base-Class" );
      return 0.0;
    };

    virtual Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idxElem, bool overwrite = true,bool overwriteDirection= true,bool debugOut = false) {
      EXCEPTION( "computeValue_vec not implemented in base-Class" );
      Vector<Double> Yout;
      return Yout;
    };

    virtual Vector<Double> computeInput_vec(Vector<Double> yVal, Integer operatorIndex, 
      Matrix<Double> mu, bool overwriteDirection = true){
      EXCEPTION("computeInput_vec not implemented in base-class");
    }
    
    virtual Vector<Double> computeInput_vec_withStatistics(Vector<Double> yVal, Vector<Double> prevYval,
      Vector<Double> prevXval, Vector<Double> prevHystval, Integer operatorIndex, 
      Matrix<Double> mu, bool overwriteDirection, 
      UInt& totalNumberOfLMIterations, UInt& totalNumberOfLinesearchIterations, 
      UInt& maximalNumberOfLinesearchIterations, UInt& successCode, 
      Double& minAlpha, Double& maxAlpha, Double& avgAlpha ){
      EXCEPTION("computeInput_vec_withStatistics not implemented in base-class");
    }
        
    // do not overwrite memory per default
    virtual Double computeInputAndUpdate(Double Yin, Double eps_mu, Integer operatorIndex, bool overwrite = false){
      EXCEPTION( "computeInputAndUpdate not implemented in base-Class" );
      return 0.0;
    };
    
    //!
    virtual Double computeValueAndUpdate(Double xVal, Integer idxElem, 
                                         bool overwrite = true ) {
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
                                  bool overwrite) {
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
    
    virtual void SetParamsForInversion(UInt maxIter, Double resTolH, Double resTolB, Double jacobiResolution,
          bool useTikhonov, Double alphaLSStart, Double alphaLSMin, Double alphaLSMax, Double angClipping){
      EXCEPTION( "Only implemented and required for VectorPreisach model");
    };

    inline Double evalAnhystPart_normalized(Double xNormalizedUnclipped){
      // returns normalized anhysteretic part
      return anhyst_A_*std::atan(anhyst_B_*xNormalizedUnclipped) + anhyst_C_*xNormalizedUnclipped;
    }
    
    Double bisectForAnhyst_normalized(Double Ytarget_normalized, 
      Double Xdown_normalized, Double Xup_normalized, Double Poffset_normalized, Double eps_mu_normalized, Double tol);
    
    Double bisectForAnhyst(Double Ytarget, Double Xdown, Double Xup, Double Poffset, Double eps_mu, Double tol);
    
  protected:
    Double anhyst_A_;
    Double anhyst_B_;
    Double anhyst_C_;
    bool anhystOnly_;
    Double XSaturated_;
    Double YSaturated_;
  private:

    Integer numElements_;

  };


} //end of namespace


#endif

