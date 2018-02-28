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

    ///
    virtual ~Hysteresis();

    //!
    virtual Double computeValue(Double& xVal, Integer idxElem, bool overwrite = true) {
      EXCEPTION( "computeValue not implemented in base-Class" );
      return 0.0;
    };

    virtual Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idxElem, bool overwrite = true,bool overwriteDirection= true) {
      EXCEPTION( "computeValue not implemented in base-Class" );
      Vector<Double> Yout;
      return Yout;
    };

    virtual Vector<Double> computeInput_vec(Vector<Double>& yVal, Integer idElem, Matrix<Double> mu, Double& alpha, bool overwrite = true,bool overwriteDirection = true) {
      EXCEPTION("computeInput_vec not implemented in base-class");
    };
    
    virtual Double computeInputAndUpdate(Double Yin, Double eps_mu, Integer idx, bool overwrite = true){
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

    virtual void setFlags(UInt performanceFlag, UInt textOutputFlag, UInt mappingFlag){
      EXCEPTION( "setFlags not inplemented in base-Class")
    };


  protected:

  private:

    Integer numElements_;
  };


} //end of namespace


#endif

