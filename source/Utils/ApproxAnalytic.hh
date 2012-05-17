// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_APPROXANALYTIC
#define FILE_APPROXANALYTIC

#include <string>

#include "ApproxData.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/vector.hh"
#include "Utils/mathParser/mathParser.hh"

namespace CoupledField {

  //! Base class for approximation of sampled data
  class ApproxAnalytic : public ApproxData
  {
  public:

    //! constructor
    ApproxAnalytic( std::string filename, MaterialType matType );

    //! destructor: nothing to do
    virtual ~ApproxAnalytic() {;};

    //! set the analytic expressions
    void SetAnalyticExpressions( std::string fncStr, std::string fncDerivStr );

    //! reads in the sampled data form file with name fncName
    void ReadNlinFunc( std::string fncName ) {;};

    //! perform checks
    void PerformChecksOnInputData() {;};

    //! performs the approximation
    virtual void CalcApproximation( bool start=true ) {;};

    //! set accuracy of measured data
    virtual void SetAccuracy( Double val ) {
      EXCEPTION(" ApproxAnalytic: SetAccuracy not implemented");
    };

    //! set maximal y-value
    virtual void SetMaxY( Double val ) {
      EXCEPTION(" ApproxAnalytic: SetMaxY not implemented");
    };

    //! evaluates the functions
    virtual Double EvaluateFunc( Double t ) {
      EXCEPTION(" ApproxAnalytic: EvaluateFunc not implemented");
    };

    //! evalutates the derivazive of the function
    virtual Double EvaluatePrime( Double t ) {
      EXCEPTION(" ApproxAnalytic: EvaluatePrime not implemented");
    };

    //! evalutates the invers of the function
    virtual Double EvaluateFuncInv( Double t ) {
      EXCEPTION(" ApproxAnalytic: EvaluateFuncInv not implemented");
    };


    //! valuates the inverse of the derivative of the function
    virtual Double EvaluatePrimeInv( Double t ) {
      EXCEPTION(" ApproxAnalytic: EvaluatePrimeInv not implemented");
    };


    //! computes the best approximation
    virtual void CalcBestParameter() {;};

    //! computes the magnetic reluctivity
    virtual Double EvaluateFuncNu(Double t);

    //! computes the derivative of magnetic reluctivity
    virtual Double EvaluatePrimeNu(Double t);

    //! prints out original and approximated function
    virtual void Print() {;};

  protected:

    //! Pointer to MathParser object
    MathParser * mParser_;

    //! math parser handle for function
    MathParser::HandleType mphFactor_;

    //! math parser handle for derivative of function
    MathParser::HandleType mphFactorDeriv_;

  };


} //end of namespace


#endif

