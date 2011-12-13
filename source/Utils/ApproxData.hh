// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ApproxData
#define FILE_ApproxData

#include <string>

#include "General/defs.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/vector.hh"

namespace CoupledField {

  //! Base class for approximation of sampled data
  class ApproxData
  {
  public:

    //! constructor
    ApproxData( std::string nlFncName, MaterialType matType );

    //! destructor: nothing to do
    virtual ~ApproxData() {;};

    //! reads in the sampled data form file with name fncName
    void ReadNlinFunc( std::string fncName );

    //! perform checks
    void PerformChecksOnInputData();

    //! performs the approximation
    virtual void CalcApproximation( bool start=true ) = 0;

    //! set accuracy of measured data
    virtual void SetAccuracy( Double val ) {
      EXCEPTION(" ApproxData: SetAccuracy not implemented");
    };

    //! set maximal y-value
    virtual void SetMaxY( Double val ) {
      EXCEPTION(" ApproxData: SetMaxY not implemented");
    };

    //! evaluates the functions
    virtual Double EvaluateFunc( Double t ) = 0;

    //! evalutates the derivazive of the function
    virtual Double EvaluatePrime( Double t ) = 0;

    //! evalutates the invers of the function
    virtual Double EvaluateFuncInv( Double t ) = 0;

    //! valuates the inverse of the derivative of the function
    virtual Double EvaluatePrimeInv( Double t ) = 0;

    //! computes the best approximation
    virtual void CalcBestParameter() = 0;

    //! computes the magnetic reluctivity
    virtual Double EvaluateFuncNu(Double t) {     
      EXCEPTION(" ApproxData: EvaluateFuncNu not implemented");
      return -1.0; 
    }

    //! computes the derivative of magnetic reluctivity
    virtual Double EvaluatePrimeNu(Double t) {
      EXCEPTION(" ApproxData: EvaluatePrimeNu not implemented");
      return -1.0; 
    }

    //! prints out original and approximated function
    virtual void Print() {;};

  protected:

    Vector<Double> x_;  //!< independent value
    Vector<Double> y_;  //!< function value
    UInt numMeas_;      //!< number of sampled points
    MaterialType matType_; //!< material parameter to be approximated
    std::string nlFileName_; //!< name of file
  };


} //end of namespace


#endif

