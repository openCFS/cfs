// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ApproxData
#define FILE_ApproxData

#include <string>
#include "General/environment.hh"

namespace CoupledField {

  class ApproxData
  {
  public:
    ApproxData(std::string nlFncName);

    ///
    virtual ~ApproxData();

    ///
    void ReadNlinFunc(std::string fncName);

    ///
    virtual void CalcApproximation(int start=1) = 0;

    ///
    virtual double EvaluateFunc(double t) = 0;

    ///
    virtual double EvaluatePrime(double t) = 0;

    ///
    virtual double EvaluateFuncInv(double t) = 0;

    ///
    virtual double EvaluatePrimeInv(double t) = 0;

    ///
    virtual void CalcBestParameter() = 0;

    ///
    virtual void CalcMonotoneParameter() = 0;

    double EvaluateFuncNu(double t)
    {return (EvaluateFuncInv(t)/t);};

    double EvaluatePrimeNu(double t)
    {return (EvaluatePrimeInv(t)*t - EvaluateFuncInv(t))/(t*t);};

  protected:

    double * x;
    double * y;
    UInt nummeas;

  };


} //end of namespace


#endif

