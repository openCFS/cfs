// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LININTERPOLATE
#define FILE_LININTERPOLATE

#include<string>
#include "ApproxData.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/vector.hh"

namespace CoupledField {

  //! \todo Replace C Arrays in constructor
  //! by normal vectors and improve error robustness,
  //! if position gets out of range
  class LinInterpolate : public ApproxData
  {
  public:
    //! constructor getting x, y(x)
    LinInterpolate(std::string nlFncName, MaterialType matType );

    //! destructor
    virtual ~ LinInterpolate();

    //computes the approximation polynom
    virtual void CalcApproximation(bool start=true) {;};

    //! computes the regularization parameter
    virtual void CalcBestParameter() {;};

    //! returns y(x)
    virtual double EvaluateFunc(double x);

    //! returns  y'(x)  
    virtual double EvaluatePrime(double x) { 
      EXCEPTION(" LinInterpolate: EvaluatePrime not implemented");
      return -1.0; 
    };

    ///
    virtual double EvaluateFuncInv(double t);

    ///
    virtual double EvaluatePrimeInv(double t);
    //  { Error(" LinInterpolate:  EvaluatePrimeInv not implemented");};

    ///
    int GetSize() {return numMeas_;};

    ///
    double EvaluateOrigB(int i) {return y_[i];};

    ///
    double EvaluateOrigNu(int i) {return x_[i]/y_[i];};


  private:


    int size;
    int node;

    int ind;

    double mu;
    double * mat;
    double * coef;
    double * rhs;
    double * h;
    double * g;


    double alpha;
    double beta;
    double theta;
    double delta;

  };

} //end of namespace


#endif
