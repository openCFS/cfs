#ifndef FILE_LININTERPOLATE
#define FILE_LININTERPOLATE

#include<string>
#include "Utils/StdVector.hh"
#include "General/environment.hh"
#include "ApproxData.hh"

namespace CoupledField {

  //! \todo Replace C Arrays in constructor
  //! by normal vectors and improve error robustness,
  //! if position gets out of range
  class LinInterpolate : public ApproxData
  {
  public:
    //! constructor getting x, y(x)
    LinInterpolate(std::string nlFncName);

    //! destructor
    virtual ~ LinInterpolate();

    //computes the approximation polynom
    virtual void CalcApproximation(int start=1) {;};

    //! computes the regularization parameter
    virtual void CalcBestParameter() {;};

    ///
    virtual void CalcMonotoneParameter() {;};


    //! returns y(x)
    virtual double EvaluateFunc(double x) {;};

    //! returns  y'(x)  
    virtual double EvaluatePrime(double x)
    { Error(" LinInterpolate: EvaluatePrime not implemented");};

    ///
    virtual double EvaluateFuncInv(double t);

    ///
    virtual double EvaluatePrimeInv(double t);
    //  { Error(" LinInterpolate:  EvaluatePrimeInv not implemented");};

    ///
    int GetSize() {return nummeas;};

    ///
    double EvaluateOrigB(int i) {return y[i];};

    ///
    double EvaluateOrigNu(int i) {return x[i]/y[i];};

    ///
    void Print(double *x, double *y, double lower, double upper);

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
