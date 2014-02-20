#ifndef FILE_TRILININTERPOLATE
#define FILE_TRILININTERPOLATE

#include<string>
#include "Utils/StdVector.hh"
#include "General/Environment.hh"
#include "ApproxData.hh"
#include "BiLinInterpolate.hh"

namespace CoupledField {

  //! \todo Replace C Arrays in constructor
  //! by normal vectors and improve error robustness,
  //! if position gets out of range
  class TriLinInterpolate : public ApproxData
  {
  public:
    //! constructor getting x, y(x)
    TriLinInterpolate(std::string nlFncName, MaterialType matType );

    //! destructor
    virtual ~ TriLinInterpolate();

    //computes the approximation polynom
    virtual void CalcApproximation(bool start=true) {;};

    //! computes the regularization parameter
    virtual void CalcBestParameter() {;};

    //! returns y(x)
    virtual double EvaluateFunc(double x);

    virtual double EvaluateFunc(double x, double y);

    virtual double EvaluateFunc(double x, double y, double z);

    //! returns  y'(x)  
    virtual double EvaluatePrime(double x) { 
      EXCEPTION(" TriLinInterpolate: EvaluatePrime not implemented");
      return -1.0; 
    };

    ///
    virtual double EvaluateFuncInv(double t);

    ///
    virtual double EvaluatePrimeInv(double t);
    //  { Error(" TriLinInterpolate:  EvaluatePrimeInv not implemented");};

    ///
    int GetSize() {return numMeas_;};

    ///
    double EvaluateOrigB(int i) {return y_[i];};

    ///
    double EvaluateOrigNu(int i) {return x_[i]/y_[i];};


  private:

    StdVector<BiLinInterpolate *> slices_;
    

  };

} //end of namespace


#endif
