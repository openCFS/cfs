#ifndef FILE_LININTERPOLATE
#define FILE_LININTERPOLATE

#include<string>
#include "Utils/StdVector.hh"
#include "General/Environment.hh"
#include "ApproxData.hh"

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

  };

} //end of namespace


#endif
