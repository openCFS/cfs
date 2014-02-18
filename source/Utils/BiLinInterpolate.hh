#ifndef FILE_BILININTERPOLATE
#define FILE_BILININTERPOLATE

#include<string>
#include "Utils/StdVector.hh"
#include "General/Environment.hh"
#include "ApproxData.hh"

namespace CoupledField {

  //! \todo Replace C Arrays in constructor
  //! by normal vectors and improve error robustness,
  //! if position gets out of range
  class BiLinInterpolate : public ApproxData
  {
  public:
    //! constructor getting x, y(x)
    BiLinInterpolate(std::string nlFncName, MaterialType matType );

    //! destructor
    virtual ~ BiLinInterpolate();

    //computes the approximation polynom
    virtual void CalcApproximation(bool start=true) {;};

    //! computes the regularization parameter
    virtual void CalcBestParameter() {;};

    //! returns y(x)
    virtual double EvaluateFunc(double x);

    void findBracketIndices(const double &x, const Vector<Double> & axis, UInt & klo, UInt & khi, double &diff);

    virtual double EvaluateFunc(double x, double y);

    //! returns  y'(x)  
    virtual double EvaluatePrime(double x) { 
      EXCEPTION(" BiLinInterpolate: EvaluatePrime not implemented");
      return -1.0; 
    };

    ///
    virtual double EvaluateFuncInv(double t);

    ///
    virtual double EvaluatePrimeInv(double t);
    //  { Error(" BiLinInterpolate:  EvaluatePrimeInv not implemented");};

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
