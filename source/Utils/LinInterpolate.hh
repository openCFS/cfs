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

    //! Copy Constructor
    LinInterpolate(const LinInterpolate& right)
      : ApproxData(right){
    }

    virtual LinInterpolate* Clone(){
      return new LinInterpolate( *this );
    }
    //! destructor
    virtual ~ LinInterpolate();

    //computes the approximation polynom
    virtual void CalcApproximation(bool start=true) {;};

    //! computes the regularization parameter
    virtual void CalcBestParameter() {;};

    //! returns y(x)
    virtual Double EvaluateFunc(Double x);

    //! returns  y'(x)  
    virtual Double EvaluatePrime(Double x) { 
      EXCEPTION(" LinInterpolate: EvaluatePrime not implemented");
      return -1.0; 
    };

    ///
    virtual Double EvaluateFuncInv(Double t);

    ///
    virtual Double EvaluatePrimeInv(Double t);
    //  { Error(" LinInterpolate:  EvaluatePrimeInv not implemented");};

    ///
    int GetSize() {return numMeas_;};

    ///
    Double EvaluateOrigB(int i) {return y_[i];};

    ///
    Double EvaluateOrigNu(int i) {return x_[i]/y_[i];};


  private:

  };

} //end of namespace


#endif
