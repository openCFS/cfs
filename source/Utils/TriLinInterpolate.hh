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

    //! Copy Constructor
    TriLinInterpolate(const TriLinInterpolate& right)
      : ApproxData(right){
      this->slices_ = right.slices_;
    }

    virtual TriLinInterpolate* Clone(){
      return new TriLinInterpolate( *this );
    }

    //! destructor
    virtual ~ TriLinInterpolate();

    //computes the approximation polynom
    virtual void CalcApproximation(bool start=true) {;};

    //! computes the regularization parameter
    virtual void CalcBestParameter() {;};

    //! returns y(x)
    virtual Double EvaluateFunc(Double x);

    virtual Double EvaluateFunc(Double x, Double y);

    virtual Double EvaluateFunc(Double x, Double y, Double z);

    //! returns  y'(x)  
    virtual Double EvaluatePrime(Double x) { 
      EXCEPTION(" TriLinInterpolate: EvaluatePrime not implemented");
      return -1.0; 
    };

    ///
    virtual Double EvaluateFuncInv(Double t);

    ///
    virtual Double EvaluatePrimeInv(Double t);
    //  { Error(" TriLinInterpolate:  EvaluatePrimeInv not implemented");};

    ///
    int GetSize() {return numMeas_;};

    ///
    Double EvaluateOrigB(int i) {return y_[i];};

    ///
    Double EvaluateOrigNu(int i) {return x_[i]/y_[i];};


  private:

    StdVector<BiLinInterpolate *> slices_;
    

  };

} //end of namespace


#endif
