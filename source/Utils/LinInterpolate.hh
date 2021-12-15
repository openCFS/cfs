#ifndef FILE_LININTERPOLATE
#define FILE_LININTERPOLATE

#include <string>
#include "ApproxData.hh"
#include "General/Environment.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {

  //! \todo Replace C Arrays in constructor
  //! by normal vectors and improve error robustness,
  //! if position gets out of range
  class LinInterpolate : public ApproxData
  {
  public:
    //! constructor getting x, f(x)
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
    virtual void CalcApproximation(const bool start=true) override {};

    //! returns f(x)
    virtual double EvaluateFunc(const double x) const override;

    //! returns f'(x)
    virtual double EvaluateDeriv(const double x) const override {
      EXCEPTION("LinInterpolate: EvaluateDeriv not implemented.");
      return 0;
    };

    ///
    double EvaluateFuncInv(const double t) const;

    ///
    double EvaluatePrimeInv(const double t) const;

    ///
    int GetSize() const {return numMeas_;};

    ///
    double EvaluateOrigB(const int i) const {return y_[i];};

    ///
    double EvaluateOrigNu(const int i) const {return x_[i]/y_[i];};


  private:

  };

} //end of namespace


#endif
