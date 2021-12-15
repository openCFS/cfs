#ifndef FILE_BILININTERPOLATE
#define FILE_BILININTERPOLATE

#include <string>
#include "ApproxData.hh"
#include "General/Environment.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {

  //! \todo Replace C Arrays in constructor
  //! by normal vectors and improve error robustness,
  //! if position gets out of range
  class BiLinInterpolate : public ApproxData
  {
  public:
    //! constructor getting x, y(x)
    BiLinInterpolate(std::string nlFncName, MaterialType matType );

    //! Copy Constructor
    BiLinInterpolate(const BiLinInterpolate& right)
      : ApproxData(right){
    }

    virtual BiLinInterpolate* Clone(){
      return new BiLinInterpolate( *this );
    }

    //! destructor
    virtual ~ BiLinInterpolate();

    //computes the approximation polynom
    virtual void CalcApproximation(const bool start=true) override {};

    //! returns f(x,y)
    virtual double EvaluateFunc(const double x, const double y) const override;

    //! returns grad f(x,y)
    Vector<double> EvaluatePrime(const double x, const double y) const {
      EXCEPTION("BiLinInterpolate: EvaluatePrime not implemented");
      Vector<double> ret(2,0);
      return ret;
    }

    ///
    int GetSize() {return numMeas_;};

  private:

  };

} //end of namespace


#endif
