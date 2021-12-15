#ifndef FILE_TRILININTERPOLATE
#define FILE_TRILININTERPOLATE

#include <string>
#include "ApproxData.hh"
#include "BiLinInterpolate.hh"
#include "General/Environment.hh"
#include "Utils/StdVector.hh"

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

    //! computes the approximation polynom
    virtual void CalcApproximation(const bool start=true) override {};

    //! returns f(x,y,z)
    virtual double EvaluateFunc(const double x, const double y, const double z) const override;

    //! returns grad f(x,y,z)
    Vector<double> EvaluatePrime(const double x, const double y, const double z) const {
      EXCEPTION("TriLinInterpolate: EvaluatePrime not implemented");
      Vector<double> ret(3,0);
      return ret;
    }

    ///
    int GetSize() const {return numMeas_;};

  private:

    StdVector<BiLinInterpolate *> slices_;
    

  };

} //end of namespace


#endif
