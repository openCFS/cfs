#include <iostream>
#include <fstream>
#include <math.h>

#include "Jiles.hh"

namespace CoupledField
{ 

  Jiles::Jiles(Integer numElem, Double ysat, Double a, 
                 Double alpha, Double k, Double c)
    : Hysteresis(numElem)
  {

    Ysaturated_  = ysat;
    a_           = a;
    alpha_       = alpha;
    k_           = k;
    c_           = c;

    Xold_.Resize(numElem);
    Xold_.Init(0);
    YirrOld_.Resize(numElem);
    YirrOld_.Init(0);
    YirrNew_.Resize(numElem);
    YirrNew_.Init(0);

    YirrPrev_.Resize(numElem);
    YirrPrev_.Init(0);

    iterMax_ = 100;
    err_     = 1e-10;
  }

  Jiles::~Jiles()
  {
  }

  Double Jiles::computeValue(Double Xin, Integer idxElem) 
  {

    Integer idx = idxElem - 1;
  
    Double err, dX, YirrOld, YirrNew, Xe, Yan, delta;
    Double a1, b1, c1, nlfunc, dnlfunc, Ynew;
    Integer iter;

    iter = 0;
    err = 1e15;


    //  std::cout << "ComputeValue, Xin=" << Xin << " Xold=" << Xold_[idx] << std::endl;

    //Starting balue
    YirrNew = YirrNew_[idx];

    //YirrNew = YirrOld_[idx];

    dX      = (Xin - Xold_[idx] ) / dt_;

    //  std::cout << "Ystart=" << YirrNew << " dE=" << dX << std::endl;

    iter = 0;
    a1   = alpha_ / dt_;
    do {
      iter++;
      YirrOld = YirrNew;

      Xe = Xin + alpha_*YirrOld;

      //    std::cout << "iter=" << iter << " Xe=" << Xe << std::endl;

      if ( abs(Xe/a_) < 0.1 ) {
        Yan = Ysaturated_ * Xe/(3.0*a_);
      }
      else {
        Yan = Ysaturated_ * ( cosh(Xe/a_) / sinh(Xe/a_) - a_/Xe );
      }

      if (dX >= 0.0 ) {
        delta = 1.0;
      }
      else {
        delta = -1.0;
      }

      b1 = (1.0/dt_) * ( k_*delta - alpha_ * ( Yan + YirrOld_[idx] ) ) + dX;
      c1 = (1.0/dt_) * ( alpha_*YirrOld_[idx]*Yan - YirrOld_[idx]*k_*delta ) - dX*Yan;

      nlfunc  = a1*YirrOld*YirrOld + b1*YirrOld + c1;
      dnlfunc = 2*a1*YirrOld + b1;

      if ( abs(dnlfunc) > 0 ) {
        YirrNew = YirrOld - nlfunc / dnlfunc;
      }
      else {
        YirrNew = 0;
      }

      //    std::cout << "nlfunc=" << nlfunc << " dnlfunc=" << dnlfunc << std::endl;
      //    std::cout << "YirrNew = " << YirrNew << std::endl;

      if ( abs(YirrNew) > 0 ) {
        err = ( YirrNew - YirrOld ) / YirrNew;
      }
      else {
        err = YirrNew - YirrOld;     
      }

    } while ( (iter < iterMax_) && (err > err_) );


    if (iter >= iterMax_ ) {
      EXCEPTION( "No convergence with Jiles hysteresis model in "
               << iterMax_ << " steps" );
    }

    YirrNew_[idx]  = YirrNew;
    YirrPrev_[idx] = YirrNew;

    //  Xold_[idx]     = Xin; 

    Ynew = (1.0 - c_)*YirrNew + c_*Yan;

    //  std::cout << "YirrSol=" << YirrNew << std::endl;

    return Ynew;
  }


  void Jiles::updateMinMaxList(Double Xin, Integer nrEl)
  {

    Integer idx = nrEl-1;


    //YirrOld_[idx] = computeValue(Xin, nrEl);

    YirrOld_[idx] = YirrNew_[idx];
    Xold_[idx]    = Xin; 
  }


}
