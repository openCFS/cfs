// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include <cmath>     


#include "ODESolver_ExplEuler.hh"



namespace CoupledField
{
  void ODESolver_ExplEuler::Solve( const Double tInit,
                                   const Double tStop,
                                   StdVector<Double> &yInitOut,
                                   BaseODEProblem &myODE,
                                   Double hInit,
                                   Double hMin,
                                   Double hMax ){

    Integer i;

    Double t;
    Double h;

    Integer nvar = yInitOut.GetSize();
    StdVector<Double> y(nvar);
    StdVector<Double> dydt(nvar);

    t = tInit;
    h = hInit;
    std::cerr <<"t " <<t << " stop " << tStop<<"\n";
    for ( i = 0; i < nvar; i++ ){
      y[i] = yInitOut[i];
    }

    while (t <= tStop){ 
         

      if (t == tStop){
        for (i=0; i<nvar; i++){
          yInitOut[i] = y[i];
          successLastSolve_ = true;
        }
        return;      // Normal exit       
      }
      else{
        if (t + h > tStop){
          h = tStop - t;
        }

        if (std::abs(h) <= hMin){
          Error("Step size too small",__FILE__,__LINE__);
          successLastSolve_ = false;
        }

        myODE.CompDeriv(t,y,dydt);

        for ( i = 0; i < nvar; i++ ){
          y[i] += dydt[i] * h; 
        }

        t += h;

        ++numStepsLastSolve_;
      }
        

    }
 
    Error( "t>tStop", __FILE__, __LINE__ );
    successLastSolve_ = false;
  }
            

} // end of namespace

