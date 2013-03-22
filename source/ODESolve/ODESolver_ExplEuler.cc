#include <fstream>
#include <iostream>
#include <string>

#include <cmath>     


#include "ODESolver_ExplEuler.hh"
#include "General/Exception.hh"


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
          EXCEPTION("Step size too small");
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
 
    EXCEPTION("t>tStop");
    successLastSolve_ = false;
  }
            

} // end of namespace

