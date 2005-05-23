#include <fstream>
#include <iostream>
#include <string>

#include <cmath>     


#include "ODESolver_RKF45.hh"
#include "ODEDescr/Gilmore.hh"
#include "ODEDescr/KellerMiksis.hh"

namespace CoupledField
{
  void ODESolver_RKF45::Solve( const Double tInit,
                               const Double tStop,
                               StdVector<Double> &yInitOut,
                               BaseODEProblem &myODE,
                               Double hInit,
                               Double hMin,
                               Double hMax ){
    const Double tiny = 1.0e-30;
    Integer i;
    Integer nstp;
    Double t;
    Double h;
    Double hNext;
    Double hDid;

    Integer nvar = yInitOut.GetSize();
    StdVector<Double> yScal(nvar);
    StdVector<Double> y(nvar);
    StdVector<Double> dydt(nvar);

    t = tInit;
    //Achtung bei nachfolgender Abfrage Problem durch default wert?
    h = ((tStop-tInit) >= 0 ? 
         (hInit>= 0 ? hInit :-hInit) : (hInit>= 0 ? -hInit : hInit));

    // Initialise status parameters
    numStepsLastSolve_    = 0;
    numBadStepsLastSolve_ = 0;
    successLastSolve_     = 0;

    for ( i = 0; i < nvar; i++ ){
      y[i] = yInitOut[i];
    }
    for (nstp=0; nstp<maxSteps_; nstp++){
      myODE.CompDeriv(t,y,dydt);
      //Scaling to monitor accuracy.
      //Can be modified if different yscal is needed.
      for ( i = 0; i < nvar; i++ ){
        yScal[i] = fabs(y[i])+fabs(dydt[i]*h)+tiny;
      }
      //If stepsize can overshoot, decrease
      if ((t+h-tStop)*(t+h-tInit) > 0.0){
        h = tStop-t;
      }
      RKAdaptiveStepsize(y,dydt,t,h,yScal,hDid,hNext,myODE);
      //  Info-> PrintF("","%e %e  %e  %e\n",t,hDid,y[0],y[1]);
      //(*data)<< t << "    " << y[0] << "    " << y[1] << std::endl;
      if (hDid == h){
        ++numStepsLastSolve_;
      }
      else{
        ++numBadStepsLastSolve_;
      } 
      //Did we get the final value for tStop?
      if ((t-tStop)*(tStop-tInit) >= 0.0){
        for (i=0; i<nvar; i++){
          yInitOut[i] = y[i];
        }
        successLastSolve_ = true;
        //Info-> PrintF("","%d Timestep %e %e %e No. of good steps %d No. of bad steps %d\n",nstp,t,y[0],y[1],numStepsLastSolve_,numBadStepsLastSolve_);
        //      if (numEl_ == 89){
        //  (*data)<< t << "    " << numEl_ <<  "    " << y[0] << "    " 
        //<< y[1];
        //}

        //      if (numEl_ == 90){
        //(*data)<< "    " << numEl_ <<  "    " << y[0] << "    " << y[1]//;
        //               << "   " << yInitOut[0] << "   " 
        //<< yInitOut[1] << std::endl;
        //}

        //if (numEl_ == 91){
        // (*data)<<"    " << numEl_ <<  "    " << y[0] << "    " 
        //<< y[1] << std::endl;
        //}
        return;      // Normal exit
      }
      if (std::abs(hNext) <= hMin){
        hNext=hMin;
        //Error("Step size too small",__FILE__,__LINE__);
        //successLastSolve_ = false;
      }
      h = hNext;
    }
    Double dummyp;
    Double dummydpdt;
    try { 
      //KellerMiksis &theODE = dynamic_cast<KellerMiksis&>(myODE);
      Gilmore &theODE = dynamic_cast<Gilmore&>(myODE);
      dummyp = theODE.GetP();
      dummydpdt = theODE.GetDpdt();
      Info->PrintF("","ElemNo. %d P %e Dpdt %e t %e h %e  R %e dRdt %e\n"
                   ,numEl_,dummyp,dummydpdt,t,h,y[0],y[1]);
      Error( "Too many steps", __FILE__, __LINE__ );
    }
    catch (...) {
      Error( "myODE is not of type Gilmore. Dynamic cast failed!", 
             __FILE__, __LINE__ );
    }
    successLastSolve_ = false;

            
  }



  void ODESolver_RKF45::RKAdaptiveStepsize (StdVector<Double> &y,
                                            StdVector<Double> &dydt,
                                            Double &t,
                                            const Double hTry,
                                            const StdVector<Double> &yScal,
                                            Double &hDid,
                                            Double &hNext,
                                            BaseODEProblem &myODE){
    // powerGrow and powerShrink are exponents in computation of step size
    Double  powerGrow = -0.2;
    Double  powerShrink = -0.25;

    // errCon = (5/safetyFac_) rasied to the power (1/powergrow)
    Double errCon = 1.89e-4;
    //Double errCon =2.48832e-5;
    //Double errCon =1.0e-5;
   
    Integer i;

    Double  errMax;
    Double  h;
    Double  hTemp;
    Double  tNew;
      
    Integer n = y.GetSize();

    h = hTry;               //Set stepsize to inital guess

    StdVector<Double> yError(n);
    StdVector<Double> yTemp(n);

    for(;;){
      RKCashKarp(y, dydt, t, h, yTemp, yError, myODE); //Take one step
      //        std::cerr<<yTemp[0]<<std::endl;
      // Info-> PrintF("","%e  %e\n",yTemp[0],yTemp[1]);
      errMax = 0.0;
      for (i=0; i<n; i++){
        errMax = ( errMax > std::abs(yError[i]/yScal[i]) 
                   ? errMax : (std::abs(yError[i]/yScal[i])));
      }
      errMax /= eps_;             //Scale relative to required tolernace

      if (errMax < 1.0) {
        break;  //Step succeeded. Compute size of next step
      }

      hTemp = safetyFac_ * h * std::pow(errMax,powerShrink);
      //Truncation error too large, reduce stepsize, no more than factor 10
      h = ( h >= 0.0 ?
            (hTemp>0.1*h ? hTemp : 0.1*h):(hTemp<0.1*h ? hTemp : 0.1*h ));
      tNew = t + h;
      if (tNew == t){
        Double dummyp;
        Double dummydpdt;
        Double testr;

        try {
          // KellerMiksis &theODE = dynamic_cast<KellerMiksis&>(myODE);
          Gilmore &theODE = dynamic_cast<Gilmore&>(myODE);
          dummyp = theODE.GetP();
          dummydpdt = theODE.GetDpdt();
          //dummyp=myODE.GetP();
          //    dummydpdt=myODE.GetDpdt();
          Info-> PrintF("","ElemNo. %d P %e Dpdt %e t %e h %e  R %e dRdt %e\n"
                        ,numEl_,dummyp,dummydpdt,t,h,yTemp[0],yTemp[1]);
          Error("Stepsize underflow",__FILE__,__LINE__);
        }
        catch (...) {
          Error( "myODE is not of type Gilmore. Dynamic cast failed!", 
                 __FILE__, __LINE__ );
        }
      }
    }

    //Compute size of next step, no more than factor 5 increase
    if ( errMax > errCon ) {
      hNext = safetyFac_ * h * std::pow(errMax, powerGrow);
    }
    else hNext = 5.0*h;
    hDid = h;
 
    t += hDid;
    for (i=0; i<n; i++) {
      y[i] = yTemp[i];
    }
      
  } 




  void ODESolver_RKF45::RKCashKarp (const StdVector<Double> &y,
                                    const StdVector<Double> &dydt,
                                    const Double &t,
                                    const Double h,
                                    StdVector<Double> &yOut,
                                    StdVector<Double> &yError,
                                    BaseODEProblem &myODE){
    static const Double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875,
      b21 = 0.2, b31 = 3.0/40.0, b32 = 9.0/40.0, b41 = 0.3, b42 = -0.9,
      b43 = 1.2, b51 = -11.0/54.0, b52 = 2.5, b53 = -70.0/27.0,
      b54 = 35.0/27.0, b61 = 1631.0/55296.0,b62 = 175.0/512.0,
      b63 = 575.0/13824.0, b64 = 44275.0/110592.0, b65 = 253.0/4096.0,
      c1 = 37.0/378.0, c3 = 250.0/621.0, c4 = 125.0/594.0, 
      c6 = 512.0/ 1771.0,
      dc1 = c1-2825.0/27648.0, dc3 = c3-18575.0/48384.0,
      dc4 = c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6 = c6-0.25;

    Integer i;

    Integer n = y.GetSize();

    //vectors contain intermediate function values
    StdVector<Double> ak2(n), ak3(n), ak4(n), ak5(n), ak6(n), yTemp(n);

    for ( i = 0; i < n; i++ ) {
      yTemp[i] = y[i] + b21 * h * dydt[i];
    }
    myODE.CompDeriv( t+a2*h, yTemp, ak2);
    for ( i = 0; i < n; i++ ) {
      yTemp[i] = y[i] + h * ( b31 * dydt[i] + b32 * ak2[i] );
    }
    myODE.CompDeriv( t+a3*h, yTemp, ak3);
    for ( i = 0; i < n; i++ ) {
      yTemp[i] = y[i] + h * ( b41 * dydt[i] + b42 * ak2[i] + b43 * ak3[i] );
    }
    myODE.CompDeriv( t+a4*h, yTemp, ak4);
    for ( i = 0; i < n; i++ ) {
      yTemp[i] = y[i] + h * ( b51 * dydt[i] + b52 * ak2[i] + b53 * ak3[i] 
                              + b54 * ak4[i] );
    }
    myODE.CompDeriv( t+a5*h, yTemp, ak5);
    for ( i = 0; i < n ; i++ ) {
      yTemp[i] = y[i] + h * ( b61 * dydt[i] + b62 * ak2[i] + b63 * ak3[i]
                              + b64 * ak4[i] + b65 * ak5[i] );
    }
    myODE.CompDeriv( t+a6*h, yTemp, ak6);
    //Accumulate intermediate values with proper weights
    for ( i = 0; i < n; i++ ) {
      yOut[i] = y[i]+ h * ( c1 * dydt[i] + c3 * ak3[i] + c4 * ak4[i] +
                            c6 * ak6[i] );
    }
    //Estimate error as difference between fourth and fifth order methods
    for ( i = 0; i < n; i++ ) {
      yError[i] = h * ( dc1 * dydt[i] + dc3 * ak3[i] + dc4 * ak4[i]
                        + dc5 * ak5[i] + dc6 * ak6[i] );
    }

  }

  


} // end of namespace

