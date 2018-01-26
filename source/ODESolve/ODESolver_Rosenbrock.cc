#include <fstream>
#include <iostream>
#include <string>

#include <cmath>     


#include "ODESolver_Rosenbrock.hh"


namespace CoupledField
{
void ODESolver_Rosenbrock::Solve( const Double tInit,
                                  const Double tStop,
                                  StdVector<Double> &yInitOut,
                                  BaseODEProblem &myODE,
                                  Double &hInit,
                                  Double hMin,
                                  Double hMax ) {
  eps_ = 1e-4;

  // TODO: Check if this is still needed
  //  const Double tiny = 1.0e-30;
  UInt i;
  UInt nstp;
  Double t;
  Double h;
  Double hNext;
  Double hDid;

  UInt nvar = yInitOut.GetSize();
  StdVector<Double> yScal(nvar);
  StdVector<Double> C(nvar);
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


  // *** Test ************************************
  if (t+h-tStop > 0.0) {
    h = tStop-t;
  }
  // *********************************************

  for (nstp=0; nstp<maxSteps_; nstp++){
    myODE.CompDeriv(t,y,dydt);
    //Scaling to monitor accuracy.
    //Can be modified if different yscal is needed.
    //Initialise the vector for the error control
    //Choice of C only good for dimensionless case
    for ( i = 0; i < nvar; i++ ){
      C[i] = 1; 
      yScal[i] = fmax(C[i],fabs(y[i]));
    }

    // *** Test ************************************
    // //If stepsize can overshoot, decrease
    //if ((t+h-tStop)*(t+h-tInit) > 0.0){
    //  h = tStop-t;
    //}
    // *********************************************
    RosenbAdaptiveStepsize(y,dydt,t,h,yScal,hDid,hNext,myODE);


    if (hDid == h){
      ++numStepsLastSolve_;
    }
    else{
      ++numBadStepsLastSolve_;
    }	
    //Did we get the final value for tStop?
    // *** MODIFIED testing case, added second condition********
    if (((t-tStop)*(tStop-tInit) >= 0.0) || (t+hNext-tStop > 0.0)){
      for (i=0; i<nvar; i++){
        yInitOut[i] = y[i];
      }
      hInit = hNext; //added compared to nr
      successLastSolve_ = true;
      return;      // Normal exit
    }
    if (fabs(hNext) <= hMin){
      //hNext=hMin;
//       std::cerr << "Step size too small, should be: " << fabs(hNext) << std::endl;
//       std::cerr << "Values of yInit:\n " << yInitOut << std::endl << std::endl;
//       EXCEPTION("Step size too small");
      successLastSolve_ = false;
    }
    h = hNext;
  }

//   EXCEPTION("Too many steps");

  successLastSolve_ = false;


}



void ODESolver_Rosenbrock::RosenbAdaptiveStepsize (StdVector<Double> &y,
    StdVector<Double> &dydt,
    Double &t,
    const Double hTry,
    const StdVector<Double> &yScal,
    Double &hDid,
    Double &hNext,
    BaseODEProblem &myODE){

  // powerGrow and powerShrink are exponents in computation of step size
  Double  powerGrow = -0.25;
  Double  powerShrink = -1.0/3.0;

  Double  Grow   = 1.5;
  Double  Shrink = 0.5;

  // errCon = (Grow/safetyFac_) rasied to the power (1/powergrow)
  Double errCon =  0.012345679;//0.1296;
  UInt maxTry = 80;

  Double GAM, A21, A31, A32, C21, C31, C32, C41, C42, C43, B1, B2, B3, B4;
  Double  E1, E2, E3, E4, A2X, A3X;

  GAM = 1.0 / 2.0;
  A21 = 2.0;
  A31 = 48.0 / 25.0;
  A32 = 6.0 / 25.0;
  C21 = -8.0;
  C31 = 372.0 / 25.0;
  C32 = 12.0 / 5.0;
  C41 = -112.0 / 125.0;
  C42 = -54.0 / 125.0;
  C43 = -2.0/5.0;
  B1  = 19.0 / 9.0;
  B2  = 1.0 / 2.0;
  B3  = 25.0 / 108.0;
  B4  = 125.0 / 108.0;
  E1  = 17.0 / 54.0;
  E2  = 7.0 / 36.0;
  E3  = 0.0;
  E4  = 125.0 / 108.0;
  A2X = 1.0;
  A3X = 3.0 / 5.0;

  UInt i, j, jtry;

  UInt n = y.GetSize();

  StdVector<UInt> indx(n);

  Double d;

  Double  errMax;
  Double  h;
  Double  xsav;

  Matrix<Double> a(n,n);
  Matrix<Double> dfdy(n,n);


  StdVector<Double> dysav(n);     
  StdVector<Double> err(n);
  StdVector<Double> g1(n);
  StdVector<Double> g2(n);
  StdVector<Double> g3(n);
  StdVector<Double> g4(n);
  StdVector<Double> ysav(n);

  xsav = t;
  for(i=0; i<n; i++){
    ysav[i]   = y[i];
    dysav[i] = dydt[i];
  }

  myODE.Jacobi(ysav,dfdy,xsav);

  h = hTry;               //Set stepsize to inital guess

  for( jtry=0; jtry<maxTry; jtry++){
    for(i=0; i<n; i++){
      for(j=0; j<n; j++){
        a(i,j)= -dfdy(i,j);
      }
      a(i,i) += 1.0/(GAM*h);
    }
    LUDecomposition(a,n,indx,d);

    for(i=0; i<n; i++){
      g1[i]= dysav[i];
    }

    LUBackSub(a,n,indx,g1);

    for(i=0; i<n; i++){
      y[i]=ysav[i]+A21*g1[i];
    }
    t=xsav+A2X*h;
    myODE.CompDeriv(t,y,dydt);

    for(i=0; i<n; i++){
      g2[i]= dydt[i]+C21*g1[i]/h;
    }

    LUBackSub(a,n,indx,g2); 
    for(i=0; i<n; i++){
      y[i]=ysav[i]+A31*g1[i]+A32*g2[i];
    }
    t=xsav+A3X*h;
    myODE.CompDeriv(t,y,dydt);
    for(i=0; i<n; i++){
      g3[i]= dydt[i]+(C31*g1[i]+C32*g2[i])/h;
    }
    LUBackSub(a,n,indx,g3);
    for(i=0; i<n; i++){
      g4[i]= dydt[i]+(C41*g1[i]+C42*g2[i]+C43*g3[i])/h;
    }

    LUBackSub(a,n,indx,g4);
    for(i=0; i<n; i++){
      y[i]=ysav[i]+B1*g1[i]+B2*g2[i]+B3*g3[i]+B4*g4[i];
      err[i]=E1*g1[i]+E2*g2[i]+E3*g3[i]+E4*g4[i];
    }

    t= xsav+h;
    if (t == xsav){
      EXCEPTION("Stepsize not significant in stiff");
    }
    errMax = 0.0;
    for(i=0; i<n; i++){
      errMax =( (fabs(err[i]/yScal[i])) > errMax 
          ? (fabs(err[i]/yScal[i])): (errMax));
      //	std::cerr << err[i] << "   " << yScal[i] << "   " <<errMax<<std::endl;
    }
    errMax /= eps_;
    if ( errMax <= 1.0){
      hDid=h;
      hNext=(errMax >errCon? safetyFac_ * h* pow(errMax, powerGrow) : Grow*h);
      //std::cout << "NumSteps: " << jtry << std::endl;
      return;
    }
    else{
      hNext = safetyFac_ *h * pow(errMax,powerShrink);
      h = ( h >= 0.0 ? 
          ((Shrink*h) > hNext ? (Shrink*h) : (hNext)):((Shrink*h) < hNext ? (Shrink*h) : (hNext) ));
    }
  }
  EXCEPTION("Exceeded maxTry in stiff");
}





// LU-Decomposition 
void ODESolver_Rosenbrock::LUDecomposition(Matrix<Double> &a,
    UInt n,
    StdVector<UInt> &indx,
    Double &d){

  const Double tiny = 1.0e-20;
  UInt i, imax, j, k;
  Double big, dum, sum, temp;
  StdVector<Double> v(n);
  d = 1.0;
  imax = 0;

  for (i=0; i<n ; i++){
    big=0.0;
    for (j=0; j<n ; j++){
      if ( (temp=fabs(a(i,j))) > big ) big=temp;
    }
    if ( big == 0){
      EXCEPTION("Singular matrix in routine LUDecomposition");
    }
    v[i]= 1.0 / big;
  }

  for (j=0; j<n ; j++){
    for (i=0; i<j ; i++){
      sum = a(i,j);
      for (k=0; k<i ; k++){
        sum -= a(i,k) * a(k,j);
      }
      a(i,j) = sum;
    }
    big = 0.0;
    for (i=j; i<n ; i++){
      sum = a(i,j);
      for (k=0; k<j ; k++){
        sum -= a(i,k) * a(k,j);
      }
      a(i,j) = sum;
      if ( (dum = v[i] * fabs(sum)) >= big){
        big = dum;
        imax = i;
      }
    }
    if( j != imax){
      for (k=0; k< n ; k++){
        dum = a(imax,k);
        a(imax,k) = a(j,k);
        a(j,k) = dum;
      }
      d = -d;
      v[imax] = v[j];
    }
    indx[j]= imax;
    if ( a(j,j) == 0.0){
      a(j,j)=tiny;
    }
    if (j!=n-1){
      dum = 1.0/(a(j,j));
      for (i=j+1; i<n; i++){
        a(i,j) *= dum;
      }
    }
  }

}

// LU-BAcksubstitution 
void ODESolver_Rosenbrock::LUBackSub(Matrix<Double> &a,
                                     UInt n,
                                     StdVector<UInt> &indx,
                                     StdVector<Double> &b){
  UInt i,ii, ip,j;
  Double sum;
  ii =0;

  for(i=0; i<n; i++){

    ip = indx[i];
    sum = b[ip];
    b[ip]= b[i];
    if( ii!=0) {
      for( j = ii-1; j<i; j++){
        sum -= a(i,j)*b[j];
      }
    }
    else{
      if ( sum!=0.0 ){
        ii = i+1;
      }
    }
    b[i] = sum;
  }

  for(Integer i= n-1; i >= 0; i--){
    sum = b[i];
    for( j= i+1; j< n; j++){
      sum -= a(i,j) *b[j];
    }
    b[i] = sum / a(i,i);
  }

}

} // end of namespace

