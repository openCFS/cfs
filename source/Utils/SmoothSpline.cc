// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

#include "SmoothSpline.hh"

namespace CoupledField
{ 
  SmoothSpline::SmoothSpline( std::string nlFileName,  ApproxCurveType curveType )
    : ApproxData( nlFileName, curveType )
  {

    delta_ = 0.01;
    mu_    = 1e-8;;
    node_  = numMeas_-2;       // nummeas = number of measurements
    size_  = node_*2;          // size of the spline system
    ind_   = 200;

    mat_.Resize(size_*size_);
    coef_.Resize(size_+4);
    rhs_.Resize(size_);
    h_.Resize(node_+1);
    g_.Resize(ind_+1);
 
    g_.Init(0);
    h_.Init(0);
    rhs_.Init(0);
    coef_.Init(0);
    mat_.Init(0);

    nuMax_     = 1.0;
    if ( curveType_ == BH ) {
      //define maximal reluctivity
      nuMax_ = 7.9577e5;
    }
  }

  SmoothSpline::~SmoothSpline()
  {
  }


  void SmoothSpline::CalcApproximation( bool start )
  {

    Integer i;

    xStart_ = x_[0];
    xEnd_  = x_[node_+1];
    yEnd_ = y_[node_+1];
    theta_ = ( yEnd_ - y_[0] ) / ind_; // y-interval

    for ( i=0; i<=node_; i++ ) {
      h_[i] = x_[i+1]-x_[i]; // -interval
    }

    // initial conditions
    coef_[0]       = y_[0];
    coef_[1]       = ( y_[1] - y_[0] ) / h_[0];
    coef_[size_+2] = y_[node_+1];
    coef_[size_+3] = ( y_[node_+1] - y_[node_] ) / h_[node_];

    // calculation of the coefficients
    ConstructMatrix();
    ConstructRHS(y_);

    CalcCoef();

    // calculation of start values
    if ( start  && curveType_ == BH ) {
      CalcStart();
    }

    // compute extrapolation parameter
    Double xEndPrime;
    xEndPrime      = (xEnd_ - x_[numMeas_-2]) / (yEnd_ - y_[numMeas_-2]);
    extrapolAlpha_ = (xEnd_ - xEndPrime * yEnd_) / (xEnd_*yEnd_ - nuMax_*yEnd_*yEnd_);
    extrapolBeta_  = ( (xEnd_/yEnd_) - nuMax_ ) * std::exp(extrapolAlpha_*yEnd_);   
  }


  void SmoothSpline::CalcBestParameter()
  {

    Integer i, j;
    bool monotone = true;

    Double z,fac;
    Double mu_old = 0.0;

    //coarse tuning of discrepancy parameter mu, so that the computed
    //approximation is monoton 
    while ( monotone ) {
      CalcApproximation(0);

      if ( MonotoneBH() == false ) {
        mu_ *= 2;
      }
      else {
        mu_old   = 0.5*mu_; 
        monotone = false;
      }
    }

    //fine tuning of discrepancy parameter mu, so that the computed
    //approximation is monoton  
    monotone = true;
    fac      = (mu_-mu_old)/10.;
    mu_       = mu_old;

    while ( monotone ) {
      CalcApproximation(0);

      if ( MonotoneBH() == false ) {
        mu_ += fac;
      }
      else {
        monotone = false;
      }
    }

    // this is the main loop to get the discrepancy parameter mu;
    Integer  iter = 0;
    Double res = 1e14;

    monotone = true;
    while ( monotone == true && iter <= 1000) {
      iter++;

      CalcApproximation(0);
      
      res = 0;
      j   = 0;
      for ( i=0; i<=node_+1; i++)  {
        z = fabs(y_[i]-coef_[j]);

        j += 2;

        if (z > res) 
          res = z;
      }

      if ( ( MonotoneBH() == true ) && ( res > 2*delta_ ) ) {
        mu_old = mu_;
        mu_ *= 2; 
      }
      else if ( ( MonotoneBH() == false ) && ( res > 2*delta_ ) ) {
        mu_ = mu_old + ( mu_ - mu_old )*0.5;
      }
      else if ( ( MonotoneBH() == true ) && ( res <= 2*delta_ ) ) {
        monotone = false;
      }
      else if ( ( MonotoneBH() == false ) && ( res <= 2*delta_ ) ) {
        mu_ = mu_old + ( mu_ - mu_old ) * 0.5;
      }

    }
    //    std::cout << "muEnd = " << mu_ << std::endl;
  }



  void SmoothSpline::ConstructMatrix()
  {

    Integer i,j,k;
    Double x1,x2,x3,x4,x5,x6,x7,x8;

    for (i=0; i<size_*size_; i++)
      {
        mat_[i] = 0;
      }

    // construct the matrix
    // D B 0 0
    // A D B 0
    // 0 A D B
    // 0 0 A D
  
    // diagonal entries D

    k = 0;
  
    for (i=0; i<size_; i += 2) {
      j  = i+1;
      x1 = 1/(h_[k]*h_[k]*h_[k]);
      x2 = 1/(h_[k]*h_[k]);
      x3 = 1/h_[k];
      x4 = h_[k];
      
      x5 = 1/(h_[k+1]*h_[k+1]*h_[k+1]);
      x6 = 1/(h_[k+1]*h_[k+1]);
      x7 = 1/h_[k+1];
      x8 = h_[k+1];
      
      mat_[i*size_+i]   = 72*x5-144*x6+96*x7 + 72*x1-144*x2+96*x3+2*mu_;
      mat_[i*size_+i+1] = 48*x6-84*x7+48 - 24*x2+60*x3-48;
      mat_[j*size_+j-1] = 48*x6-84*x7+48 - 24*x2+60*x3-48;
      mat_[j*size_+j]   = 32*x7-48+24*x8 + 8*x3-24+24*x4;
      
      k++;
    }

    // super off diagonal B

    k = 1;

    for (i=0; i<size_-2; i += 2) {
      j  = i+1;
      x1 = 1/(h_[k]*h_[k]*h_[k]);
      x2 = 1/(h_[k]*h_[k]);
      x3 = 1/h_[k];
      x4 = h_[k];
      
      mat_[i*size_+i+2] = -72*x1+144*x2-96*x3;
      mat_[i*size_+i+3] = 24*x2-60*x3+48;
      mat_[j*size_+j+1] = -48*x2+84*x3-48;
      mat_[j*size_+j+2] = 16*x3-36+24*x4;
      
      k++;
    }
    
    // sub off diagonal A
    k = 1;
    for (i=2; i<size_; i += 2) {
      j  = i+1;
      x1 = 1/(h_[k]*h_[k]*h_[k]);
      x2 = 1/(h_[k]*h_[k]);
      x3 = 1/h_[k];
      x4 = h_[k];
      
      mat_[i*size_+i-2] = -72*x1+144*x2-96*x3;
      mat_[i*size_+i-1] = -48*x2+84*x3-48;
      mat_[j*size_+j-3] = 24*x2-60*x3+48;
      mat_[j*size_+j-2] = 16*x3-36+24*x4;
      
      k++;
    }
  }

  void SmoothSpline::ConstructRHS(Vector<Double>& y)
  {

    Integer i,j;

    for (i=1; i<size_; i += 2) {
      rhs_[i] = 0;
    }

    j = 1;
    for (i=0; i<size_; i += 2) {
      rhs_[i] = 2*mu_*y[j];
      j++;
    }

    // incorporate initial conditions
    Double x1,x2,x3,x4;
  
    x1 = 1/(h_[0]*h_[0]*h_[0]);
    x2 = 1/(h_[0]*h_[0]);
    x3 = 1/h_[0];
    x4 = h_[0];

    rhs_[0] -= (((-72*x1+144*x2-96*x3)*coef_[0]+(-48*x2+84*x3-48)*coef_[1]));
    rhs_[1] -= (((24*x2-60*x3+48)*coef_[0]+(16*x3-36+24*x4)*coef_[1]));

    x1 = 1/(h_[node_]*h_[node_]*h_[node_]);
    x2 = 1/(h_[node_]*h_[node_]);
    x3 = 1/h_[node_];
    x4 = h_[node_];

    rhs_[size_-2] -= (((-72*x1+144*x2-96*x3)*coef_[size_+2]+(24*x2-60*x3+48)*coef_[size_+3]));
    rhs_[size_-1] -= (((-48*x2+84*x3-48)*coef_[size_+2]+(16*x3-36+24*x4)*coef_[size_+3]));

  }


  void SmoothSpline::CalcCoef()
  {

    Vector<Double> y(size_+1);
    Vector<Double> c(size_*size_+1);
      
    Integer i,j,k;
    Double h;

    // Solve the system by LU-factorization
      
    // factor

    for ( i=1; i<= size_; i++ ) {
      for ( j=i; j<=size_; j++ ) {
        h = 0;
        
        for ( k=1; k<=i-1; k++ ) {
          h += c[(i-1)*size_+k]*c[(k-1)*size_+j];
        }
          
        c[(i-1)*size_+j] = mat_[(i-1)*size_+j-1]-h;
      }
      
      for ( j=i+1; j<=size_; j++ ) {
        h = 0;
          
        for (k=1; k<=i-1; k++) {
          h += c[(j-1)*size_+k]*c[(k-1)*size_+i];
        }
          
        c[(j-1)*size_+i] = (mat_[(j-1)*size_+i-1]-h)/c[(i-1)*size_+i];
      }
    }

    // solver rs = v, ui = u
    y[1] = rhs_[0];
    for ( i=2; i<=size_; i++ ) {
      h = 0;
      
      for ( j=1; j<=i-1; j++ ) {
        h += c[(i-1)*size_+j]*y[j];
      }
      
      y[i]= rhs_[i-1] - h;
    }
  
    coef_[size_+1] = y[size_]/c[size_*size_];

    for ( i=size_-1; i>=1; i-- ) {
        h = 0;
      
        for (j=i+1; j<=size_; j++) {
          h += c[(i-1)*size_+j]*coef_[j+1];
        }

        coef_[i+1]=(y[i]-h)/c[(i-1)*size_+i];
    }
  }


  Double SmoothSpline::EvaluateFunc(Double t)
  {

    Integer i,j,k;
    Double f0,f1,f2,f3,p0,p1,p2,p3;
    Double z, val, p;

    i = GetInterval(t);

    if (i == -1) {
      i = ind_;
    }

    j = 2*i;

    f0 = coef_[j];
    f1 = coef_[j+2];
    f2 = coef_[j+1];
    f3 = coef_[j+3];

    p = xStart_;
    for (k=0; k<i; k++) {
      p += h_[k];
    }

    z = ( t - p ) / h_[i];
              
    // function value
  
    p0 = HermiteFunc(z,0);
    p1 = HermiteFunc(z,1);
    p2 = HermiteFunc(z,2)*h_[i];
    p3 = HermiteFunc(z,3)*h_[i];
  
    val= f0*p0+f1*p1+f2*p2+f3*p3;
              
    return val;
  }


  Double SmoothSpline::EvaluatePrime(Double t)
  {

    Integer i,j,k;
    Double f0,f1,f2,f3,p0,p1,p2,p3;
    Double z, val, p;

    i = GetInterval(t);

    if ( i == -1 ) {
      i = ind_; 
    }

    j  = 2*i;

    f0 = coef_[j];
    f1 = coef_[j+2];
    f2 = coef_[j+1];
    f3 = coef_[j+3];

    p = xStart_;

    for (k=0; k<i; k++) {
      p += h_[k];
    }

    z = ( t - p ) / h_[i];

    // prime value         
    p0 = HermitePrime(z,0)/h_[i];
    p1 = HermitePrime(z,1)/h_[i];
    p2 = HermitePrime(z,2);
    p3 = HermitePrime(z,3);
  
    val= f0*p0+f1*p1+f2*p2+f3*p3;
  
    return val;
  }


  Double SmoothSpline::EvaluateFuncInv( Double f )
  {

    Double z,k,d;
    Integer i;

    if ( f <= y_[numMeas_-1] ) {
      i = Integer( f / theta_ );
      z = Newton(f,g_[i]);
    }
    else {
      k = (y_[numMeas_-1] - y_[numMeas_-2])/(x_[numMeas_-1] - x_[numMeas_-2]);
      d = y_[numMeas_-1] - k*x_[numMeas_-1];
      z = (f - d)/k;
    }

    return z;
  }


  Double SmoothSpline::EvaluatePrimeInv( Double f )
  {

    Double z,p;
    Integer i;

    i = Integer ( f / theta_);
    z = Newton( f,g_[i] );
    p = 1.0 / EvaluatePrime( z );

    return p;
  }


  void SmoothSpline::EvaluateInv( Double v, Double& f, Double& p )
  {

    Integer i;
  
    i = Integer ( v / theta_ );
    f = Newton( v, g_[i] );
    p = 1.0 / EvaluatePrime( f ); 
  }


  Double SmoothSpline::HermiteFunc( Double t, Integer i )
  {

    Double x = 0.0;

    if ( i == 0 ) {
      x = (1-t)*(1-t)*(2*t+1);
    }
    else if ( i == 1 ) {
      x = t*t*(3-2*t);
    }
    else if ( i == 2 ) {
      x = t*(1-t)*(1-t);
    }
    else if ( i == 3 ) {
      x = -(1-t)*t*t;
    }

    return x;
  }


  Double SmoothSpline::HermitePrime( Double t, Integer i )
  {

    Double x = 0.0;

    if ( i == 0 ) {
      x = -6*t+6*t*t;
    }
    else if ( i == 1 ) {
      x = 6*t-6*t*t;
    }
    else if ( i == 2 ) {
      x = 1-4*t+3*t*t;
    }
    else if ( i == 3 ) {
      x = -2*t+3*t*t;
    }
  
    return x;
  }


  Integer SmoothSpline::GetInterval( Double t )
  {

    Integer i;
    Double theta;

    if (t < xStart_ || t > xEnd_) {
      std::cerr << "x-Value is too small -> no convergence!\n";
      //return -1;
      return 0;
    }

    theta = xStart_;
    i     = 0;

    while ( i <= node_ ) {
      if ( t >= theta && t <= theta + h_[i] ) {
        return i;
      }
      
      theta += h_[i];
      i++;
    }

    return i;
  }


  Double SmoothSpline::Newton( Double f, Double start )
  {

    Double za,zn,rel,eps;
    Integer k;

    za  = start; // start value
    zn  = xEnd_;
    eps = 1e-1;
    k   = 1;

    if (za == 0) {
      rel = 1;
    }
    else {
      rel = za;
    }

    while ( fabs((za-zn)/rel) > eps ) {
      zn  = za;
      za -= (EvaluateFunc(za)-f)/EvaluatePrime(za);
      
      k++;
    }

    return za;
  }


  void SmoothSpline::CalcStart()
  {

    Integer i;
    Double start = 0;

    for (i=0; i<ind_; i++) {
      g_[i] = Newton( i*theta_, start );
      start = g_[i];
    }

    g_[ind_] = xEnd_;
  }


  bool SmoothSpline::MonotoneBH()
  {

    Integer i,j;
    bool monotone = false;
    Double f0,f1,f2,f3;
    Double c1,c2,c3;
    Double x1,x2,x3;

    for ( i=0; i<=node_; i++) {
      j  = 2*i;
      
      f0 = coef_[j] / h_[i];
      f1 = coef_[j+2] / h_[i];
      f2 = coef_[j+1];
      f3 = coef_[j+3];
      
      c1 = 6*f0-6*f1+3*f2+3*f3;
      c2 = -6*f0+6*f1-4*f2-2*f3;
      c3 = f2;
      
      if ((c2*c2) < (4*c1*c3) ) {
        // no zero of the quadratic polynomial
        monotone = true;
      }
      else if ((c2*c2) == (4*c1*c3)) {
        // one zero of the quadratic polynomial
        x1 = -c2/(2*c1);
        x2 = -c2/(2*c1);

        if ((2*c1 > 0) && (x1 < 0 || x1 > 1)) {
          monotone = true;
        }
        else {
          return false;
        }
      }
      else {
        // two zeros of the quadratic polynomial
        x1 = (-c2+sqrt(c2*c2-4*c1*c3))/(2*c1);
        x2 = (-c2-sqrt(c2*c2-4*c1*c3))/(2*c1);

        if (x1 > x2) {
          x3 = x1;
          x1 = x2;
          x2 = x3;
        }

        if (x1 >= 1 && x2 >= 1 && 2*c1 > 0) {
          monotone = true;
        }
        else if (x1 <= 0 && x2 <= 0 && 2*c1 > 0) {
          monotone = true;
        }
        else if (x1 <= 0 && x2 >= 1 && 2*c1 < 0) {
          monotone = true;
        }
        else {
          return false;
        }
      }
    }

    return monotone; // BH curve is monotone
  }


  bool SmoothSpline::MonotoneNu()
  {

    Integer i,j;
    Double eps = 1e-6;

    for (i=0; i<=node_; i++) {
      j  = 2*i;

      if ( ( coef_[j+1] - coef_[j+3]) / h_[i] < eps ) {
        return false; // no monotone nu
      }
    }

    return true; // nu curve is monotone
  }



  //================================================ just for testing =================================
  //

  void SmoothSpline::Read()
  {

    UInt i;

    delta_ = 0.01;

    std::ifstream infile;
  
    infile.open("bhorig.fnc");
  
    if (infile.good() != 1){
      std::string str = "Input file for BH-curve with name 'bhorig.dat' not available";
      EXCEPTION( str );
    }
  
    infile >> numMeas_;
  
    x_.Resize(numMeas_);
    y_.Resize(numMeas_);
  
    for (i=0; i<numMeas_; i++) {
      infile >> x_[i];
      infile >> y_[i];
    }

    infile.close();
  }


  void SmoothSpline::Print( )
  {

    MakeOutput( x_, y_ );
    MakeOutputInv( x_, y_ );
    MakeOutputNu();
  }

  /////////////////////// private functions /////////////////////////////////////////

  void SmoothSpline::MakeOutput( Vector<Double>& x, Vector<Double>& y )
  {

    Integer i,j;
    Double t,z,delta,val;
    Double f0,f1,f2,f3,p0,p1,p2,p3;

    std::ofstream out_orig;
    std::ofstream out_func;
    std::ofstream out_prime;

    out_orig.open("orig.dat");
    out_func.open("func.dat");
    out_prime.open("prime.dat");

    // output of the data
    for (i=0; i<node_+2; i++) {
      out_orig << x[i] << " " << y[i] << std::endl;
    }
    

    // output function 
    j = 0;
    for (i=0; i<node_+1; i++) {
      f0 = coef_[j];
      f1 = coef_[j+2];
      f2 = coef_[j+1];
      f3 = coef_[j+3];
      
      delta = h_[i]/100.;
      t      = x_[i];

      while (t <= x_[i+1]-delta) {
        z = (t-x_[i])/h_[i];
          
        // function value
        p0 = HermiteFunc(z,0);
        p1 = HermiteFunc(z,1);
        p2 = HermiteFunc(z,2)*h_[i];
        p3 = HermiteFunc(z,3)*h_[i];
          
        val= f0*p0+f1*p1+f2*p2+f3*p3;
        
        out_func << t << " " << val << std::endl;
        
        // prime value
        p0 = HermitePrime(z,0)/h_[i];
        p1 = HermitePrime(z,1)/h_[i];
        p2 = HermitePrime(z,2);
        p3 = HermitePrime(z,3);
        
        val= f0*p0+f1*p1+f2*p2+f3*p3;
        
        out_prime << t << " " << val << std::endl;
        
        t += delta;
      }
      
      j += 2;
    }
    
    j -= 2;
    t  = xEnd_;

    z  = 1;
    
    f0 = coef_[j];
    f1 = coef_[j+2];
    f2 = coef_[j+1];
    f3 = coef_[j+3];
    
    val= yEnd_;
    
    out_func << t << " " << val << std::endl;
    
    p0 = HermitePrime(z,0)/h_[node_];
    p1 = HermitePrime(z,1)/h_[node_];
    p2 = HermitePrime(z,2);
    p3 = HermitePrime(z,3);
    
    val = f0*p0+f1*p1+f2*p2+f3*p3;
    out_prime << t << " " << val << std::endl;
    
    out_orig.close();
    out_func.close();
    out_prime.close();
  }


  void SmoothSpline::MakeOutputInv( Vector<Double>& x, Vector<Double>& y)
  {

    Integer i;
    Double t,delta;
    
    delta = ( yEnd_ - y_[0] ) / 500.;
    t     = y_[0];
    
    std::ofstream out_orig;
    std::ofstream out_func;
    std::ofstream out_prime;

    out_orig.open("originv.dat");
    out_func.open("funcinv.dat");
    out_prime.open("primeinv.dat");

    // output of the data

    for (i=0; i<node_+2; i++) {
      out_orig << y_[i] << " " << x_[i] << std::endl;
    }
        
    Double f,p;
    while (t <= yEnd_ + delta) {
      EvaluateInv(t,f,p);
      
      out_func << t << " " << f << std::endl;
      out_prime << t << " " << p << std::endl;
      
      t += delta;
    }
    
    out_orig.close();
    out_func.close();
    out_prime.close();
  }



  void SmoothSpline::MakeOutputNu()
  {

    std::ofstream out_nu;
    out_nu.open("nu_B.dat");

    UInt numPoints = 500;
    Double maxB = yEnd_ * 1.5;
    Double dB = maxB / ( (Double)numPoints );

    Double actB = 0;
    // output of the data
    for ( UInt i=0; i<numPoints; i++) {
      out_nu << actB << "  " << EvaluateFuncNu(actB) << "  " << EvaluatePrimeNu(actB) << std::endl;
      actB += dB;
    }

    out_nu.close();
  
  }
}
