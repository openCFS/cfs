#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>

#include "SmoothSpline.hh"

namespace CoupledField
{ 
  SmoothSpline::SmoothSpline( std::string nlFileName,  MaterialType matType )
    : ApproxData( nlFileName, matType )
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
    if ( matType_ == MAG_PERMEABILITY_SCALAR ) {
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
    if ( start  && matType_ == MAG_PERMEABILITY_SCALAR ) {
      CalcStart();
    }

    // compute extrapolation parameter
    double xEndPrime;
    xEndPrime      = (xEnd_ - x_[numMeas_-2]) / (yEnd_ - y_[numMeas_-2]);
    extrapolAlpha_ = (xEnd_ - xEndPrime * yEnd_) / (xEnd_*yEnd_ - nuMax_*yEnd_*yEnd_);
    extrapolBeta_  = ( (xEnd_/yEnd_) - nuMax_ ) * std::exp(extrapolAlpha_*yEnd_);   
  }


  void SmoothSpline::CalcBestParameter()
  {
    Integer i, j;
    bool monotone = true;

    double z,fac;
    double mu_old = 0.0;

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
    double res = 1e14;

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
  }



  void SmoothSpline::ConstructMatrix()
  {
    Integer i,j,k;
    double x1,x2,x3,x4,x5,x6,x7,x8;

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

  void SmoothSpline::ConstructRHS(const Vector<double>& y)
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
    double x1,x2,x3,x4;
  
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
    Vector<double> y(size_+1);
    Vector<double> c(size_*size_+1);
      
    Integer i,j,k;
    double h;

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


  double SmoothSpline::EvaluateFunc(double t) const
  {

    Integer i,j,k;
    double f0,f1,f2,f3,p0,p1,p2,p3;
    double z, val, p;

    if ( t <= xEnd_ ) {
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
    }
    else {
      // linear extrapolation to compute (in case of magnetics) magnetic flux density B
      Double k = (y_[numMeas_-1] - y_[numMeas_-2])/(x_[numMeas_-1] - x_[numMeas_-2]);
      Double d = y_[numMeas_-1] - k*x_[numMeas_-1];
      Double Bval = k*t + d;
      val = t / ( nuMax_ + extrapolBeta_*std::exp(-extrapolAlpha_*Bval) );
    }

    return val;
  }


  double SmoothSpline::EvaluateDeriv(double t) const
  {
    Integer i,j,k;
    double f0,f1,f2,f3,p0,p1,p2,p3;
    double z, val, p;

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


  double SmoothSpline::EvaluateFuncInv( double f ) const
  {
    double z,k,d;
    Integer i;

    if ( f <= y_[numMeas_-1] ) {
      i = Integer( f / theta_ );
      z = Newton(f,g_[i]);
    }
    else {
      // linear extrapolation to compute (in case of magnetics) magnetic flux density B
      k = (y_[numMeas_-1] - y_[numMeas_-2])/(x_[numMeas_-1] - x_[numMeas_-2]);
      d = y_[numMeas_-1] - k*x_[numMeas_-1];
      z = (f - d)/k;
    }

    return z;
  }


  double SmoothSpline::EvaluatePrimeInv( double f ) const
  {
    double z,p;
    Integer i;

    i = Integer ( f / theta_);
    z = Newton( f,g_[i] );
    p = 1.0 / EvaluateDeriv( z );

    return p;
  }


  double SmoothSpline::HermiteFunc( double t, Integer i ) const
  {
    double x = 0.0;

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


  double SmoothSpline::HermitePrime( double t, Integer i ) const
  {

    double x = 0.0;

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


  Integer SmoothSpline::GetInterval( double t ) const
  {
    Integer i;
    double theta;

    if (t < xStart_ ) { 
      //x-Value is smaller than smallest measured H-value!
      return 0;
    }
    if ( t > xEnd_ ) {
      //x-value is larger than karges measured H-value!
      return -1;
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


  double SmoothSpline::Newton( double f, double start ) const
  {
    double za  = start; // start value
    double zn  = xEnd_;
    double eps = 1e-2;

    double rel = za == 0 ? 1 : za;

    while ( fabs((za-zn)/rel) > eps ) {
      zn  = za;
      za -= (EvaluateFunc(za)-f)/EvaluateDeriv(za);
    }

    return za;
  }


  void SmoothSpline::CalcStart()
  {
    double start = 0;

    for (int i=0; i<ind_; i++) {
      g_[i] = Newton( i*theta_, start );
      start = g_[i];
    }

    g_[ind_] = xEnd_;
  }


  bool SmoothSpline::MonotoneBH()
  {
    bool monotone = false;
    double f0,f1,f2,f3;
    double c1,c2,c3;
    double x1,x2,x3;

    for (int i=0; i<=node_; i++) {
      int j  = 2*i;
      
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
    double eps = 1e-6;

    for (i=0; i<=node_; i++) {
      j  = 2*i;

      if ( ( coef_[j+1] - coef_[j+3]) / h_[i] < eps ) {
        return false; // no monotone nu
      }
    }

    return true; // nu curve is monotone
  }


  void SmoothSpline::Print( ) const
  {
    //write result of approximation to ascii file
    std::ofstream out_B;
    out_B.open(std::string(nlFileName_+std::string(".approxBH.dat")).c_str());
    out_B << "#H      B       nu       dnu/dB " << std::endl;

    UInt numPoints = 500;
    double maxH = xEnd_* 2.0; //extend the range of provided max. H-value by factor of 2
    double dH = maxH / ( (double)numPoints );
    double actH, actB;

    // output of the data
    actH = 0.1;
    for ( UInt i=0; i<numPoints; i++) {
      actB = EvaluateFunc(actH) ;
      out_B << actH << "  " << actB << "  " << EvaluateFuncNu(actB) << "  " << EvaluatePrimeNu(actB) << std::endl;
      actH += dH;
    }    
    
    out_B.close(); 
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



}
