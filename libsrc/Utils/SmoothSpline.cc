#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

#include "SmoothSpline.hh"

namespace CoupledField
{ 
SmoothSpline :: SmoothSpline(std::string nlFileName)
  : ApproxData(nlFileName)
{

  delta = 0.01;
  mu    = 1e-14;;
  node  = nummeas-2;       // n = number of measurements
  size  = node*2;    // size of the spline system
  ind  = 200;

  mat  = new double[size*size];
  coef = new double[size+4];
  rhs  = new double[size];
  h    = new double[node+1];
  g    = new double[ind+1];


}

SmoothSpline :: ~SmoothSpline()
{
  delete [] mat;
  delete [] coef;
  delete [] rhs;
  delete [] h;
  delete [] g;

}


void SmoothSpline :: CalcApproximation(int start)
{

  int i;

  alpha = x[0];
  beta  = x[node+1];

  theta = (y[node+1]-y[0])/ind;

  for (i=0; i<=node; i++)
    {
      h[i] = x[i+1]-x[i]; // interval
    }

  // initial conditions

  coef[0]      = y[0];
  coef[1]      = (y[1]-y[0])/h[0];
  coef[size+2] = y[node+1];
  coef[size+3] = (y[node+1]-y[node])/h[node];

  // calculation of the coefficients

  ConstructMatrix();
  ConstructRHS(y);

  CalcCoef();

  // calculation of start values

  if (start == 1)
    {
      CalcStart();
    }
}


void SmoothSpline :: CalcBestParameter()
{
  int i,j,iter = 0,monotone = 1;

  double error = delta;
  double res = 1e14;
  double z,fac;
  double mu_old;

  while (monotone == 1)
    {
      CalcApproximation(0);

      if (MonotoneBH() == -1)
	{
	  mu *= 2;
	}
      else
	{
	  mu_old   = 0.5*mu; 
	  monotone = 0;
	}
    }

  monotone = 1;
  fac      = (mu-mu_old)/10.;
  mu       = mu_old;

  while (monotone == 1)
    {
      CalcApproximation(0);

      if (MonotoneBH() == -1)
	{
	  mu += fac;
	}
      else 
	{
	  monotone = 0;
	}
    }

//  std::cout << "first monotone mu BH " << mu << std::endl;

  monotone = 1;

  while (monotone == 1 && iter <= 1000)
    {
      iter++;

      CalcApproximation(0);

      res = 0;
      j   = 0;
      
      for (i=0; i<=node+1; i++)
	{
	  z = fabs(y[i]-coef[j]);

	  j += 2;

	  if (z > res)
	    {
	      res = z;
	    }
	}

      if ((MonotoneBH() == 1) && (res > 2*error))
	{
	  mu_old = mu;
	  mu *= 2; 
	}
      else if ((MonotoneBH() == -1) && (res > 2*error))
	{
	  mu = mu_old +(mu - mu_old)*0.5;
	}
      else if ((MonotoneBH() == 1) && (res <= 2*error))
	{
	  monotone = 0;
	}
      else if ((MonotoneBH() == -1) && (res <= 2*error))
	{
	  mu = mu_old +(mu - mu_old)*0.5;
	}

    }

//  std::cout << "monotone BH " << MonotoneBH() << std::endl;
//  std::cout << "monotone nu " << MonotoneNu() << std::endl;
}



void SmoothSpline :: ConstructMatrix()
{
  int i,j,k;
  double x1,x2,x3,x4,x5,x6,x7,x8;

  for (i=0; i<size*size; i++)
    {
      mat[i] = 0;
    }

  // construct the matrix
  // D B 0 0
  // A D B 0
  // 0 A D B
  // 0 0 A D
  
  // diagonal entries D

  k = 0;
  
  for (i=0; i<size; i += 2)
    {
      j  = i+1;
      x1 = 1/(h[k]*h[k]*h[k]);
      x2 = 1/(h[k]*h[k]);
      x3 = 1/h[k];
      x4 = h[k];

      x5 = 1/(h[k+1]*h[k+1]*h[k+1]);
      x6 = 1/(h[k+1]*h[k+1]);
      x7 = 1/h[k+1];
      x8 = h[k+1];

      mat[i*size+i]   = 72*x5-144*x6+96*x7 + 72*x1-144*x2+96*x3+2*mu;
      mat[i*size+i+1] = 48*x6-84*x7+48 - 24*x2+60*x3-48;
      mat[j*size+j-1] = 48*x6-84*x7+48 - 24*x2+60*x3-48;
      mat[j*size+j]   = 32*x7-48+24*x8 + 8*x3-24+24*x4;

      k++;
    }

  // super off diagonal B

  k = 1;

  for (i=0; i<size-2; i += 2)
    {
      j  = i+1;
      x1 = 1/(h[k]*h[k]*h[k]);
      x2 = 1/(h[k]*h[k]);
      x3 = 1/h[k];
      x4 = h[k];

      mat[i*size+i+2] = -72*x1+144*x2-96*x3;
      mat[i*size+i+3] = 24*x2-60*x3+48;
      mat[j*size+j+1] = -48*x2+84*x3-48;
      mat[j*size+j+2] = 16*x3-36+24*x4;

      k++;
    }

  // sub off diagonal A

  k = 1;

  for (i=2; i<size; i += 2)
    {
      j  = i+1;
      x1 = 1/(h[k]*h[k]*h[k]);
      x2 = 1/(h[k]*h[k]);
      x3 = 1/h[k];
      x4 = h[k];

      mat[i*size+i-2] = -72*x1+144*x2-96*x3;
      mat[i*size+i-1] = -48*x2+84*x3-48;
      mat[j*size+j-3] = 24*x2-60*x3+48;
      mat[j*size+j-2] = 16*x3-36+24*x4;

      k++;
    }
}

void SmoothSpline :: ConstructRHS(double *y)
{
  int i,j;

  for (i=1; i<size; i += 2)
    {
      rhs[i] = 0;
    }

  j = 1;

  for (i=0; i<size; i += 2)
    {
      rhs[i] = 2*mu*y[j];
      j++;
    }

  // incorporate initial conditions

  double x1,x2,x3,x4;
  
  x1 = 1/(h[0]*h[0]*h[0]);
  x2 = 1/(h[0]*h[0]);
  x3 = 1/h[0];
  x4 = h[0];

  rhs[0] -= (((-72*x1+144*x2-96*x3)*coef[0]+(-48*x2+84*x3-48)*coef[1]));
  rhs[1] -= (((24*x2-60*x3+48)*coef[0]+(16*x3-36+24*x4)*coef[1]));

  x1 = 1/(h[node]*h[node]*h[node]);
  x2 = 1/(h[node]*h[node]);
  x3 = 1/h[node];
  x4 = h[node];

  rhs[size-2] -= (((-72*x1+144*x2-96*x3)*coef[size+2]+(24*x2-60*x3+48)*coef[size+3]));
  rhs[size-1] -= (((-48*x2+84*x3-48)*coef[size+2]+(16*x3-36+24*x4)*coef[size+3]));

}

void SmoothSpline :: CalcCoef()
{
  double * y = new double[size+1];
  double * c = new double[size*size+1];
      
  int i,j,k;
  double h;

  // Solve the system by LU-factorization
      
  // factor

  for (i=1; i<= size; i++)
    {
      for (j=i; j<=size; j++)
	{
	  h = 0;
	  
	  for (k=1; k<=i-1; k++)
	    {
	      h += c[(i-1)*size+k]*c[(k-1)*size+j];
	    }
	  
	  c[(i-1)*size+j] = mat[(i-1)*size+j-1]-h;
	}
      
      for (j=i+1; j<=size; j++)
	{
	  h = 0;
	  
	  for (k=1; k<=i-1; k++)
	    {
	      h += c[(j-1)*size+k]*c[(k-1)*size+i];
	    }
	  
	  c[(j-1)*size+i] = (mat[(j-1)*size+i-1]-h)/c[(i-1)*size+i];
	}
    }

  // solver rs = v, ui = u
  
  y[1] = rhs[0];
  
  for (i=2; i<=size; i++)
    {
      h = 0;
      
      for (j=1; j<=i-1; j++)
	{
	  h += c[(i-1)*size+j]*y[j];
	}
      
      y[i]=(rhs[i-1]-h);
    }
  
  coef[size+1] = y[size]/c[size*size];

  for (i=size-1; i>=1; i--)
    {
      h = 0;
      
      for (j=i+1; j<=size; j++)
	{
	  h += c[(i-1)*size+j]*coef[j+1];
	}

      coef[i+1]=(y[i]-h)/c[(i-1)*size+i];
    }

  // calc the error vector

//   for (i=1; i<=size; i++)
//     {
//       h = 0;

//       for (j=1; j<=size; j++)
// 	{
// 	  h += mat[(i-1)*size+j-1]*coef[j+1];
// 	}

//       y[i] = h;
//     }

//   h = 0;

//   for (i=1; i<=size; i++)
//     {
//       h += (rhs[i-1]-y[i])*(rhs[i-1]-y[i]);
//     }
  
//   std::cout  << "res norm " << sqrt(h) << std::endl;

  delete [] y;
  delete [] c;
}


void SmoothSpline :: CalcMonotoneParameter()
{
  int i,j,iter = 0,monotone = 1;

  double error = delta;
  double res = 1e14;
  double z,fac;
  double mu_old;

  while (monotone == 1)
    {
      CalcApproximation(0);

      if (MonotoneBH() == -1)
	{
	  mu *= 2;
	}
      else
	{
	  mu_old   = 0.5*mu; 
	  monotone = 0;
	}
    }

  monotone = 1;
  fac      = (mu-mu_old)/10.;
  mu       = mu_old;

  while (monotone == 1)
    {
      CalcApproximation(0);

      if (MonotoneBH() == -1)
	{
	  mu += fac;
	}
      else 
	{
	  monotone = 0;
	}
    }

//  std::cout << "first monotone mu BH " << mu << std::endl;

  monotone = 1;

  if (MonotoneNu() == -1)
    {
      std::cout << "+++ ERROR: no approximation can be found" << std::endl;
      return;
    }

  while (monotone == 1 && iter <= 1000)
    {
      iter++;

      CalcApproximation(0);

      res = 0;
      j   = 0;
      
      for (i=0; i<=node+1; i++)
	{
	  z = fabs(y[i]-coef[j]);

	  j += 2;

	  if (z > res)
	    {
	      res = z;
	    }
	}

      if ((MonotoneBH() == 1) && (res > 2*error) && (MonotoneNu() == 1))
	{
	  mu_old = mu;
	  mu *= 2; 
	}
      else if ((MonotoneBH() == -1) && (res > 2*error) && (MonotoneNu() == 1))
	{
	  mu = mu_old +(mu - mu_old)*0.5;
	}
      else if ((MonotoneBH() == 1) && (res <= 2*error) && (MonotoneNu() == 1))
	{
	  monotone = 0;
	}
      else if ((MonotoneBH() == -1) && (res <= 2*error) && (MonotoneNu() == 1))
	{
	  mu = mu_old +(mu - mu_old)*0.5;
	}

      if (MonotoneNu() == -1)
	{
	  mu /= 2;
	  CalcApproximation(0);
	  monotone = 0;
	}

      //std::cout << res << " " << mu << std::endl;
    }

  //  std::cout << "res and mu  " << res << " " << mu << std::endl;
  //  std::cout << "monotone BH " << MonotoneBH() << std::endl;
  //std::cout << "monotone nu " << MonotoneNu() << std::endl;
}


double SmoothSpline :: EvaluateFunc(double t)
{
  int i,j,k;
  double f0,f1,f2,f3,p0,p1,p2,p3;
  double z, val, p;

  i = GetInterval(t);

  if (i == -1)
    {
      std::cout << "no function evaluation possible" << std::endl;
      return 1e+32;
    }

  j = 2*i;

  f0 = coef[j];
  f1 = coef[j+2];
  f2 = coef[j+1];
  f3 = coef[j+3];

  p = alpha;

  for (k=0; k<i; k++)
    {
      p += h[k];
    }

  z = (t-p)/h[i];
	      
  // function value
  
  p0 = HermiteFunc(z,0);
  p1 = HermiteFunc(z,1);
  p2 = HermiteFunc(z,2)*h[i];
  p3 = HermiteFunc(z,3)*h[i];
  
  val= f0*p0+f1*p1+f2*p2+f3*p3;
	      
  return val;
}

double SmoothSpline :: EvaluatePrime(double t)
{
  int i,j,k;
  double f0,f1,f2,f3,p0,p1,p2,p3;
  double z, val, p;

  i = GetInterval(t);

  if (i == -1)
    {
      std::cout << "no prime evaluation possible" << std::endl;
      return 1e+32;
    }

  j  = 2*i;

  f0 = coef[j];
  f1 = coef[j+2];
  f2 = coef[j+1];
  f3 = coef[j+3];

  p = alpha;

  for (k=0; k<i; k++)
    {
      p += h[k];
    }

  z = (t-p)/h[i];

  // prime value
	      
  p0 = HermitePrime(z,0)/h[i];
  p1 = HermitePrime(z,1)/h[i];
  p2 = HermitePrime(z,2);
  p3 = HermitePrime(z,3);
  
  val= f0*p0+f1*p1+f2*p2+f3*p3;
  
  return val;
}

double SmoothSpline :: EvaluateFuncInv(double f)
{
  double z,k,d;
  int i;

  if (f <= y[nummeas-1])
    {
      i = int(f/theta);
      z = Newton(f,g[i]);
    }
  else
    {
      k = (y[nummeas-1] - y[nummeas-2])/(x[nummeas-1] - x[nummeas-2]);
      d = y[nummeas-1] - k*x[nummeas-1];
      z = (f - d)/k;
    }

  return z;
}

double SmoothSpline :: EvaluatePrimeInv(double f)
{
  double z,p;
  int i;

  i = int(f/theta);
  z = Newton(f,g[i]);
  p = 1./EvaluatePrime(z);

  return p;
}

void SmoothSpline :: EvaluateInv(double v, double *f, double *p)
{
  int i;
  
  i  = int(v/theta);
  *f = Newton(v,g[i]);
  *p = 1./EvaluatePrime(*f); 
}




double SmoothSpline :: HermiteFunc(double t, int i)
{
  double x;

  if (i == 0)
    {
      x = (1-t)*(1-t)*(2*t+1);
    }
  else if (i == 1)
    {
      x = t*t*(3-2*t);
    }
  else if (i == 2)
    {
      x = t*(1-t)*(1-t);
    }
  else if (i == 3)
    {
      x = -(1-t)*t*t;
    }

  return x;
}

double SmoothSpline :: HermitePrime(double t, int i)
{
  double x;

  if (i == 0)
    {
      x = -6*t+6*t*t;
    }
  else if (i == 1)
    {
      x = 6*t-6*t*t;
    }
  else if (i == 2)
    {
      x = 1-4*t+3*t*t;
    }
  else if (i == 3)
    {
      x = -2*t+3*t*t;
    }
  
  return x;
}


int SmoothSpline :: GetInterval(double t)
{
  int i;
  double theta;

  if (t < alpha && t > beta)
    {
      return -1;
    }

  theta = alpha;
  i     = 0;

  while (i <= node)
    {
      if (t >= theta && t <= theta+h[i])
	{
	  return i;
	}
      
      theta += h[i];
      i++;
    }

  return i;
}

double SmoothSpline :: Newton(double f,double start)
{
  double za,zn,rel,eps;
  int k;

  za  = start; // start value
  zn  = beta;
  eps = 1e-1;
  k   = 1;

  if (za == 0)
    {
      rel = 1;
    }
  else
    {
      rel = za;
    }

  while (fabs((za-zn)/rel) > eps)
    {
      zn  = za;
      za -= (EvaluateFunc(za)-f)/EvaluatePrime(za);

      k++;
    }

  return za;
}

void SmoothSpline :: CalcStart()
{
  int i;
  double start = 0;

  for (i=0; i<ind; i++)
    {
      g[i] = Newton(i*theta,start);
      start = g[i];
    }

  g[ind] = beta;
}

int SmoothSpline :: MonotoneBH()
{
  int i,j,monotone;
  double f0,f1,f2,f3;
  double c1,c2,c3;
  double x1,x2,x3;

  for (i=0; i<=node; i++)
    {
      j  = 2*i;

      f0 = coef[j]/h[i];
      f1 = coef[j+2]/h[i];
      f2 = coef[j+1];
      f3 = coef[j+3];

      c1 = 6*f0-6*f1+3*f2+3*f3;
      c2 = -6*f0+6*f1-4*f2-2*f3;
      c3 = f2;
      
      if (c2*c2 < 4*c1*c3) // no zero of the quadratic polynomial
	{
	  monotone = 1;
	}
      else if (c2*c2 == 4*c1*c3) // one zero of the quadratic polynomial
	{
	  x1 = -c2/(2*c1);
	  x2 = -c2/(2*c1);

	  if ((2*c1 > 0) && (x1 < 0 || x1 > 1))
	    {
	      monotone = 1;
	    }
	  else
	    {
	      return -1;
	    }
	}
      else // two zeros of the quadratic polynomial
	{
	  x1 = (-c2+sqrt(c2*c2-4*c1*c3))/(2*c1);
	  x2 = (-c2-sqrt(c2*c2-4*c1*c3))/(2*c1);

	  if (x1 > x2)
	    {
	      x3 = x1;
	      x1 = x2;
	      x2 = x3;
	    }

	  if (x1 >= 1 && x2 >= 1 && 2*c1 > 0)
	    {
	      monotone = 1;
	    }
	  else if (x1 <= 0 && x2 <= 0 && 2*c1 > 0)
	    {
	      monotone = 1;
	    }
	  else if (x1 <= 0 && x2 >= 1 && 2*c1 < 0)
	    {
	      monotone = 1;
	    }
	  else
	    {
	      return -1;
	    }
	}
    }

  return monotone; // BH curve is monotone
}

int SmoothSpline :: MonotoneNu()
{
  int i,j;

  for (i=0; i<=node; i++)
    {
      j  = 2*i;

      if ((coef[j+1] - coef[j+3])/h[i] < -1e-6)
	{
	  return -1; // no monotone nu
	}
    }

  return 1; // nu curve is monotone
}



//================================================ just for testing =================================
//

void SmoothSpline :: Read()
{
  int i;

  delta = 0.01;

 std:: ifstream infile;
  
  //infile.open("Vm165_2.fnc");
  //infile.open("nlhfo.fnc");
  infile.open("bhorig.fnc");
  
  if (infile.good() != 1)
    {
      std::cout << "input file not available" << std::endl;
    }
  
  infile >> nummeas;
  
  x = new double[nummeas];
  y = new double[nummeas];
  
  for (i=0; i<nummeas; i++)
    {
      infile >> x[i];
      infile >> y[i];
    }

  infile.close();
}


void SmoothSpline :: Print(double * x, double * y, double lower, double upper)
{
  MakeOutput(x, y, 0, 10000);
  MakeOutputInv(x, y, 0, 2);
}

/////////////////////// private functions /////////////////////////////////////////

void SmoothSpline :: MakeOutput(double * x, double * y, double left, double right)
{
  int i,j;
  double t,z,delta,val;
  double f0,f1,f2,f3,p0,p1,p2,p3;

  std::ofstream out_orig;
  std::ofstream out_func;
  std::ofstream out_prime;

  if (left >= right)
    {
      std::cout << "error" << std::endl;
      return;
    }

  out_orig.open("orig.dat");
  out_func.open("func.dat");
  out_prime.open("prime.dat");

  // output of the data

  for (i=0; i<node+2; i++)
    {
      if (x[i] >= left && x[i] <= right)
	{
	  out_orig << x[i] << " " << y[i] << std::endl;
	}
    }

  // output function 

  j = 0;

  for (i=0; i<node+1; i++)
    {
      f0 = coef[j];
      f1 = coef[j+2];
      f2 = coef[j+1];
      f3 = coef[j+3];

      delta = h[i]/100.;
      t     = x[i];

      while (t <= x[i+1]-delta)
	{
	  if ( t >= left && t <= right)
	    {
	      z = (t-x[i])/h[i];
	      
	      // function value
	      
	      p0 = HermiteFunc(z,0);
	      p1 = HermiteFunc(z,1);
	      p2 = HermiteFunc(z,2)*h[i];
	      p3 = HermiteFunc(z,3)*h[i];
	      
	      val= f0*p0+f1*p1+f2*p2+f3*p3;
	      
	      out_func << t << " " << val << std::endl;
	      
	      // prime value
	      
	      p0 = HermitePrime(z,0)/h[i];
	      p1 = HermitePrime(z,1)/h[i];
	      p2 = HermitePrime(z,2);
	      p3 = HermitePrime(z,3);
	      
	      val= f0*p0+f1*p1+f2*p2+f3*p3;
	      
	      out_prime << t << " " << val << std::endl;
	    }

	  t += delta;
	}
      
      j += 2;
    }

  j -= 2;
  t  = x[node+1];

  if (t >= left && t <= right)
    {
      z  = 1;
      
      f0 = coef[j];
      f1 = coef[j+2];
      f2 = coef[j+1];
      f3 = coef[j+3];
      
      val= y[node+1];
      
      out_func << t << " " << val << std::endl;
      
      p0 = HermitePrime(z,0)/h[node];
      p1 = HermitePrime(z,1)/h[node];
      p2 = HermitePrime(z,2);
      p3 = HermitePrime(z,3);
      
      val = f0*p0+f1*p1+f2*p2+f3*p3;
      
      out_prime << t << " " << val << std::endl;
    }

  out_orig.close();
  out_func.close();
  out_prime.close();
}

void SmoothSpline :: MakeOutputInv(double *x, double *y, double left, double right)
{
  int i;
  double t,delta;

  if (left < y[0])
    {
      left = y[0];
    }

  if (right > y[node+1])
    {
      right = y[node+1];
    }

  delta = (right-left)/100.;
  t     = left;

  std::ofstream out_orig;
  std::ofstream out_func;
  std::ofstream out_prime;

  if (left >= right)
    {
      std::cout << "error" << std::endl;
      return;
    }

  out_orig.open("originv.dat");
  out_func.open("funcinv.dat");
  out_prime.open("primeinv.dat");

  // output of the data

  for (i=0; i<node+2; i++)
    {
      if (y[i] >= left && y[i] <= right)
	{
	  out_orig << y[i] << " " << x[i] << std::endl;
	}
    }

  double f,p;

  while (t <= right+delta)
    {
      EvaluateInv(t,&f,&p);

      out_func << t << " " << f << std::endl;
      out_prime << t << " " << p << std::endl;

      t += delta;
    }

  out_orig.close();
  out_func.close();
  out_prime.close();
}

}
