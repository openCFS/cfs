//#include <stdlib.h>
#include <fstream>
#include <iostream>
//#include <string>
#include <algorithm>

#include <general_head.hh>
#include <utils_head.hh>
#include <datainout_head.hh>
#include <elements_head.hh>
#include <forms_head.hh>
#include <linalg_head.hh>
#include <domain_head.hh>

namespace CoupledField
{

template<class T, class T_Matrix>
LinSystem<T, T_Matrix>::LinSystem(const Double tol)
{
#ifdef TRACE
   (*trace) << "entering LinSystem::LinSystem" << std::endl; 
#endif
 eps=tol;

}

/// This function is only for test methods of this class
template<class Dim, class T_Matrix>
void LinSystem<Dim, T_Matrix>::Set()
{
#ifdef TRACE
   (*trace) << "entering LinSystem::SetAb" << std::endl;
#endif

  Integer n=3;
  A.Resize(n,n);
  b.Resize(n);
 
  Integer i,j;
  for (i=0; i<n; i++)
   for (j=0; j<n; j++)
     A[i][j]=1;
  A[0][0]=3;
  A[1][1]=2;
  A[2][2]=4;
//  for (i=0; i<n; i++) AA[i][i]=1.77778e+8;
//  AA[0][1]=8.88889e+7;
//  AA[1][0]=AA[0][1];
//  A[0][0]=6.0102; A[0][1]=-2.0000;
//  A[1][0]=-2.000; A[1][1]=4.0051;
  cout << A; 
  for (i=0; i<n; i++) b[i]=100;
  cout << b;
}

template<class T, class T_Matrix>
LinSystem<T, T_Matrix>::~LinSystem()
{
#ifdef TRACE
   (*trace) << "entering LinSystem::~LinSytem" << std::endl;
#endif
  ;
}

template<class T, class T_Matrix>
Boolean LinSystem<T, T_Matrix>::CG(const Integer maxIter,  enum precond typePrecond) 
{
//***************************************************************************
// x - on entry: initial guess
//   - on return: approximate solution
// b - right hand side
// maxIter - max number of iteration
// typePrecond - type of precondition {}
//***************************************************************************

#ifdef TRACE
   (*trace) << "entering LinSystem::CG" << std::endl;
#endif
   Clock oClock;
   oClock.ClockCount(Clock::beg);

  Integer iter;
  T alpha, beta, pAp, zr, zrOld;
  Double stop;
  Integer dim=b.size();
  x.Resize(dim);


  if (x.size()!=b.size())
       Error("Size of vector for solution is not equal size of system matrix",
              __FILE__,__LINE__);

  Vector<T> p(dim), r(dim), Ap(dim), z(dim);

  r=b-A*x;
  A.precond(z,r,typePrecond);
  p=z;
  zr=z*r;
  stop=eps*b.norm_2();  

  for (iter=0; iter < maxIter; iter++)
{ 
  Ap=A*p;
  pAp=Ap*p;
  alpha=zr/pAp;
 
  x+=p*alpha;
  r-=Ap*alpha;

  if ( r.norm_2()<=stop) break;

  A.precond(z,r,typePrecond);

  zrOld=zr;
  zr=z*r;
  
  beta=zr/zrOld;
  p=z+p*beta;   
}

if (InfoPrint)
 (*infofile) << "--------------  Conjugate Gradiate Method --------------" <<
  std::endl << "Number of max iteration is " << maxIter << endl 
       << "Number of iterations is " << iter << std::endl
       << "Precondition is " << typePrecond << std::endl
  << " ----------------------------------------------------- " << std::endl;

  oClock.ClockCount(Clock::end, " CG ");
  if (iter==maxIter) return FALSE;
  else return TRUE;
}

template<class T, class T_Matrix>
Boolean LinSystem<T, T_Matrix>::GMRes_m(const Integer maxIter, enum precond typePrecond, const Integer m)
{
//***************************************************************************
// x - on entry: initial guess
//   - on return: approximate solution 
// b - right hand side
// maxIter - max number of iteration
// typePrecond - type of precondition {}
// m - number of iteration for restarting GMRes_m
//***************************************************************************

#ifdef TRACE
  (*trace) << "entering LinSystem::GMRes_m" << std::endl;
#endif
   Double stop, normaW;

   Integer n=b.size();
  
   x.Resize(n); ///////////////777777777777777777777777777
   Integer iter,i,j,k;

   Vector<T> help(n),r(n),w(n), cs(m), sn(m);
   Vector<T> g(m+1); // 

   // orthonormal basis for Krylov space
   // v[i] - pointer to Vector(n)
   // We delete it by default with deconstructor

    Vector<T> * v=new Vector<T>[m+1];
    for (i=0; i<=m; i++) v[i].Allocate(n);

   // Hessenburg matrix h, (m+1) by m
   // h[i] stores ith column of h that has i+2 entries
   T** h=new T*[m];         
   for (i=0; i<m; i++) h[i]=new T[i+2];

   Boolean conv=FALSE;
   T tmp;
   
   // Define stopping criterion

   A.precond(help,b,typePrecond);
   stop=help.norm_2()*eps;

   // Calculate first residual
   
   A.precond(r,b-A*x,typePrecond);
   Double beta=r.norm_2();

   /// Check: Is our initial guess suited as approximate solution
   if (beta<=stop) { cout << " 777 You guessed solution "; return 0;}

   iter=1;
   while (iter<=maxIter)
{
    cout << r << " r " << std::endl;
    cout << beta << " beta " << std::endl;
    v[0]=r*(1.0/beta);
    g[0]=beta;

    for (k=0; k<m && iter <=maxIter; k++, iter++)
    {
      cout << std::endl <<  "NEW ITERATION " << k << endl;
      A.precond(w,A*v[k],typePrecond);
      cout << w << "w" << std::endl;  /// VVVV
      cout << v[0] << "v[0]" << std::endl;
      normaW=w.norm_2();
      for (i=0; i<=k; i++) { h[k][i]=w*v[i]; 
                             cout << h[k][i] << " h[k][i]" << std::endl;
                               //// VVVV
                             cout << v[i]*h[k][i] << std::endl;
                             w-=v[i]*h[k][i];  //// VVVV
                             cout << "\t" << i << " i"  << w << " w " << std::endl;
                            }
      cout << w << " w " << std::endl;
      cout << w.norm_2() << " norma w" << std::endl;
      h[k][k+1]=w.norm_2();
      cout << v[0] << std::endl;
      cout << h[k][k+1] << " HKK+1 " << i << " i " << k << " k " <<  std::endl;
     
   /// If h[k][k+1] is small, do reorthogonalization

      if (normaW+1.0e-4*h[k][k+1] == normaW)
   { for (i=0; i<=k; i++) 
         {
           tmp=w*v[i];
           h[k][i]+=tmp;
           w-=v[i]*tmp;
         }
       h[k][k+1]=w.norm_2();
    }
     
    cout << h[k][k+1] << " after reorthogonalization " << std::endl; 
    if (h[k][k+1] == 0) Error("zero-divisor in GMRes_m()");
    cout << " w " << w << " w " << std::endl;
    v[k+1]=w/h[k][k+1];   //// VVVV
    cout << " v[k+1] " << v[k+1] << std::endl;
    /// Apply rotations G0,G1, ... ,Gk-1 to column k of h

    cout << __LINE__ << std::endl;

    for (i=0; i<k; i++)
{   
    tmp=h[k][i]*cs[i]+h[k][i+1]*sn[i];
    h[k][i+1]=-sn[i]*h[k][i]+cs[i]*h[k][i+1];
    h[k][i]=tmp;
} 

      cout << __LINE__ << std::endl; 
    /// Generate rotations

    if (h[k][k+1] == 0) { cs[k]=1; sn[k]=0;}
       else if (abs(h[k][k+1])>abs(h[k][k])) 
        {  tmp=h[k][k]/h[k][k+1];
            sn[k]=1.0/sqrt(1.0+tmp*tmp);
            cs[k]= tmp*sn[k];}
       else 
        {  tmp=h[k][k+1]/h[k][k];
            sn[k]=1.0/sqrt(1.0+tmp*tmp);
            cs[k]= tmp*sn[k];}
  
    /// Apply rotation Gk to column k of h and to g       

      cout << __LINE__ << std::endl; 
    tmp=h[k][k]*cs[k]+h[k][k+1]*sn[k];
    h[k][k+1]=-sn[k]*h[k][k]+cs[k]*h[k][k+1];
    h[k][k]=tmp;

      cout << __LINE__ << std::endl; 
    tmp=g[k]*cs[i]+g[k+1]*sn[i];
    g[k+1]=-sn[i]*g[k]+cs[i]*g[k+1];
    g[k]=tmp;

(*trace) << g[k+1] << " g[k+1] " << std::endl << stop << " stop " << endl;

    if (abs(g[k+1])<=stop) { k++; cout << " break " << std::endl; break;}
      cout << __LINE__ << std::endl; 
   } /// end of loop

   // back solve the upper triangular system
     cout << __LINE__ << std::endl; 
   for (i=k-1; i>=0; i--) {
     for (j=i+1; j<k; j++) g[i]-=h[j][i]*g[j];
                     g[i]/=h[i][i];
                           }
      cout << "g[i]" << g << std::endl;
      cout << __LINE__ << std::endl;  
    for (i=0; i <k; i++) x+=v[i]*g[i];  /// VVVV

    A.precond(r, b-A*x, typePrecond);
    cout << r << " r " << std::endl;
    beta=r.norm_2();
    cout << " beta " << beta << std::endl;
    if (abs(beta)<=stop) { conv=TRUE; break;}

}
   Double tol=(b-A*x).norm_2();
  
//   for (i=0; i<=m; i++) delete v[i];   //// VVVVVV
   delete [] v;
   
   
   for (i=0; i<m; i++) delete [] h[i];
   delete [] h;

   cout << " solution " << x << std::endl;

#ifdef TRACE
 (*trace) << " leaving LinSystem::GMRes_m " << std::endl;
#endif

   return conv;

}

template<class T, class T_Matrix>
Boolean LinSystem<T, T_Matrix>::BiCGSTAB(const Integer maxIter, enum precond typePrecond)
{
//***************************************************************************
// x - on entry: initial guess
//   - on return: approximate solution
// b - right hand side
// maxIter - max number of iteration
// typePrecond - type of precondition {}
//***************************************************************************
 
#ifdef TRACE
  (*trace) << "entering LinSystem::BiCGSTAB" << std::endl;
#endif

 if (InfoPrint)
  (*infofile) <<   "--------------  BiCGSTAB  --------------" <<
  std::endl << "Number of max iteration is " << maxIter << endl
       << "Precondition is " << typePrecond << std::endl
  << " ----------------------------------------------------- " << std::endl;

 Double tol;
 Double rho_1, rho_2=1, alpha=1, beta, omega=1;
 Vector<Double> p, phat, s, shat, t, v;

 Integer numiter=0;
 Integer n=b.size();
 phat.Allocate(n);
 shat.Allocate(n); 

 x.Resize(n); 
 Vector<Double> r=b-A*x;
 Vector<Double> rtilde=r;

 Double normb = b.norm_2();
 if (normb == 0) normb=1;
 
 Double stop=eps*normb;
 if (r.norm_2() <= stop ) { tol=r.norm_2()/normb; return TRUE; }

 Integer i;
 for (i=0; i<=maxIter; i++)
  { rho_1=rtilde*r;

    if (i==0) p=r;
    else { beta=(rho_1/rho_2)*(alpha/omega);
           p=r+(p-v*omega)*beta;
         }
     
    A.precond(phat,p,typePrecond); 
    v=A*phat;
    alpha=rho_1/(rtilde*v);
    s=r-v*alpha;
    if (s.norm_2()<= stop) { x+=phat*alpha; tol=r.norm_2()/normb; return FALSE;} 
    A.precond(shat,s,typePrecond);
    t=A*shat;
    omega=(t*s)/(t*t);
    x+=phat*alpha+shat*omega;
    r=s-t*omega;
   
    rho_2=rho_1;
    if (r.norm_2()<= stop) { tol=r.norm_2()/normb; return TRUE;} 
    if ( omega == 0) { tol=r.norm_2()/b.norm_2();
                       return FALSE;}
}
   tol=r.norm_2()/b.norm_2(); 
   return TRUE;
} 
} // end of namespace
