// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <cmath>
#include <iostream>

#include "mathfunctions.hh"
#include "MatVec/vector.hh"
#include "Utils/tools.hh"
#include "General/environment.hh"
#include "MatVec/matrix.hh"

namespace CoupledField {

  Double gammaln(Double xx)
  {
    // Internal arithmetic will be done in double precision
    // A nicety that you can omit if five-figure accuracy is good enough.
    Double x,y,tmp,ser;
    static Double cof[6]={76.18009172947146,-86.50532032941677, 24.01409824083091,
                          -1.231739572450155, 0.1208650973866179e-2,-0.5395239384953e-5};
    y=xx;
    x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (Integer j=0; j<=5; j++)
      ser += cof[j]/++y;

    return -tmp+log(2.5066282746310005*ser/x);
  }
 
  /* this method determines an approximation of the eigenvalues lambda of a 
     symmetric positive definite matrix. (lambda * A = lambda x)
     Input: 
     mat: the matrix A whose eigenvalues are wanted. 
     Remark that, after the computation, the result of Givens rotations 
     is overwritten on the array "mat". 
     eps: error tolerance. 
     If the sum of the square of the non-diagonal elements is less than 
     the one of the initial matrix, the iteration is terminated.  
  */
  void eigenValues(Matrix<Double> & mat, Double eps, Vector<Double> & eigen)
  {
    Integer i, j;
    Integer k, kmax;
		// initialize or receive compiler warning
    Integer l_conv = 0;
    Double a2, eps2, dkmax;
    Double n, n2;
    Integer l_sort=1;

    Integer ndim = mat.GetNumRows();
  
    a2 = 0.0;
    for (i=0; i<ndim; ++i)
      for (j=0; j<i-1; ++j)
        a2 += mat[i][j] * mat[i][j];
    a2 *= 2.0;
   
    n = (Double)ndim;
    n2 = n * n;
    eps2 = eps * eps;
    dkmax = std::log(eps)/std::log((n2-n-2)/(n2-n));
    kmax = (Integer)std::ceil(dkmax);

   
    for (k=1; k<kmax; ++k)
      {
        givensRotation(ndim, mat);
        terminationCriterion(ndim, mat, a2, eps2, &l_conv);

        if (l_conv==1)
          break;
      }
    if (l_conv==0) 
      std::cout<<"\n Note: The Jacobi method did not converge!.\n";
  
    for (i=0; i<ndim; ++i)
      eigen[i] = mat[i][i];
    sortArray(ndim, l_sort, eigen);
  }


  void givensRotation(Integer ndim, Matrix<Double> & mat)
  {
    Matrix<Double> b(ndim,ndim);
    
    Double a2, max_a2;
		// initialize or receive compiler warning
    Integer i, j, p = 0, q = 0;
    Double z, t, c, s, u;
  
    max_a2 = 0.0;
    for (i=0; i<ndim; ++i)
      {
        for (j=0; j<=i-1; ++j)
          {
            a2 = mat[i][j];
            a2 = a2 * a2;
            if (a2 > max_a2)
              {
                p = i; 
                q = j;
                max_a2 = a2;
              }
          }
      }
   
    z = 0.5 * (mat[q][q] - mat[p][p]) / mat[p][q];
    t = std::fabs(z) + std::sqrt(1.0 + z*z);
    if (z < 0.0) t = - t;
    t = 1.0 / t;
    c = 1.0 / std::sqrt(1.0 + t*t);
    s = c * t;
    u = s / (1.0 + c);
   
    for (i=0; i<ndim; ++i)
      for (j=0; j<ndim; ++j)
        b[i][j] = mat[i][j];
  
    b[p][p] = mat[p][p] - t * mat[p][q];
    b[q][q] = mat[q][q] + t * mat[p][q];
    b[p][q] = b[q][p] = 0.0;
    for (j=0; j<ndim; ++j)
      if ((j!=p) && (j!=q))
        {
          b[p][j] = b[j][p] = mat[p][j] - s * (mat[q][j] + u*mat[p][j]);
          b[q][j] = b[j][q] = mat[q][j] + s * (mat[p][j] - u*mat[q][j]);
        }
  
    for (i=0; i<ndim; ++i)
      for (j=0; j<ndim; ++j){
        mat[i][j] = b[i][j];


      }
  
  }
  /*
    termination criterion
  
    The iteration is terminated if F(k) <= eps^2 * F(0), 
    where "eps" is the relative error tolerance, 
    F(k) is the sum of all the non-diagonal elements squared of the matrix 
    after the k-th iteration (Given's rotation).

    (input) 
    ndim:    the matrix size
    mat:     double [][NDIM] : the matrix after each Givens rotation. 
    a2:      double : the sum of all the non-diagonal elements squared 
    of the matrix A.
    eps2:    the square of the error tolerance "eps", i.e., eps^2.

    (output) 
    l_conv int * : l_conv = 1 if converged and 
    l_conv = 1 if not yet converged. 
  */
  void terminationCriterion(Integer ndim, Matrix<Double> & mat, Double a2, Double eps2, Integer *l_conv)
  {
    Double a_nd2;
    Integer i, j;
  
    *l_conv = 0;
  
    a_nd2 = 0.0;
    for (i=0; i<ndim; ++i)
      for (j=0; j<i-1; ++j)
        a_nd2 += mat[i][j] * mat[i][j];
    a_nd2 *= 2.0;
  
    if (a_nd2/a2 < eps2) *l_conv = 1;

  }
  /*
    sorting of a 1-dimensional array d(1), d(2), ..., d(n)
    Input: l_sort decides about the order, i.e. l_sort==1 -> sorting downwards
    lsort == 0 -> sorting upwards   */

  void sortArray(Integer ndim, Integer l_sort, Vector<Double> & d)
  {
    Double dv;
    Integer k, i;
 
    if (l_sort == 0)
      {
        for (k=0; k<ndim-1; ++k)
          for (i=k+1; i<ndim; ++i)
            if (d[i] > d[k])
              {
                dv = d[k];
                d[k] = d[i];
                d[i] = dv;
              }
      }
    if (l_sort == 1)
      {
        for (k=0; k<ndim-1; ++k)
          for (i=k+1; i<ndim; ++i)
            if (d[i] < d[k])
              {
                dv = d[k];
                d[k] = d[i];
                d[i] = dv;
              }
      }
  }

  Double Normalize(Vector<Double>& vec)
  {
    Double norm = vec.NormL2();
    if(norm < 1e-20)
        vec *= 0.0;
    else
        vec *= 1.0 / norm;

    return norm;
  }
    
  void CrossProd(const Vector<Double>& a,
                 const Vector<Double>& b,
                 Vector<Double> &result)
  {
#ifdef CHECK_INDEX
    if (a.GetSize() != 3)
      EXCEPTION("Incompatible vector dimensions for CrossProd()");
    if (b.GetSize() != 3)
      EXCEPTION("Incompatible vector dimensions for CrossProd()");
#endif

    result.Resize(3);

    result[0] =  a[1] * b[2] - a[2] * b[1];
    result[1] = -a[0] * b[2] + a[2] * b[0];
    result[2] =  a[0] * b[1] - a[1] * b[0];
  }
    
  bool CoLinear(const Vector<Double> &a, const Vector<Double> &b)
  {
      Vector<Double> c;
      CrossProd(a, b, c);
      return (fabs(c.NormL2()) < 1e-12);
  }

  void GetBarycentricCoords(const Vector<Double>& p1,
                            const Vector<Double>& p2,
                            const Vector<Double>& p3,
                            const Vector<Double>& p,
                            Vector<Double>& b)
  {
    Vector<Double> u, w, v, diff, normal, n;
    Double area_factor;

    w = p2 - p1;
    //    v = p1 - p3;
    u = p3 - p2;

    CrossProd(w, u, normal);
    area_factor = 1.0 / Normalize(normal);

    diff = p - p1;
    CrossProd(w, diff, n);

    normal.Inner(n, b[2]);
    b[2] *= area_factor;

    diff = p - p2;
    CrossProd(u, diff, n);
    normal.Inner(n, b[0]);
    b[0] *= area_factor;

    /*
    diff = p - p3;
    v.Cross(diff, n);
    normal.Inner(n, b[1]);
    b[1] *= area_factor;
    */

    b[1] = 1 - b[0] - b[2];
  }

} // end of namespace
