// ================================================================================================
/*!
 *       \file     CoefFunction.hh
 *       \brief    Base class for describing coefficients
 *
 *       \date     30/10/2011
 *       \author   Andreas Hueppe
 */
//================================================================================================

#ifndef OLAS_PALM_EIGENSOLVER_HH
#define OLAS_PALM_EIGENSOLVER_HH

#include "OLAS/solver/BaseEigenSolver.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "MatVec/StdMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"

#include "slu_zdefs.h"

#include "MatVec/Vector.hh"

#include <math.h>

#define ABS std::abs

namespace CoupledField {
// =========================================================================
//   PALM SOLVER
// =========================================================================
//
//! Class for interfacing to PALM solvers for Quadratic Eigenvalue problems
//  using Padé approximation linearisation from
//  "A Padé approximate linearization algorithm for solving the quadratic eigenvalue problem with low-rank damping
//   by Ding Lu, Xin Huang, Zhaojun Bai, and Yangfeng Su
//   Int. J. Numer. Methods Eng., 2015. 103(11): 840–858."
//
class PALMEigenSolver : public BaseEigenSolver {

public:

  //! Default Constructor
  PALMEigenSolver( shared_ptr<SolStrategy> strat,
      PtrParamNode xml,
      PtrParamNode solverList,
      PtrParamNode precondList,
      PtrParamNode eigenInfo );

  //! Default Destructor
  ~PALMEigenSolver();



  //! Setup for a standard EVP
  void Setup(const BaseMatrix & A, bool isHermitian=false)
  {
    EXCEPTION("PALM is for quadratic EVPs only");
  }
  //! Setup for a generalised EVP
  void Setup(const BaseMatrix & A, const BaseMatrix & B, bool isHermitian=false)
  {
    EXCEPTION("PALM is for quadratic EVPs only");
  }
  //! setup the quadratic EVP
  void Setup(const BaseMatrix & K, const BaseMatrix & C, const BaseMatrix & M);

  void Setup(const BaseMatrix & mat,  UInt numFreq, double freqShift, bool sort)
  {
    EXCEPTION("old interface - deprecated - do not use");
  }

  void Setup( const BaseMatrix & stiffMat, const BaseMatrix & massMat,
      UInt numFreq, double freqShift, bool sort, bool bloch)
  {
    EXCEPTION("old interface - deprecated - do not use");
  }

  /** Solve the linear generalized eigenvalue problem
   * @see BaseEigenSolver::CalcEigenFrequencies() */
  void CalcEigenFrequencies(BaseVector& sol, BaseVector& err){
    EXCEPTION("obsolete - should be removed from interface");
  }

  void CalcEigenValues(BaseVector& sol, BaseVector& err, Double minVal, Double maxVal){
    EXCEPTION("not implemented yet");
  }
  void CalcEigenValues(BaseVector& sol, BaseVector& err){
    EXCEPTION("not implemented yet");
  }
  void CalcEigenValues(BaseVector& sol, BaseVector& err, UInt N, Double shiftPoint);
  void CalcEigenValues(BaseVector& sol, BaseVector& err, UInt N, Complex shiftPoint){
    EXCEPTION("not implemented yet");
  }

  /**Calculate a particular eigenmode as a postprocessing solution
   * @see BaseEigenSolver::GetEigenMode() */
  void GetEigenMode(unsigned int modeNr, Vector<Complex>& mode, bool right=true);


  void GetComplexEigenMode(unsigned int modeNr, Vector<Complex>& mode)
  {
    EXCEPTION("not implemented yet");
  }


  /** @see BaseEigenSolver::CalcConditionNumber() */
  void CalcConditionNumber( const BaseMatrix& mat, double& condNumber, Vector<double>& evs, Vector<double>& err)
  {
    //EXCEPTION("not implemented yet");
  }

  /** Setup routine for a quadratic eigenvalue problem
   * @see BaseEigenSolver::Setup() */
  void Setup( const BaseMatrix& stiffMat,
      const BaseMatrix& massMat,
      const BaseMatrix& dampMat,
      unsigned int numFreq, double freqShift, bool sort ){
    EXCEPTION("old interface - do not use");
  };

private:

  /* a) Problem variables */

  int   p;       // Pade approximant order.
  int rk;      // Rank of the damping matrix.
  Complex sigma;    // shifts

  const StdMatrix* m_;
  const StdMatrix* c_;
  const StdMatrix* k_;
  StdVector<int> ic_, jc_, im_, jm_, ik_, jk_;
  int    nev;  // Number of eigenvalues to be computed.



  /* b) Internal variables */

  /* b.1) SUPERLU variables */

  double threshold; // Row pivoting threshold, default = 0.1.
  int  order;  // Ordering parameters, default = 1, natural order.
  int  *permc, *permr; // Array of row and column permutation
  SuperMatrix L, U;  // L, U factors
  GlobalLU_t     Glu; // facilitate multiple factorizations with SamePattern_SameRowPerm
  SuperLUStat_t stat;  // Variables.
  bool  isfactor; // Indicates whether LU factorization is done.

  /* b.2) PAL variables */

  Matrix<double> E;
  Matrix<double>   F;     // store the factorization D=EF

  int  n;       // Dimension of the QEP
  int  nl;   // Dimension of the LEP
  double delta;  // scaling factor.
  double  *b;         // E = -diag(b), size p-by-p.
  double  *a;         // length p vector.
  double  normm, normd, normk;   // matrix 1-norms.
  int    nzmax;  // Max number of nonzero elements of Q(sigma)=sigma^2*M+sigma*D+K.
  // (default by nnzm+nnzd+nnzk)
  // (default by eps.)
  bool   isscaling;  // true for applying scaling; false, otherwise.
  bool   isscaled;  // true for applying scaling; false, otherwise.
  bool   iseigv;  // true if eigvectors required; false, otherwise.
  bool   isresid; // true if relative residual norms required; false, otherwise.
  bool   resdone; // true if relative residual norms required; false, otherwise.
  int  nconv;  // number of converged eigenvalues
  // set to -1 at the beginning
  //int  *idxe;  // index of converged eigenvalues and eigenvectors

  /* b.3) Arpack variables */

  double artol;  // (default=eps) Stopping criteria of Arnoldi iteration.

  /* b.4) Data variables */

  double  *res;  // Eigenpair relative residual norms
  Complex* EigVal; // Eigenvalues
  Complex* EigVec; // Eigenvectors

  /* c. 1) Hilfsvariablen */
  UInt* ridx;
  UInt* cidx;

  const SCRS_Matrix<double>* k_CRS; //Stiffness Matrix CSR Format
  const SCRS_Matrix<double>* m_CRS; //Mass Matrix CSR Format
  const SCRS_Matrix<double>* c_CRS; //Damping Matrix CSR Format

  /* Methoden ursprüngliche Implementierung */
  int  np() {return nl;}
  // function that returns the dimension of the problem

  double scaling();
  // Compute matrix norms and scaling

  void  Qsigma(Complex* val, int* colptr, int* rowind,  int& nnz, int& info);
  // Compute Q(sigma) = sigma^2M + sigma*D + K

  void  lufactor();
  // compute lu factorization of Q(sigma)

  void  MultMv(Complex* v, Complex* w);
  // Matrix vector multiplication w <- M*v.
  void  Mult(const SCRS_Matrix<double>* M, Complex *v, Complex *w, Complex alpha);
  void  MultM(Matrix<double>* M, Complex *v, Complex *w, Complex alpha, char trans);
  // Matrix Vector product w=M^T*v

  void  EigComp();
  // compute eigenvalues by ARPACK

  //  void eigvec();
  //  compute eigenvectors
  Complex Eigenvalue( int i);
  // return i-th eigenvalue

  //  void* Eigenvector();
  // return eigenvectors, not support

  void* Eigenvector( int i);
  // return i-th eigenvector

  double Resid( int i);
  // return backward error of i-th eigenpair.

  //void PrintEig( int prec );
  //void PrintEig( );
  // Print eigenvalues and backward error.

  //  void geteig();
  // retrieve eigenvalue and eigenvectors

  void  ResidComp();
  // compute residuals
  void LowRankDcomp();

  double Norm(const SCRS_Matrix<double>* Matr);
  //Computes matrix norm

  int Calculate_DataArrayPosition(UInt, UInt,const UInt*,const UInt*); //Calculates the Position in the CRS_Matrix. If Element is not found returns -1;
  /* Variablen für EigenSolver Class*/
  /** print setup information */
  void ToInfo();

  /** Attribute for xml paramnode of <solver> section */
  PtrParamNode xml_;

  int numFreq_;

  double freqShift_;
  /** set by Setup() */
  bool generalized_;

  /** set by Setup() */
  bool bloch_;

  /* Parameter handler. */
  //void SetTol(double x) { rtol = (x>0)? x: 1.0E-16;};  // backward error tolerance: CAN BE IMPROVED
  //void SetARTol(double x) { artol = (x>0)? x: 1.0E-16;}; // ARPACK stopping criteria
  //void SetARMaxit(int x) { maxit = (x>300||x<0)? 300: x;}; // ARPACK stopping criteria
  void SetLUThresh(double x) { threshold = (x>1||x<0)? 1 : x;};
  void NoScaling() { isscaling = false;}
  void NoEigVec() { iseigv = false;}

  void SetStat(); // output static

};
}


#endif //OLAS_PALM_EIGENSOLVER_HH

