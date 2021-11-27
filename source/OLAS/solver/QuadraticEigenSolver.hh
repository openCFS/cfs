// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;// ================================================================================================

#ifndef OLAS_QUADRATIC_EIGENSOLVER_HH
#define OLAS_QUADRATIC_EIGENSOLVER_HH

#include "OLAS/solver/BaseEigenSolver.hh"

#include "MatVec/StdMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include <math.h>

#define ABS std::abs


namespace CoupledField {
// =========================================================================
//   Quadratic Eigenvalue Solver
// =========================================================================
//
//! Class for solving quadratic eigenvalue problems
//  Solution is done by linearisation and tranformation to a generalised EVP,
//  which is then solved using a standard eigen solver.
//
class QuadraticEigenSolver : public BaseEigenSolver {

public:

  //! Default Constructor

  QuadraticEigenSolver( shared_ptr<SolStrategy> strat,
      PtrParamNode xml,
      PtrParamNode solverList,
      PtrParamNode precondList,
      PtrParamNode eigenInfo,
      PtrParamNode eSolverList);

  //! Default Destructor
  ~QuadraticEigenSolver();



  //! Setup for a standard EVP
  void Setup(const BaseMatrix & A, bool isHermitian=false)
  {
    EXCEPTION("QuadraticEigenSolver is for quadratic EVPs only");
  }
  //! Setup for a generalised EVP
  void Setup(const BaseMatrix & A, const BaseMatrix & B, bool isHermitian=false)
  {
    //EXCEPTION("QuadraticEigenSolver is for quadratic EVPs only");
    //TODO reinstall the exception once it works
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

  void CalcEigenValues(BaseVector& sol, BaseVector& err, Double minVal, Double maxVal);
  void CalcEigenValues(BaseVector& sol, BaseVector& err, UInt N, Double shiftPoint);
  void CalcEigenValues(BaseVector& sol, BaseVector& err, UInt N, Complex shiftPoint){
    solverForGeneralisedEVP->CalcEigenValues(sol,err,N,shiftPoint);
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

  template<typename T>
  void doLinearisation(T dummy, const BaseMatrix & K, const BaseMatrix & C, const BaseMatrix & M, bool isSymmetric=false);

  UInt GetNnz(const BaseMatrix & M){
    bool isReal, isStoredSymmetric;
    CheckMatrix(isReal, isStoredSymmetric, M);
    UInt nzz;
    if(isStoredSymmetric)
      if(isReal) {
        nzz = dynamic_cast<const SCRS_Matrix<double>&>(M).GetNnz();
      } else {
        nzz = dynamic_cast<const SCRS_Matrix<Complex>&>(M).GetNnz();
      }
    else{
      if(isReal) {
        nzz = dynamic_cast<const CRS_Matrix<double>&>(M).GetNnz();
      } else {
        nzz = dynamic_cast<const CRS_Matrix<Complex>&>(M).GetNnz();
      }
    }
    return nzz;
  };

  // matrices for the linearisation, which are passed to the generalised eigensolver
  BaseMatrix* A_;
  BaseMatrix* B_;

  string  linmode;//mode of linearisation, either first companion form or second companion form

  //! Pointer to eigenvalue solver
  BaseEigenSolver* solverForGeneralisedEVP;
  // Attribute for xml linearisatiin
  PtrParamNode xml_;
  PtrParamNode myParam_;
  PtrParamNode myInfo_;

  void ToInfo();
};
}


#endif // OLAS_QUADRATIC_EIGENSOLVER_HH

