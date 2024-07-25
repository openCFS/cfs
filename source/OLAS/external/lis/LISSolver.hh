// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     LISSolver.hh
 *       \brief    This implements the interface to the LIS solver package
 *
 *       \date     Jul 28, 2012
 *       \author   ahueppe
 */
//================================================================================================

#ifndef LISSOLVER_HH_
#define LISSOLVER_HH_


#include "General/Environment.hh"
#include "OLAS/solver/BaseSolver.hh"


#include "General/Enum.hh"


#include "lis_config.h"
#include "lis.h"
// note that lis_precon.h is by default not in lis/install/include - we make it manually via cfsdeps
#include "lis_precon.h"

namespace CoupledField
{

  class BaseMatrix;
  class BaseVector;
  class Flags;

  class LISSolver : public BaseIterativeSolver
  {
  public:
    typedef enum {NOSOLVER = 0,
                  CG = 1,
                  BICG = 2,
                  CGS = 3,
                  BICGSTAB = 4,
                  BICGSTABL = 5,
                  GPBICG = 6,
                  TFQMR = 7,
                  ORTHOMIN = 8,
                  GMRES = 9,
                  JACOBI = 10,
                  GS = 11,
                  SOR = 12,
                  BICGSAFE = 13,
                  CR = 14,
                  BICR = 15,
                  CRS = 16,
                  BICRSTAB = 17,
                  GPBICR = 18,
                  BICRSAFE = 19,
                  FGMRES = 20,
                  IDRS = 21,
                  MINRES = 22
            } LISSolverType;
   static Enum<LISSolverType> lisSolverType;

   typedef enum {NONE = 0,
                 JACOBI_PRE = 1,
                 ILU = 2,
                 SSOR = 3,
                 HYBRID = 4,
                 IS = 5,
                 SAINV = 6,
                 SAAMG = 7,
                 ILUC = 8,
                 ILUT = 9,
                 ADDS = 10
           } LISPrecondType;
  static Enum<LISPrecondType> lisPrecondType;

   LISSolver(PtrParamNode param, PtrParamNode olasInfo, BaseMatrix::EntryType type);

   virtual ~LISSolver();

   /** Every call sets up a new preconditionier. */
   void Setup(BaseMatrix &sysmat);

    //! Dummy method: Notify the solver that a new matrix pattern has been set
    void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented for LIS. GetRidOfZeros for NCIs will not work.");};

   /** To satisfy the compiler
    * @param sysmat shall be the one Setup() is called with */
   void Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol);

   BaseSolver::SolverType GetSolverType() { return BaseSolver::LIS; }

  private:

   ///Method to read xml definition and create the configuration string
   void CreateConfigString(PtrParamNode configNode, std::string& output);

   ///Reads information for solver and creates an appropriate string
   void CreateSolverString(PtrParamNode solverNode, std::string& output);

   ///Reads information for solver and creates an appropriate string
   void CreatePrecondString(PtrParamNode precondNode, std::string& output);

   //LIS_SOLVER solver_;

   LIS_PRECON precond_;

   LIS_SOLVER solver_;

   //stores the system Matrix
   LIS_MATRIX A_;
   // this is our working copy we can safely delete without disturbing cfs
   LIS_MATRIX A0_;

   //stores the current solution
   LIS_VECTOR x_;

   ///stores the current RHS
   LIS_VECTOR b_;

   ///pointer to xml node
   PtrParamNode xml_;

   ///internal status flag for updated matrix computations
   bool firstSetup_;

   ///reset solution vector to zero if Setup is called
   bool resetXZero_;

   bool ownMatrixA_;

   /** with throw exception when exceeded */
   int maxIter_ = -1;

   /** not that LIS has three ways to calculate the residuum */
   double tolerance_ = -1.0;

   /** throws only an exception on maxIter exceeded if also minTol_ is not met */
   double minTol_ = -1.0;

   /** activates a lis feature */
   bool logging_ = false;
  };
}
#endif
