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
#include <def_expl_templ_inst.hh>

#include "General/Environment.hh"
#include "OLAS/solver/BaseSolver.hh"


#include "General/Enum.hh"


// include the original lis header
//extern "C"

#include "lis_config.h"
#include "lis.h"
#include "lis_precon.h"
   //#include <lisf.h>
//}

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

   /** Every call sets up a new preconditionier.
    * @param analysis_id shall be the current info/analysis/progress/step entry and contain an "analysis_id" element */
   void Setup(BaseMatrix &sysmat, PtrParamNode analysis_id);

   /** To satisfy the compiler
    * @param sysmat shall be the one Setup() is called with
    * @param analysis_id @see Setup() */
   void Solve( const BaseMatrix &sysmat,
               const BaseVector &rhs, BaseVector &sol, PtrParamNode analysis_id);

   BaseSolver::SolverType GetSolverType() { return BaseSolver::LIS; }

  private:

   ///Method to read xml definition and create the configuration string
   void createConfigString(PtrParamNode configNode, std::string& output);

   ///Reads information for solver and creates an appropriate string
   void createSolverString(PtrParamNode solverNode, std::string& output);

   ///Reads information for solver and creates an appropriate string
   void createPrecondString(PtrParamNode precondNode, std::string& output);

   //LIS_SOLVER solver_;

   LIS_PRECON precond_;

   LIS_SOLVER solver_;

   //stores the system Matrix
   LIS_MATRIX A_;

   //stores the current solution
   LIS_VECTOR x_;

   ///stores the current RHS
   LIS_VECTOR b_;

   ///pointer to xml node
   PtrParamNode xml_;
  };
}
#endif
