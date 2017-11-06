// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     PETSCSolver.hh
 *       \brief    This implements the interface to the PETSC solver package
 *        
 *        The interface consist of  a solver class (ie master)  and a worker class. The sysmat and rhs
 *        obtained is set into PetscMat only using the master class ,all other commands must be called
 *        by all workers
 *              
 *
 *       \date     oct 30, 2017
 *       \author   sri
 */
//================================================================================================

#ifndef PETSCSOLVER_HH_
#define PETSCSOLVER_HH_
#include <def_expl_templ_inst.hh>

#include "General/Environment.hh"
#include "OLAS/solver/BaseSolver.hh"


#include "General/Enum.hh"


// include the original PETSC header
//extern "C"

#include "petsc.h"

 
namespace CoupledField
{

  class BaseMatrix;
  class BaseVector;
  class Flags;

  class PETSCSolver : public BaseIterativeSolver
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
            } PETSCSolverType;
   static Enum<PETSCSolverType> petscSolverType;

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
           } PETSCPrecondType;
  static Enum<PETSCPrecondType> petscPrecondType;

   PETSCSolver(PtrParamNode param, PtrParamNode olasInfo, BaseMatrix::EntryType type);

   virtual ~PETSCSolver();

   /** Every call sets up a new preconditionier. */
   void Setup(BaseMatrix &sysmat);

   /** To satisfy the compiler
    * @param sysmat shall be the one Setup() is called with */
   void Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol);

   void SendWorkerCommand(int TAG);
   
   BaseSolver::SolverType GetSolverType() { return BaseSolver::PETSC; }

  private:

   ///Method to read xml definition and create the configuration string
   //void CreateConfigString(PtrParamNode configNode, std::string& output);

   ///Reads information for solver and creates an appropriate string
   //void CreateSolverString(PtrParamNode solverNode, std::string& output);

   ///Reads information for solver and creates an appropriate string
   //void CreatePrecondString(PtrParamNode precondNode, std::string& output);

   //PETSC_SOLVER solver_;

   PC precond_;//ksp linear precond context

   KSP solver_; //ksp linear solver context

   //stores the system Matrix
   Mat A_;
   // this is our working copy we can safely delete without disturbing cfs
   Mat A0_;

   //stores the current solution
   Vec x_;

   ///stores the current RHS
   Vec b_;

   ///pointer to xml node
   PtrParamNode xml_;

   ///internal status flag for updated matrix computations
   bool firstSetup_;

   ///reset solution vector to zero if Setup is called
   bool resetXZero_;

   bool ownMatrixA_;

   /** with throw exception when exceeded */
   PetscInt maxIter_ = -1;

   /** not that PETSC has three ways to calculate the residuum */
   PetscScalar tolerance_ = -1.0;

   /** throws only an exception on maxIter exceeded if also minTol_ is not met */
   PetscScalar minTol_ = -1.0;

   /** activates a PETSC feature */
   bool logging_ = false;
  
   //To find the rank of the processor its currently in
   PetscMPIInt mpi_rank;

  };
  
  
  class PETSCWorker
  {
  public:
   
    
    void run();
    PETSCWorker();
    ~PETSCWorker();
    
  private:
   
    int InitPetscWorker();
    void SetupPetscWorker(int);
    void GetSol();

    PC precond_;//ksp linear precond context

    KSP solver_; //ksp linear solver context

    //stores the system Matrix
    Mat A_;
    // this is our working copy we can safely delete without disturbing cfs
    Mat A0_;

    //stores the current solution
    Vec x_;

    ///stores the current RHS
    Vec b_;

    ///pointer to xml node
    PtrParamNode xml_;

    ///internal status flag for updated matrix computations
    bool firstSetup_;

    ///reset solution vector to zero if Setup is called
    bool resetXZero_;

    bool ownMatrixA_;

    /** with throw exception when exceeded */
    PetscInt maxIter_ = 10000;

    /** not that PETSC has three ways to calculate the residuum */
    PetscScalar tolerance_ = 1e-12;

    /** throws only an exception on maxIter exceeded if also minTol_ is not met */
    PetscScalar minTol_ =  1e-11;

    /** activates a PETSC feature */
    bool logging_ = false;
  };
  
}
#endif
