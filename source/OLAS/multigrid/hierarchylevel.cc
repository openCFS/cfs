// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */
#include <def_use_lapack.hh>
#include <def_use_pardiso.hh>

#include "OLAS/multigrid/hierarchylevel.hh"
#include "OLAS/multigrid/gaussseidel.hh"
#include "OLAS/multigrid/AFWsmoother.hh"
#include "OLAS/multigrid/jacobi.hh"
#include "OLAS/multigrid/transfer.hh"

#include "OLAS/solver/generatesolver.hh"
#include "OLAS/solver/BaseSolver.hh"
#include "OLAS/precond/IdPrecondStd.hh"

#ifdef USE_LAPACK
#include "OLAS/external/lapack/lapack.hh"
#endif


#ifndef  LOGLINE
#define  LOGLINE  " --------------------------------------\n"
#endif
/**********************************************************/

namespace CoupledField {
/**********************************************************/

template <typename T>
HierarchyLevel<T>::HierarchyLevel()
: LevelID_( 0 ),
  SysMatrix_( NULL ),
  AuxMatrix_( NULL ),
  DirSysMatrix_( NULL ),
  dirSysMatPPool_( NULL ),
  Topology_( NULL ),
  Transfer_( NULL ),
  Agglomerate_(NULL),
  Smoother_( NULL ),
  directSolver_( NULL ),
  directSolverParams_( NULL ),
  vh1_( NULL ),
  vH1_( NULL ),
  vH2_( NULL ),
  nextLevel_( NULL ),
  logging_( true ),
  deleteDirSysMatrix_( true ),
  amgType_(AMGType::SCALAR)
  {
  }

/**********************************************************/

template <typename T>
HierarchyLevel<T>::HierarchyLevel( const Integer level_id, const AMGType amgType )
: LevelID_( level_id ),
  SysMatrix_( NULL ),
  AuxMatrix_( NULL ),
  DirSysMatrix_( NULL ),
  dirSysMatPPool_( NULL ),
  Topology_( NULL ),
  Transfer_( NULL ),
  Agglomerate_(NULL),
  Smoother_( NULL ),
  directSolver_( NULL ),
  directSolverParams_( NULL ),
  vh1_( NULL ),
  vH1_( NULL ),
  vH2_( NULL ),
  nextLevel_( NULL ),
  logging_( true ),
  deleteDirSysMatrix_( true ),
  amgType_(amgType)
  {
  }

/**********************************************************/

template <typename T>
HierarchyLevel<T>::~HierarchyLevel()
{
  Reset();
}

/**********************************************************/

template <typename T>
inline const HierarchyLevel<T>* HierarchyLevel<T>::
GetNextLevel() const
{
  return nextLevel_;
}

/**********************************************************/

template <typename T>
inline Integer HierarchyLevel<T>::GetLevelID() const
{
  return LevelID_;
}

/**********************************************************/

template <typename T>
void HierarchyLevel<T>::Reset()
{

  if( !GetLevelID() && (SysMatrix_ || AuxMatrix_) ) {
    std::cerr
    << "called \033[31mHierarchyLevel::Reset\033[0m"
    << std::endl << "on level "<<GetLevelID()<<" with still"
    " present system or auxiliary matrix" << std::endl
    << "might not be intended -> eventually crash"
    << std::endl;
  }

  delete SysMatrix_;     SysMatrix_ = NULL;
  delete AuxMatrix_;     AuxMatrix_ = NULL;

  if( deleteDirSysMatrix_ ) {
    delete DirSysMatrix_;
    DirSysMatrix_ = NULL;
  }
  deleteDirSysMatrix_ = true;
  // It is crucial to delete the pattern pool AFTER
  // DirSysMatrix_, since the destructor of DirSysMatrix_
  // deregisters itself at the pool. Otherwise the pool
  // would not free the pattern memory.
  delete dirSysMatPPool_;  dirSysMatPPool_ = NULL;

  delete Topology_;      Topology_     = NULL;
  delete Transfer_;      Transfer_     = NULL;
  delete Agglomerate_;   Agglomerate_  = NULL;
  delete Smoother_;      Smoother_     = NULL;
  delete directSolver_;  directSolver_ = NULL;
  delete directSolverParams_;  directSolverParams_ = NULL;

  // delete temporary vectors
  delete vh1_;  vh1_ = NULL;
  delete vH1_;  vH1_ = NULL;
  delete vH2_;  vH2_ = NULL;

  // remove all coarser levels
  delete nextLevel_;  nextLevel_ = NULL;
}

/**********************************************************/

template <typename T>
void HierarchyLevel<T>::
InsertSysMatrix( const CRS_Matrix<T>* const sysmat )
{

  delete SysMatrix_;
  SysMatrix_ = (CRS_Matrix<T>*)sysmat;
}


template <typename T>
void HierarchyLevel<T>::
InsertAuxMatrix( const CRS_Matrix<Double>* const auxmat )
{

  delete AuxMatrix_;
  AuxMatrix_ = (CRS_Matrix<Double>*)auxmat;
}

/**********************************************************/

template <typename T>
const CRS_Matrix<T>* HierarchyLevel<T>::GetSysMatrixPtr() const
{

  return (const CRS_Matrix<T>*) SysMatrix_;
}


template <typename T>
const CRS_Matrix<Double>* HierarchyLevel<T>::GetAuxMatrixPtr() const
{

  return (const CRS_Matrix<Double>*) AuxMatrix_;
}

/**********************************************************/

template <typename T>
CRS_Matrix<T>* HierarchyLevel<T>::UnhookSysMatrix()
{

  CRS_Matrix<T> *temp = SysMatrix_;
  SysMatrix_ = NULL;
  return temp;
}


template <typename T>
CRS_Matrix<Double>* HierarchyLevel<T>::UnhookAuxMatrix()
{

  CRS_Matrix<Double> *temp = AuxMatrix_;
  AuxMatrix_ = NULL;
  return temp;
}

/**********************************************************/

template <typename T>
bool HierarchyLevel<T>::Setup( Settings* const settings,
                              Integer numBadCoarsenings){


  if( SysMatrix_ == NULL )  return false;
  if( AuxMatrix_ == NULL )  return false;
  SetLogging( settings->logging );
  /*********************************************************/
  /**************** SETUP of THIS level ********************/
  /*********************************************************/
  /* Consists of the following steps
   * 1) Setup the topology
   * 2) Calculate Coarse Fine Splitting
   * 3) Setup auxliary data
   * 4) Setup a direct solver (PARDISO)
   * 5) Setup a smoother (GS or Jacobi)
   * 6) Create operators (P, R)
   * 7) Calculate Galerkin Product
   *
   * Perform recursively on ever deeper levels until
   * the coarsest one is reached
   */


  /**************** Step 1) ********************************/
  // setup the topology (also creates the dependency graphs)
  if( ! SetupTopology(settings->alpha,
      settings->strongDiagRatio,
      settings->forceFineRatio) ) {
    return false;
  }

  /**************** Step 2) ********************************/
  // compute the coarse-fine splitting.
  const Integer SizeH =
      Topology_->CalcCoarseFineSplitting( settings->maxCoarseDepend,
          settings->gamma );


  /**************** Step 3) ********************************/
  // create auxiliary data. Here we need the sizes of both, the
  // fine and the coarse system (this is the reason, why we stored it)
  // In case SizeH > settings->minSystemSize, this level will stay
  // the bottom level, and we do not need the temporary vectors of
  // size nH. We prevent their creation by passing 0 as second
  // parameter to SetupAuxData.
  UInt dim = SysMatrix_->GetNumRows() / AuxMatrix_->GetNumRows();

  if( ! SetupAuxData(Topology_->GetSizeh(),
      SizeH < settings->minSystemSize ? 0 : SizeH,
          dim) ) {
    return false;
  }

  // check, if the coarsening gets slow
  if( settings->badCoarseningRatio * Topology_->GetSizeh() <= SizeH ) {
    numBadCoarsenings++;
  } else {
    numBadCoarsenings = 0;
  }

  /**************** Step 4) ********************************/
  // If the number of coarse unknowns falls below the minimum, create
  // a direct solver and interrupt the setup. (note, that we still
  // keep the temporary vectors of size h).
  if( settings->minSystemSize     >  SizeH             ||
      settings->numBadCoarsenings <  numBadCoarsenings ||
      Topology_->GetSizeh()       == (UInt)SizeH        ) {

    // remove topology object and eventually the penalty flags
    delete Topology_;  Topology_ = NULL;
    // setup direct solver
    const bool rValue = SetupDirectSolver( settings );

    return rValue;

  }


  /**************** Step 5) ********************************/
  // setup the smoother
  if( ! SetupSmoother(settings) ) {
    return false;
  }

  /**************** Step 6) ********************************/
  // create transfer operator
  if( NULL  == (Transfer_ = new TransferOperator<T>()) )  return false;
  if( ! Transfer_->CreateOperators(*SysMatrix_,
                                   *AuxMatrix_,
                                   *Topology_,
                                   *Agglomerate_,
                                   settings->InterpolationType,
                                   amgType_,
                                   true) ) {
    return false;
  }




  // now we do not need the topology any more
  if( false == settings->keepTopology ){delete Topology_;  Topology_ = NULL; }


  /**************** Step 7) ********************************/
  StdVector<UInt> A_H_rP, A_H_cP, B_H_rP, B_H_cP;
  StdVector<T> A_H_dP;
  StdVector<Double> B_H_dP;
  // calc the Galerkin product
  Transfer_->GalerkinProduct( A_H_rP, A_H_cP, A_H_dP, B_H_rP, B_H_cP, B_H_dP,
                              *SysMatrix_,
                              *AuxMatrix_,
                              amgType_ );

  UInt rowsA = A_H_rP.GetSize() - 1;
  UInt rowsB = (UInt)B_H_rP.GetSize() - 1;
  UInt nnzA = (UInt)A_H_rP.Last();
  UInt nnzB = (UInt)B_H_rP.Last();

  // create coarse matrices
  CRS_Matrix<T>* coarseSysMatrix = new CRS_Matrix<T>( rowsA, rowsA, nnzA);
  CRS_Matrix<Double>* coarseAuxMatrix = new CRS_Matrix<Double>( rowsB, rowsB, nnzB);

  // setup the sparsity graph
  coarseSysMatrix->SetSparsityPatternData(A_H_rP, A_H_cP, A_H_dP);
  coarseAuxMatrix->SetSparsityPatternData(B_H_rP, B_H_cP, B_H_dP);

  // sort the rows lexicographically
  coarseSysMatrix->ChangeLayout(CRS_Matrix<T>::LEX);
  coarseAuxMatrix->ChangeLayout(CRS_Matrix<Double>::LEX);


  ////////////////////////////////////////////////////
  // continue setup process recursively on next level
  // create level object for the next level
  delete nextLevel_;
  nextLevel_ = new HierarchyLevel( GetLevelID()+1, amgType_ );
  nextLevel_->SetLogging( logging_ );
  nextLevel_->InsertSysMatrix( coarseSysMatrix );
  nextLevel_->InsertAuxMatrix( coarseAuxMatrix );
  const bool setupResult = nextLevel_->Setup( settings,
                              numBadCoarsenings );
  return setupResult;
    }

/**********************************************************/


template <typename T>
bool HierarchyLevel<T>::SetupEdge( Settings* const settings,
                              Integer numBadCoarsenings){

  if( SysMatrix_ == NULL )  return false;
  if( AuxMatrix_ == NULL )  return false;
  SetLogging( settings->logging );
  if( SysMatrix_->GetCurrentLayout() != CRS_Matrix<T>::LEX_DIAG_FIRST){
	  SysMatrix_->SetCurrentLayout(CRS_Matrix<T>::LEX_DIAG_FIRST);
  }
  /*********************************************************/
  /**************** SETUP of THIS level ********************/
  /*********************************************************/
  /* Consists of the following steps
   * 1) Define the Agglomerates
   * 2) Setup patches and extractions for Arnold-Falk-Winther-Smoother (AFW-Smoother)
   * 3) Create Operators (P, R) for auxiliary matrix Bh
   * 4) Create Galerkin product for auxiliary matrix Bh
   * 5) Create coarse edges
   * 6) Setup a direct solver (PARDISO)
   * 7) Create Operators (P, R) for system matrix Ah
   * 8) Create Galerkin product for system matrix Bh
   * 9) Create auxiliary data
   * 10) Create coarse edge index node
   *
   * Note: For this algorithm we do not have a Topology of strong
   * neighbourhoods, we use the less-complex agglomeration. According
   * to diss. Stefan Reitzinger this is the better approach for edge elements.
   *
   * Perform recursively on ever deeper levels until
   * the coarsest one is reached
   */


  /**************** Step 1) ********************************/
  // setup the agglomerates
  SetupAgglomerates();

  /**************** Step 2) ********************************/
  // define the patches for the AFW-Smoother
  if( ! SetupSmoother(settings) ) {
    return false;
  }


  /**************** Step 3) ********************************/
  // create transfer operator
  if( NULL  == (Transfer_ = new TransferOperator<T>()) )  return false;
  if( ! Transfer_->CreateProlongationOperatorEdgeAux(*AuxMatrix_,
                                                     *Agglomerate_) ) {
    return false;
  }


    /**************** Step 4) ********************************/
  // calc the Galerkin product for auxiliary-matrix, needed for
  // the construction of prolongation operator of system matrix
  StdVector<UInt> B_H_rP, B_H_cP;
  StdVector<Double> B_H_dP;

  Transfer_->GalerkinProductEdgeAux(B_H_rP, B_H_cP, B_H_dP, *AuxMatrix_);

  UInt rowsB = (UInt)B_H_rP.GetSize() - 1;
  UInt nnzB = (UInt)B_H_rP.Last();
  // create coarse matrices
  CRS_Matrix<Double>* coarseAuxMatrix = new CRS_Matrix<Double>( rowsB, rowsB, nnzB);

  // setup the sparsity graph
  coarseAuxMatrix->SetSparsityPatternData(B_H_rP, B_H_cP, B_H_dP);

  // sort the rows lexicographically
  coarseAuxMatrix->ChangeLayout(CRS_Matrix<Double>::LEX_DIAG_FIRST);


  /**************** Step 5) ********************************/
  // Now we can define the coarse edges, because BH is needed
  UInt SizeH = Agglomerate_->CreateCoarseEdges(*coarseAuxMatrix);

  // check, if the coarsening gets slow
  if( settings->badCoarseningRatio * Agglomerate_->GetSizeh() <= SizeH ) {
    numBadCoarsenings++;
  } else {
    numBadCoarsenings = 0;
  }


  /**************** Step 6) ********************************/
  // If the number of coarse unknowns falls below the minimum, create
  // a direct solver and interrupt the setup. (note, that we still
  // keep the temporary vectors of size h).
  if( (UInt)settings->minSystemSize     >  SizeH             ||
      settings->numBadCoarsenings <  numBadCoarsenings ||
      Agglomerate_->GetSizeh()       == (UInt)SizeH        ) {
    // remove agglomerate object
    delete Agglomerate_;  Agglomerate_ = NULL;
    // setup direct solver
    const bool rValue = SetupDirectSolver( settings );

    return rValue;

  }

  /**************** Step 7) ********************************/
  if( ! Transfer_->CreateProlongationOperatorEdgeSys(*SysMatrix_,
                                                     edgeIndNode_,
                                                     nodeNumIndex_,
                                                     *Agglomerate_) ) {
    return false;
  }

  /**************** Step 8) ********************************/
  StdVector<UInt> A_H_rP, A_H_cP;
  StdVector<T> A_H_dP;

  Transfer_->GalerkinProductEdgeSys(A_H_rP, A_H_cP, A_H_dP, *SysMatrix_);

  UInt rowsA = (UInt)A_H_rP.GetSize() - 1;
  UInt nnzA = (UInt)A_H_rP.Last();

  // create coarse matrices
  CRS_Matrix<T>* coarseSysMatrix = new CRS_Matrix<T>( rowsA, rowsA, nnzA);


  // setup the sparsity graph
  coarseSysMatrix->SetSparsityPatternData(A_H_rP, A_H_cP, A_H_dP);

  // sort the rows lexicographically
  coarseSysMatrix->ChangeLayout(CRS_Matrix<T>::LEX_DIAG_FIRST);


  /**************** Step 9) ********************************/
  // create auxiliary data. Here we need the sizes of both, the
  // fine and the coarse system (this is the reason, why we stored it)
  // In case SizeH > settings->minSystemSize, this level will stay
  // the bottom level, and we do not need the temporary vectors of
  // size nH. We prevent their creation by passing 0 as second
  // parameter to SetupAuxData.
  UInt dim = 1;
  //UInt dim = SysMatrix_->GetNumRows() / AuxMatrix_->GetNumRows();
  if( !SetupAuxData(Agglomerate_->GetNumFineEdges(),
                    SizeH < (UInt)settings->minSystemSize ? 0 : SizeH,
                    dim) )
  {
    return false;
  }


  /**************** Step 10) ********************************/
  StdVector< StdVector< Integer> > cEdgeIndNode;
  Transfer_->GetCoarseEdgeIndNode(cEdgeIndNode);

  SysMatrix_->ChangeLayout(CRS_Matrix<T>::LEX);
  ////////////////////////////////////////////////////
  // continue setup process recursively on next level
  // create level object for the next level
  delete nextLevel_;
  nextLevel_ = new HierarchyLevel( GetLevelID()+1, amgType_ );
  nextLevel_->SetLogging( logging_ );
  nextLevel_->InsertSysMatrix( coarseSysMatrix );
  nextLevel_->InsertAuxMatrix( coarseAuxMatrix );
  nextLevel_->SetEdgeIndNode(cEdgeIndNode);
  nextLevel_->SetNodeNumIndex(Agglomerate_->GetCNodeNumIndex());
  // also insert the new edge-index-orientation maps
  const bool setupResult = nextLevel_->SetupEdge( settings,
                              numBadCoarsenings );
  return setupResult;
  }

/**********************************************************/





template <typename T>
void HierarchyLevel<T>::
Cycle( const Vector<T>& rhs,
             Vector<T>& sol,
       const Integer    num_presmoothings,
       const Integer    num_postsmoothings,
       const Integer    cycle_parameter     )
{

    // if there exists a coarser level, apply the usual
    // AMG-circus, i.e. smoothing, restriction, ...
    if( nextLevel_ ) {
        // If this is the top level, we use the passed initial
        // solution. On all coarser levels we start with the zero
        // initial guess.
        if( GetLevelID() > 0 ) {
            sol.Init();
        }

        // use the cycle parameter only on coarser levels
        const Integer lambda = GetLevelID() == 0 ? 1 : cycle_parameter;
        for( Integer i = 0; i < lambda; i++ ) {
            // v1 presmoothing steps (forward, suppress new setup)

            for( Integer v1 = 0; v1 < num_presmoothings; v1++ ) {
                Smoother_->Step( *SysMatrix_, rhs, sol,
                                 Smoother<T>::FORWARD, false );
            }
            // calculate residual (vh1_ = r = b - Ax) and restrict it
            SysMatrix_->CompRes( *vh1_, sol, rhs );


            Transfer_->Restrict( *vh1_, *vH1_, false );

            // AMG cycle on coarser level
            nextLevel_->Cycle( *vH1_, *vH2_,
                               num_presmoothings, num_postsmoothings,
                               cycle_parameter );
            // interpolate coarse correction, including addition
            Transfer_->Prolongate( *vH2_, sol, true );

            // v2 postsmoothing steps (backward, suppress new setup)
            for( Integer v2 = 0; v2 < num_postsmoothings; v2++ ) {
                Smoother_->Step( *SysMatrix_, rhs, sol,
                                 Smoother<T>::BACKWARD, false );
            }

        }
    // if there does not exist a coarser level, we must
    // solve exactly on this level
    } else {
        if( directSolver_ ) {
            //IdPrecond idprecond;
            //directSolver_->Solve( *DirSysMatrix_, idprecond, rhs, sol );

            directSolver_->Solve( *DirSysMatrix_, rhs, sol );
        } else {
        	EXCEPTION("No direct solver provided...this shouldn't happen...Implementation Error!!");

        }


    }
}

/**********************************************************/

template <typename T>
Double HierarchyLevel<T>::GetOperatorComplexity() const
{

  const HierarchyLevel<T>* level      = this;
  Double                   complexity = 0.0;

  while( level ) {
    complexity += level->SysMatrix_->GetNnz();
    level = level->GetNextLevel();
  }
  complexity /= SysMatrix_->GetNnz();

  return complexity;
}

/**********************************************************/

template <typename T>
std::ostream& HierarchyLevel<T>::Print( std::ostream& out,
    const bool color ) const
    {

  if( SysMatrix_ == NULL ) {
    out << "no matrix present in hierarchy level "
        << LevelID_ << std::endl;
    return out;
  }

  out
  << (color ? "\033[0;1m" : "")
  << "  hierarchy level "<<LevelID_
  << (color ? "\033[0m" : "")<<std::endl
  << "  system matrix : size : " << SysMatrix_->GetNumRows() << 'x'
  << SysMatrix_->GetNumCols() << std::endl
  << "                : nnz  : " << SysMatrix_->GetNnz()   << std::endl
  << "                : dens : " << (((double)SysMatrix_->GetNnz()  /
      (double)SysMatrix_->GetNumRows() )/
      (double)SysMatrix_->GetNumRows()   )
      << std::endl
      << "  complexity    : " << GetOperatorComplexity() << std::endl;

  if( nextLevel_ )  nextLevel_->Print( out );

  return out;
    }

/**********************************************************/

template <typename T>
bool HierarchyLevel<T>::SetupTopology( const Double alpha,
                                      const Double strongDiagRatio,
                                      const Double forceFineRatio)
                                      {

  // Create topology object. The constructor, we use here
  // also creates the graphs of strong dependencies S and ST.
  if( Topology_ ) { delete Topology_;  Topology_ = NULL; }

  // create new topology object
  if( NULL == (Topology_ = new Topology<Double>) ) {
    EXCEPTION(" AMG: creation of topology object failed\n");
    return false;
  }

  // setup dependency graphs
  if( 0 > Topology_->CreateDependencyGraphs(
      *AuxMatrix_,       // auxiliary matrix
      alpha,            // alpha for coarsening
      strongDiagRatio,
      forceFineRatio) ) {
    EXCEPTION(" AMG: topology setup failed\n");
    return false;
  }

  return true;
    }

/**********************************************************/


template <typename T>
UInt HierarchyLevel<T>::SetupAgglomerates(){

  // Create topology object. The constructor, we use here
  // also creates the graphs of strong dependencies S and ST.
  if( Agglomerate_ ) { delete Agglomerate_;  Agglomerate_ = NULL; }

  // create new agglomerate object
  if( NULL == (Agglomerate_ = new Agglomerate<Double>) ) {
    EXCEPTION(" AMG: creation of agglomerate object failed\n");
    return false;
  }

  // insert edges with real node numbers (correctly oriented)
  Agglomerate_->InsertEdgeList(edgeIndNode_);
  // insert the map between real node-number and index in auxiliary matrix
  Agglomerate_->InsertNodeIndex(nodeNumIndex_);

  // setup agglomerates graphs
  Integer sizeH = Agglomerate_->CreateAgglomerateGraphs(*AuxMatrix_);
  if( 0 > sizeH ) {
    EXCEPTION(" AMG: agglomerate setup failed\n");
    return false;
  }

  return sizeH;
    }
/**********************************************************/



template <typename T>
bool HierarchyLevel<T>::SetupSmoother( const Settings *const settings)
{

  // create smoother object
  if( Smoother_ ) { delete Smoother_; Smoother_ = 0x0; }
  if( amgType_ != EDGE){
    switch( settings->SmootherType ) {
    case AMG_SMOOTHER_GAUSSSEIDEL: {
      GaussSeidel<T> *gss = new GaussSeidel<T>;
      if( gss ) {
        gss->SetOmega( settings->smootherOmega );
        Smoother_ = gss;
      }
      break;
    }
    case AMG_SMOOTHER_DAMPED_JACOBI: {
      Jacobi<T> *js = new Jacobi<T>;
      if( js ) {
        js->SetOmega( settings->smootherOmega );
        Smoother_ = js;
      }
      break;
    }
    default: {
      EXCEPTION( " AMG: not supported smoother type ID "
          << settings->SmootherType<<std::endl;)
    }
    }
  }else{
    AFWSmoother<T> *afw = new AFWSmoother<T>;
    if( afw ){
      afw->SetOmega(settings->smootherOmega);
      afw->CreatePatches(*AuxMatrix_, edgeIndNode_, nodeNumIndex_);
      Smoother_ = afw;
    }
  }

  if( ! Smoother_ ) {
    EXCEPTION( " AMG: creation of smoother object failed");
    return false;
  } else {
    // setup the smoother
    if( ! Smoother_->Setup(*SysMatrix_) ) {
      EXCEPTION(" AMG: setup of smoother failed\n");
      delete Smoother_;
      Smoother_ = 0x0;
      return false;
    }
  }

  return true;
}


/**********************************************************/

template <typename T>
bool HierarchyLevel<T>::
SetupDirectSolver( const Settings* const settings )
{


  switch( settings->directSolver) {

  // NOSOLVER as setting for the direct solver in AMG means to
  // use the smoother as direct solver on the coarsest level
  case BaseSolver::NOSOLVER: {
	  EXCEPTION("AMG: Please select a direct solver!")
    //if( ! SetupSmoother(settings) )  return false;
    break;
  }


  case BaseSolver::LU_SOLVER: {
    EXCEPTION("AMG: Please use Pardiso as direct solver, LU not yet used");
    DirSysMatrix_ = SysMatrix_;
    // remember that DirSysMatrix_ must not be deleted
    // separately
    deleteDirSysMatrix_ = false;
    // create new parameter object for the direct solver
    if ( directSolverParams_ == NULL ) {
      directSolverParams_ = new PtrParamNode;
    }

    directSolver_ = GenerateDirSolverObjectAMG( *DirSysMatrix_, *directSolverParams_);
    directSolver_->Setup( *DirSysMatrix_);
    break;
  }

  case BaseSolver::PARDISO_SOLVER: {
#ifdef USE_PARDISO

    // The PARDISO solver works with CRS matrices. Therefore
    // we do not need to convert the matrix, but only point
    // to it with the pointer to the general direct solver
    // matrix.
    if(SysMatrix_->GetCurrentLayout() != CRS_Matrix<T>::LEX){
      SysMatrix_->ChangeLayout(CRS_Matrix<T>::LEX);
    }
    DirSysMatrix_ = SysMatrix_;
    // remember that DirSysMatrix_ must not be deleted
    // separately
    deleteDirSysMatrix_ = false;
    // create new parameter object for the direct solver
    if ( directSolverParams_ == NULL ) {
      directSolverParams_ = new PtrParamNode;
    }

    /* ROK: I don't think we can use the "standard" solver object because we don't have
     * SolStrategy with a pointer to a direct solver, furthermore the solver would be
     * set-up with the full matrix, here we only need the solution of the coarse problem
     *
     * Therefore use an additional method in generatesolver
     */
    directSolver_ = GenerateDirSolverObjectAMG( *DirSysMatrix_, *directSolverParams_);
    directSolver_->Setup( *DirSysMatrix_);
#else // USE_PARDISO
    std::cout<<"To use PARDISO as direct solver on the coarsest "
        "AMG level the code must be compiled with the "
        "preprocessor flags USE_PARDISO"<<std::endl;
#endif // USE_PARDISO
    break;
  }

  default: {
    EXCEPTION( "Solver type not supported as direct solver in AMG");
    break;
  }
  }

  return true;
}


template <typename T>
bool HierarchyLevel<T>::SetupAuxData( const Integer sizeh,
    const Integer sizeH, const UInt dim )
    {

  vh1_ = new Vector<T>( sizeh * dim );

  if( sizeH ) {
    vH1_ = new Vector<T>( sizeH * dim );
    vH2_ = new Vector<T>( sizeH * dim);
  }

  return true;
    }


/**********************************************************
 * HierarchyLevel<T>::Settings
 **********************************************************/

template <typename T>
HierarchyLevel<T>::Settings::Settings(  PtrParamNode& params,
     PtrParamNode& info)
: maxCoarseDepend( 6 ),
  minSystemSize( 200 ),
  numPreSmoothings( 3 ),
  numPostSmoothings( 3 ),
  CycleParameter( 1 ),
  numBadCoarsenings( 5 ),
  gamma( -1 ),
  alpha( 0.01 ),
  strongDiagRatio( 0.0 ),
  forceFineRatio( 1e-10 ),
  badCoarseningRatio( 0.6 ),
  keepTopology( false ),
  logging( true ),
  directSolver( BaseSolver::PARDISO_SOLVER ),
  InterpolationType( AMG_INTERPOLATION_CONSTANT ),
  SmootherType( AMG_SMOOTHER_DAMPED_JACOBI ),
  smootherOmega( 2.0/3.0 )
{

if( params->Has("maxCoarseDepend")) maxCoarseDepend = params->Get("maxCoarseDepend")->As<UInt>();
if( params->Has("minSystemSize")) minSystemSize = params->Get("minSystemSize")->As<UInt>();
if( params->Has("numPreSmooth")) numPreSmoothings = params->Get("numPreSmooth")->As<UInt>();
if( params->Has("numPostSmooth")) numPostSmoothings = params->Get("numPostSmooth")->As<UInt>();
if( params->Has("cycleParam")) CycleParameter = params->Get("cycleParam")->As<UInt>();
if( params->Has("numBadCoarsenings")) numBadCoarsenings = params->Get("numBadCoarsenings")->As<UInt>();
if( params->Has("gamma")) gamma = params->Get("gamma")->As<UInt>();
if( params->Has("alpha")) alpha = params->Get("alpha")->As<Double>();
if( params->Has("strongDiagRatio")) strongDiagRatio = params->Get("strongDiagRatio")->As<Double>();
if( params->Has("forceFineRatio")) forceFineRatio = params->Get("forceFineRatio")->As<Double>();
if( params->Has("badCoarseningRatio")) badCoarseningRatio = params->Get("badCoarseningRatio")->As<Double>();
if( params->Has("logging")) logging = params->Get("logging")->As<bool>();
//if( params->Has("directSolver")) directSolver = params->Get("directSolver")->As<CoupledField::BaseSolver::SolverType>();
//if( params->Has("interpolationType")) InterpolationType = params->Get("interpolationType")->As<CoupledField::AMGInterpolationType>();
//if( params->Has("smoother")) SmootherType = params->Get("smoother")->As<CoupledField::AMGSmootherType>();
      if( params->Has("gs_jacobi_omega") ) smootherOmega = params->Get("gs_jacobi_omega")->As<Double>();
}


/**********************************************************/

template <typename T>
std::ostream& HierarchyLevel<T>::Settings::
Print( std::ostream& out ) const
{
std::string interpolType;
std::string smootherType;


Enum2String<AMGInterpolationType>(InterpolationType, interpolType);
Enum2String<AMGSmootherType>(SmootherType, smootherType);


    out
    << " AMG:    MaxCoarseDependency : "<<maxCoarseDepend     <<std::endl
    << " AMG:    MinSystemSize       : "<<minSystemSize       <<std::endl
    << " AMG:    NumPreSmoothing     : "<<numPreSmoothings    <<std::endl
    << " AMG:    NumPostSmoothing    : "<<numPostSmoothings   <<std::endl
    << " AMG:    CycleParameter      : "<<CycleParameter      <<std::endl
    << " AMG:    NumBadCoarsenings   : "<<numBadCoarsenings   <<std::endl
    << " AMG:    CoarseningAlpha     : "<<alpha               <<std::endl
    << " AMG:    CoarseningGamma     : "<<gamma               <<std::endl
    << " AMG:    StrongDiagRatio     : "<<strongDiagRatio     <<std::endl
    << " AMG:    ForceFineRatio      : "<<forceFineRatio      <<std::endl
    << " AMG:    BadCoarseningRatio  : "<<badCoarseningRatio  <<std::endl
    << " AMG:    logging             : "<<(logging?"yes":"no")<<std::endl
    << " AMG:    DirectSolver        : "<< "Pardiso"          <<std::endl
    << " AMG:    InterpolationType   : "<< interpolType       <<std::endl
    << " AMG:    SmootherType        : "<< smootherType       <<std::endl
    << " AMG:    SmootherOmega       : "<<smootherOmega       <<std::endl;

  return out;
}

/**********************************************************/
} // namespace CoupledField

/**********************************************************
 * print-out operators
 **********************************************************/

namespace std {
template <typename T> std::ostream& operator <<
    ( std::ostream& out, const HierarchyLevel<T>& level ) {
  return level.Print( out ); }
}

/**********************************************************/
#ifdef DEBUG_TO_CERR
#undef debug
#endif // DEBUG_TO_CERR
/**********************************************************/
