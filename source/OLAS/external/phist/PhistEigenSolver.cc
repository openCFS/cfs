#include <limits>
#include <string.h>

#include "MatVec/StdMatrix.hh"
#include "MatVec/generatematvec.hh"
#include "Utils/Timer.hh"
#include "Domain/Domain.hh"
#include "Driver/BaseDriver.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "OLAS/algsys/SolStrategy.hh"

#include "OLAS/precond/generateprecond.hh"
#include "OLAS/precond/BasePrecond.hh"
#include "OLAS/solver/generatesolver.hh"
#include "OLAS/solver/BaseSolver.hh"

#include "PhistEigenSolver.hh"

#include <ghost.h>
#include <phist_config.h>
#include <phist_tools.h>


DECLARE_LOG(pes)
DEFINE_LOG(pes, "phistEigenSolver")

namespace CoupledField {

  PhistEigenSolver::PhistEigenSolver( shared_ptr<SolStrategy> strat,
                                        PtrParamNode xml,
                                        PtrParamNode solverList,
                                        PtrParamNode precondList,
                                        PtrParamNode eigenInfo)
    : BaseEigenSolver( strat, xml, solverList, precondList, eigenInfo )
  {
    matrixA_      = NULL;
    matrixB_      = NULL;
    A_ = NULL;
    B_ = NULL;
    x_ = NULL;
    y_ = NULL;
    xml_          = xml;
  }

  PhistEigenSolver::~PhistEigenSolver()
  {
    ghost_sparsemat_destroy(A_);
    ghost_densemat_destroy(x_);
    ghost_densemat_destroy(y_);

    matrixA_ = NULL;
    matrixB_ = NULL;


    ghost_finalize();


  }
  void PhistEigenSolver::Setup(const BaseMatrix & mat, UInt numFreq, Double freqShift, bool sort)
  {
    LOG_DBG(pes) << "PES:S(stiff)";
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    isBloch_ = false;
    sort_ = false;

    // Save frequency parameters
    numFreq_ = numFreq;
    freqShift_ = freqShift;

    ToInfo();
  }

  /** taken from essex/ghost/test/minimal.c */
  int Diag(ghost_gidx row, ghost_lidx *rowlen, ghost_gidx *col, void *val, __attribute__((unused)) void *arg)
  {
      *rowlen = 1;
      col[0] = row;
      ((double *)val)[0] = (double)(row+1);

      return 0;
  }
  
  void PhistEigenSolver::Setup(const BaseMatrix& stiffMat, const BaseMatrix& massMat, UInt numFreq, Double freqShift, bool sort, bool bloch)
  {
    LOG_DBG(pes) << "PES:S(stiff, mass)";
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    isBloch_ = bloch;
    sort_ = sort;

    // Save frequency parameters
    numFreq_ = numFreq;
    freqShift_ = freqShift;

    // Copy matrix references and determine size of system
    matrixA_ = & dynamic_cast<const StdMatrix&>(stiffMat);

    // bloch works only for non-symmetric matrices as the stiffness matrix needs to be Hermitian.
    // At least Pardiso can be used with <pardiso> <hermitean>yes</hermitean> </pardiso>
    matrixB_ = & dynamic_cast<const StdMatrix&>(massMat);

    unsigned int size = matrixA_->GetNumRows();


    ghost_densemat_traits vtraits = GHOST_DENSEMAT_TRAITS_INITIALIZER;

    ghost_sparsemat_traits mtraits = GHOST_SPARSEMAT_TRAITS_INITIALIZER;
    mtraits.datatype = (ghost_datatype)(GHOST_DT_REAL|GHOST_DT_DOUBLE);
    mtraits.flags = GHOST_SPARSEMAT_SAVE_ORIG_COLS; // needed for printing if more than one process

    //GHOST_CALL_RETURN(ghost_init(argc,argv));

    // create matrix source
    ghost_sparsemat_src_rowfunc matsrc = GHOST_SPARSEMAT_SRC_ROWFUNC_INITIALIZER;
    matsrc.func = Diag;
    matsrc.maxrowlen = 1;
    matsrc.gnrows = 4; // N

    ghost_error err = GHOST_SUCCESS;
    // create sparse matrix A from row-wise source function
    err = ghost_sparsemat_create(&A_, NULL, &mtraits, 1);
    assert(err == GHOST_SUCCESS);
    err = ghost_sparsemat_create(&A_, NULL, &mtraits, 1);
    err = ghost_sparsemat_init_rowfunc(A_,&matsrc,MPI_COMM_WORLD,0.);

    // create and initialize input vector x and output vector y
    err = ghost_densemat_create(&x_, ghost_context_max_map(A_->context), vtraits);
    err = ghost_densemat_create(&y_, ghost_context_max_map(A_->context), vtraits);
    err = ghost_densemat_init_rand(x_);      // x = random
    double zero = 0;
    err = ghost_densemat_init_val(y_,&zero); // y = 0

    // compute y = A*x
    err = ghost_spmv(y_,A_,x_,GHOST_SPMV_OPTS_INITIALIZER);

    // print y, A and x
    char *Astr, *xstr, *ystr;
    err = ghost_sparsemat_string(&Astr,A_,1);
    err = ghost_densemat_string(&xstr,x_);
    err = ghost_densemat_string(&ystr,y_);
    printf("%s\n=\n%s\n*\n%s\n",ystr,Astr,xstr);
    free(Astr);
    free(xstr);
    free(ystr);

    ToInfo();
  }

  void PhistEigenSolver::Setup(const BaseMatrix & stiffMat, const BaseMatrix & massMat, const BaseMatrix & dampMat,
                                UInt numFreq, double freqShift, bool sort)
  {
    LOG_DBG(pes) << "PES:S(stiff, mass, damp)";

    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = true;
    isBloch_ = false;
    sort_ = false;

  }

  void PhistEigenSolver::ToInfo()
  {
    PtrParamNode setup = info_->Get(ParamNode::HEADER);
    setup->Get("quadratic")->SetValue(isQuadratic_);
    setup->Get("bloch")->SetValue(isBloch_);
  }




  void PhistEigenSolver::CalcEigenFrequencies(BaseVector &sol, BaseVector &err)
  {
    assert(false);
  }

  void PhistEigenSolver::CalcConditionNumber(const BaseMatrix& mat, Double& condNumber, Vector<Double>& evs, Vector<Double>& err )
  {
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    // NOTE: Hard coded as true!!!

  }

  void PhistEigenSolver::GetEigenMode(UInt modeNr, Vector<Complex> & mode)
  {

  }

  void PhistEigenSolver::GetComplexEigenMode(UInt modeNr, Vector<Complex>& mode)
  {
    // in bloch mode case the same as GetEigenMode,
    // in quadratic case the modes have internally double size and we want the upper half
  }

}
