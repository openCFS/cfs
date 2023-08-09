#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <DataInOut/Logging/LogConfigurator.hh>
#include <DataInOut/ParamHandling/ParamNode.hh>
#include <General/defs.hh>
#include <General/Exception.hh>
#include <MatVec/SparseOLASMatrix.hh>
#include <OLAS/algsys/SolStrategy.hh>
#include <OLAS/solver/generateEigensolver.hh>
#include <OLAS/solver/QuadraticEigenSolver.hh>
#include <Utils/StdVector.hh>
#include <Utils/Timer.hh>
#include <cassert>
#include <string>
#include <time.h>
#include "generateEigensolver.hh"

DEFINE_LOG(quadEigenSolver, "quadraticEigenSolver")

namespace CoupledField
{
  QuadraticEigenSolver::QuadraticEigenSolver(shared_ptr<SolStrategy> strat, PtrParamNode xml, PtrParamNode solverList, PtrParamNode precondList, PtrParamNode eigenInfo, PtrParamNode eSolverList)
      : BaseEigenSolver(strat, xml, solverList, precondList, eigenInfo)
  {
    this->numFreq_ = 0;
    this->freqShift_ = 0;
    xml_ = xml;

    // read XML settings
    string solver_str = xml_->Get("generalisedEigenSolver")->Get("id")->As<std::string>(); //string of the chosen eigensolver
    xml_->GetValue("linearisation", linmode, ParamNode::INSERT);

    //generate a solver object according to the chosen GEVP-solver
    solverForGeneralisedEVP = NULL;
    solverForGeneralisedEVP = GenerateEigenSolverObject(solStrat_,  eSolverList, solverList, precondList, eigenInfo, solver_str);

    A_ = B_ = NULL;
  }

  QuadraticEigenSolver::~QuadraticEigenSolver()
  {
    LOG_DBG(quadEigenSolver) << "destructor";
    xml_ = NULL;
    solverForGeneralisedEVP = NULL;
    A_ = B_ = NULL;
  }

  void QuadraticEigenSolver::ToInfo()
  {
    //todo implement something
  }

  void QuadraticEigenSolver::Setup(const BaseMatrix &K, const BaseMatrix &C, const BaseMatrix &M)
  {
    // The timer shows how long the linearisation of the problem took.
    clock_t t;
    t = clock();
    // check the matrix properties
    bool isRealK, isStoredSymmetricK, isRealC, isStoredSymmetricC, isRealM, isStoredSymmetricM;
    CheckMatrix(isRealK, isStoredSymmetricK, K);
    CheckMatrix(isRealC, isStoredSymmetricC, C);
    CheckMatrix(isRealM, isStoredSymmetricM, M);
    assert(isRealK == isRealC);
    assert(isRealK == isRealM);
    bool isReal = isRealK && isRealC && isRealM;

    assert(isStoredSymmetricK == isStoredSymmetricM);
    assert(isStoredSymmetricK == isStoredSymmetricC);
    bool isSymmetric = isStoredSymmetricK && isStoredSymmetricC && isStoredSymmetricM;

    if (isReal)
    {
      double dummy = 0.0;
      doLinearisation(dummy, K, C, M, isSymmetric);
    }
    else
    {
      Complex dummy = Complex(0.0);
      doLinearisation(dummy, K, C, M, isSymmetric);
    }
    this->SetProblemType(*A_);
    solverForGeneralisedEVP->Setup(*A_, *B_);

    // Writes the time the linearisation took in seconds in the console.
    t = clock() - t;
    std::cout << "\n";
    std::cout << "Time needed for the linearisation: " << (float)t / CLOCKS_PER_SEC << " sec"
              << "\n";
    std::cout << "\n";
    //A_->Export("A");
    //B_->Export("B");
  }

  template <typename T>
  void QuadraticEigenSolver::doLinearisation(T dummy, const BaseMatrix &K, const BaseMatrix &C, const BaseMatrix &M, bool isSymmetric)
  {

    // Determine dimension of the problem:
    UInt n = K.GetNumRows();
    UInt nn = 2 * n;

    // check if all matrices have the same number of non-zeros (shouls always be the case)
    UInt nnzK = GetNnz(K);
    UInt nnzC = GetNnz(C);
    UInt nnzM = GetNnz(M);
    assert(nnzK == nnzM);
    assert(nnzK == nnzC);
    UInt nnz = nnzK;

    // assign the correct size to row, column and data vectors:
    StdVector<UInt> M_r, M_c, C_r, C_c, K_r, K_c;
    M_r.Resize(n + 1);
    C_r.Resize(n + 1);
    K_r.Resize(n + 1);
    M_c.Resize(nnz);
    C_c.Resize(nnz);
    K_c.Resize(nnz);
    StdVector<T> M_v, C_v, K_v; // data
    M_v.Resize(nnz);
    C_v.Resize(nnz);
    K_v.Resize(nnz);

    //define row, column and data vectors for the A and B matrices created in the linearisation process:
    StdVector<UInt> A_r;
    StdVector<UInt> A_c;
    StdVector<T> A_v;
    StdVector<UInt> B_r;
    StdVector<UInt> B_c;
    StdVector<T> B_v;

    // extract the row, column and data vectors of the system matrices
    if (isSymmetric)
    {
      //dynamic_cast< SCRS_Matrix<T>* >(K)
      CRS_Matrix<T> M_CRS = CRS_Matrix<T>(dynamic_cast<const SCRS_Matrix<T> &>(M));
      M_CRS.GetVectors(&M_c, &M_r, &M_v);
      CRS_Matrix<T> C_CRS = CRS_Matrix<T>(dynamic_cast<const SCRS_Matrix<T> &>(C));
      C_CRS.GetVectors(&C_c, &C_r, &C_v);
      CRS_Matrix<T> K_CRS = CRS_Matrix<T>(dynamic_cast<const SCRS_Matrix<T> &>(K));
      K_CRS.GetVectors(&K_c, &K_r, &K_v);
    }
    else
    {
      dynamic_cast<const CRS_Matrix<T> &>(M).GetVectors(&M_c, &M_r, &M_v);
      dynamic_cast<const CRS_Matrix<T> &>(C).GetVectors(&C_c, &C_r, &C_v);
      dynamic_cast<const CRS_Matrix<T> &>(K).GetVectors(&K_c, &K_r, &K_v);
    }
    LOG_DBG(quadEigenSolver) << "K: " << n << "x" << n << " with " << nnz << " non-zero elements\n";
    LOG_DBG(quadEigenSolver) << "Rows   :" << K_r.ToString() << "\n";
    LOG_DBG(quadEigenSolver) << "Columns:" << K_c.ToString() << "\n";
    LOG_DBG(quadEigenSolver) << "Data   :" << K_v.ToString() << "\n";

    LOG_DBG(quadEigenSolver) << "M: " << n << "x" << n << " with " << nnz << " non-zero elements\n";
    LOG_DBG(quadEigenSolver) << "Rows   :" << M_r.ToString() << "\n";
    LOG_DBG(quadEigenSolver) << "Columns:" << M_c.ToString() << "\n";
    LOG_DBG(quadEigenSolver) << "Data   :" << M_v.ToString() << "\n";

    // case differentiation depending on which linearisation scheme was chosen:
    UInt A_nnz, B_nnz;
    if (linmode == "firstCompanion")
    {
      // MATRIX A= [ 0  I]
      //           [-K -C]:

      // total number of non zero elements in matrix A and B:
      A_nnz = nnzC + nnzK + n;
      B_nnz = nnzM + n;

      // total number of non zero elements in matrix A:
      A_r.Resize(nn + 1);
      A_c.Resize(A_nnz);
      A_v.Resize(A_nnz);
      //A.SetSize(nn,nn,A_nnz);

      // Set I:
      A_r[0] = 0;
      for (UInt i = 0; i < n; i++)
      {
        A_v[i] = 1;
        A_r[i + 1] = i + 1;
        A_c[i] = i + n;
      }

      // Set -K and -C: Since there are two matrices in the same row, counters are needed to
      // determine the number of values in each row for the row vector:

      UInt kvc = 0;  // K value counter
      UInt tvc = n;  // total value counter
      UInt cvc = 0;  // C value counter
      UInt tvpr = 0; //total values per row counter
      UInt nkv = 0;  //number of K values per row
      UInt ncv = 0;  //number of C values per row
      for (UInt i = n; i < nn; i++)
      {
        tvpr = 0;
        nkv = K_r[i - n + 1] - K_r[i - n]; //number of K values per row
        for (UInt j = 0; j < nkv; j++)
        {
          A_v[tvc] = -K_v[kvc];
          A_c[tvc] = K_c[kvc];
          kvc++;
          tvc++;
          tvpr++;
        }
        ncv = C_r[i - n + 1] - C_r[i - n]; //number of C values per row
        for (UInt k = 0; k < ncv; k++)
        {
          A_v[tvc] = -C_v[cvc];
          A_c[tvc] = C_c[cvc] + n;
          cvc++;
          tvc++;
          tvpr++;
        }
        A_r[i + 1] = A_r[i] + tvpr;
      }

      // MATRIX B = [I 0]
      //            [0 M]:

      // set correct dimensions:
      B_r.Resize(nn + 1);
      B_c.Resize(B_nnz);
      B_v.Resize(B_nnz);
      // Set I:
      B_r[0] = 0;
      for (UInt i = 0; i < n; i++)
      {
        B_v[i] = 1;
        B_r[i + 1] = i + 1;
        B_c[i] = i;
      }
      // Set M:
      for (UInt i = n; i < nn; i++)
      {
        B_r[i + 1] = M_r[i - n + 1] - M_r[i - n] + B_r[i];
      }
      for (UInt i = 0; i < nnz; i++)
      {
        B_v[i + n] = M_v[i];
        B_c[i + n] = M_c[i] + n;
      }
    }
    //=====================================================================================
    else if (linmode == "secondCompanion")
    {
      // The layout of A and B are different, building the matrices works the same.

      // total number of non zero elements in matrix A and B:
      A_nnz = nnzK + n;
      B_nnz = nnzM + n + nnzC;

      // MATRIX A = [-K  0]
      //            [ 0 -I]:

      A_r.Resize(nn + 1);
      A_c.Resize(A_nnz);
      A_v.Resize(A_nnz);

      A_r[0] = 0;
      for (UInt i = 0; i < n; i++)
      {
        A_r[i + 1] = K_r[i + 1] - K_r[i] + A_r[i];
      }
      for (UInt i = 0; i < nnz; i++)
      {
        A_v[i] = -K_v[i];
        A_c[i] = K_c[i];
      }

      for (UInt i = n; i < nn; i++)
      {
        A_r[i + 1] = A_r[i] + 1;
        A_c[i + nnz - n] = i;
        A_v[i + nnz - n] = -1;
      }

      // MATRIX B = [ C M]
      //            [-I 0]:

      B_r.Resize(nn + 1);
      B_c.Resize(B_nnz);
      B_v.Resize(B_nnz);

      UInt mvc = 0;  // M value counter
      UInt tvc = 0;  // total value counter
      UInt cvc = 0;  // C value counter
      UInt tvpr = 0; //total values per row counter
      UInt nmv = 0;  //number of m values per row
      UInt ncv = 0;  //number of c values per row
      B_r[0] = 0;
      for (UInt i = 0; i < n; i++)
      {
        tvpr = 0;
        ncv = C_r[i + 1] - C_r[i];
        for (UInt j = 0; j < ncv; j++)
        {
          B_v[tvc] = C_v[mvc];
          B_c[tvc] = C_c[mvc];
          mvc++;
          tvc++;
          tvpr++;
        }
        nmv = M_r[i + 1] - M_r[i];
        for (UInt k = 0; k < nmv; k++)
        {
          B_v[tvc] = M_v[cvc];
          B_c[tvc] = M_c[cvc] + n;
          cvc++;
          tvc++;
          tvpr++;
        }
        B_r[i + 1] = B_r[i] + tvpr;
      }
      for (UInt i = n; i < nn; i++)
      {
        B_v[i + B_nnz - nn] = -1;
        B_r[i + 1] = B_r[i] + 1;
        B_c[i + B_nnz - nn] = i - n;
      }
    }
    else
    {
      EXCEPTION("This case should not be possible according to the xml-schema");
    }

    LOG_DBG(quadEigenSolver) << "A:\n";
    LOG_DBG(quadEigenSolver) << "Rows   :" << A_r.ToString() << "\n";
    LOG_DBG(quadEigenSolver) << "Columns:" << A_c.ToString() << "\n";
    LOG_DBG(quadEigenSolver) << "Data   :" << A_v.ToString() << "\n";

    LOG_DBG(quadEigenSolver) << "B:\n";
    LOG_DBG(quadEigenSolver) << "Rows :" << B_r.ToString() << "\n";
    LOG_DBG(quadEigenSolver) << "Columns:" << B_c.ToString() << "\n";
    LOG_DBG(quadEigenSolver) << "Data   :" << B_v.ToString() << "\n";
    // create matrix and set data
    //CRS_Matrix<T> A(nn,nn,A_nnz);
    //A.SetSparsityPatternData(A_r,A_c,A_v);

    A_ = new CRS_Matrix<T>(nn, nn, A_nnz);
    dynamic_cast<CRS_Matrix<T> *>(A_)->SetSparsityPatternData(A_r, A_c, A_v);
    LOG_DBG(quadEigenSolver) << "A = " << A_->ToString() << std::endl;

    B_ = new CRS_Matrix<T>(nn, nn, B_nnz);
    dynamic_cast<CRS_Matrix<T> *>(B_)->SetSparsityPatternData(B_r, B_c, B_v);
    LOG_DBG(quadEigenSolver) << "B = " << B_->ToString() << std::endl;
  }

  void QuadraticEigenSolver::
      CalcEigenValues(BaseVector &sol, BaseVector &err, UInt N, Double shiftPoint)
  {
    LOG_DBG(quadEigenSolver) << "CalcEigenValues: START";
    shared_ptr<Timer> timer = info_->Get("QuadraticEigenSolver/timer")->AsTimer();
    timer->Start();

    //The timer shows how long it took the chosen Solver to calculate the Eigenvalues
    clock_t t;
    t = clock();

    // Pass the Problem to the "CalcEigenValues" function of the chosen solver:
    solverForGeneralisedEVP->CalcEigenValues(sol, err, N, shiftPoint);

    t = clock() - t;
    std::cout << "Time needed to solve the Eigenvalue problem: " << (float)t / CLOCKS_PER_SEC << " sec\n";

    timer->Stop();
    LOG_DBG(quadEigenSolver) << "CalcEigenValues: END";
  }

  void QuadraticEigenSolver::
      CalcEigenValues(BaseVector &sol, BaseVector &err, Double minVal, Double maxVal)
  {
    LOG_DBG(quadEigenSolver) << "CalcEigenValues: START";
    shared_ptr<Timer> timer = info_->Get("QuadraticEigenSolver/timer")->AsTimer();
    timer->Start();

    //The timer shows how long it took the chosen Solver to calculate the Eigenvalues
    clock_t t;
    t = clock();

    // Pass the Problem to the "CalcEigenValues" function of the chosen solver:
    solverForGeneralisedEVP->CalcEigenValues(sol, err, minVal, maxVal);

    t = clock() - t;
    std::cout << "Time needed to solve the Eigenvalue problem: " << (float)t / CLOCKS_PER_SEC << " sec\n";

    timer->Stop();
    LOG_DBG(quadEigenSolver) << "CalcEigenValues: END";
  }

  void QuadraticEigenSolver::
      CalcEigenValues(BaseVector &sol, BaseVector &err)
  {
    LOG_DBG(quadEigenSolver) << "CalcEigenValues: START";
    shared_ptr<Timer> timer = info_->Get("QuadraticEigenSolver/timer")->AsTimer();
    timer->Start();

    //The timer shows how long it took the chosen Solver to calculate the Eigenvalues
    clock_t t;
    t = clock();

    // Pass the Problem to the "CalcEigenValues" function of the chosen solver:
    solverForGeneralisedEVP->CalcEigenValues(sol, err);

    t = clock() - t;
    std::cout << "Time needed to solve the Eigenvalue problem: " << (float)t / CLOCKS_PER_SEC << " sec\n";

    timer->Stop();
    LOG_DBG(quadEigenSolver) << "CalcEigenValues: END";
  }

  void QuadraticEigenSolver::GetEigenMode(unsigned int modeNr, Vector<Complex> &mode, bool right)
  {
    // Pass the Problem to the "GetEigenMode" function of the chosen solver
    solverForGeneralisedEVP->GetEigenMode(modeNr, mode, right);
  }
}
