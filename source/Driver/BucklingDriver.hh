#ifndef FILE_BUCKLING_DRIVER_HH
#define FILE_BUCKLING_DRIVER_HH

#include <fstream>
#include <iterator>

#include "SingleDriver.hh"
#include "Domain/Domain.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "OLAS/solver/BaseEigenSolver.hh"

namespace CoupledField {

template<class TYPE> class Vector;
class SingleVector;
//class SimState;

//! Driver class for calculating a buckling problem
class BucklingDriver: public virtual SingleDriver {

  public:

    //! constructor
    //! \param sequenceStep current step in multisequence simulation
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    BucklingDriver(UInt sequenceStep,
                  bool isPartOfSequence,
                  shared_ptr<SimState> state,
                  Domain *domain,
                  PtrParamNode paramNode,
                  PtrParamNode infoNode);

    //! Destructor
    ~BucklingDriver();

    //! Initialization method
    void Init(bool restart);

    //! This method constitutes the actual driving method which controls the
    //! solution process for the problem.
    void SolveProblem();

    //! Return current time / frequency step of simulation
    UInt GetActStep(const std::string &pdename) { return 1; }

    // for eigenSolver always true, otherwise eigenmodes return with bad cast.
    virtual bool IsComplex() { return true; }

    /** Return the number of eigenmodes to be calculated.
     * @see BaseDriver::GetNumSteps() */
    unsigned int GetNumSteps() { return numEigenValues_; }

    bool IsInverseProblem() { return isInverseProblem_; }

    /** @see BaseDriver::StoreResults()
     *  step_val has no effect! */
    unsigned int StoreResults(UInt stepNum, double step_val);

    void SetToStepValue(UInt stepNum, Double stepVal) {
      // ensure that this method is only called if simState has input
      if (! simState_->HasInput()) {
        EXCEPTION("Can only set external time step, if simulation state " << "is read from external file");
      }
      // Set current eigenvalue in the mathParser
      domain_->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "f", stepVal);
      domain_->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "step", stepNum);

    }

    SingleVector* eigenValues;
    SingleVector* errBounds;

  private:

    // print eigenValues to the console
    void PrintResult();

    // calculate eigen values and (implicitly) eigen modes
    void CalcValues();

    // store modes in algebraic system solution
    void StoreModes();

    // export eigenmodes
    void ExportModes();

    // sort modes with ascending eigenvalues
    void SortModes(bool inAbs);

    Vector<Double> GetRealPartOfLoadFactors();

    BaseEigenSolver::EigenSolverType solverType_;

    bool isInverseProblem_;

    //! input parameter for input method 2, number of modes to be calculated
    unsigned int numMode_;

    //! input parameter for input method 2, shift for eigenValues
    Double valueShift_;

    //! input parameter for input method 1, minimum eigenvalue
    Double minVal_;

    //! input parameter for input method 1, maximum eigenvalue
    Double maxVal_;

    //! set input methode
    UInt inputMethod_;

    //! Matrix storage type
    bool isStoredSymmetric_;

    //! input parameter, true if eigenmodes shall be calculated
    bool calcModes_;

    //! = eigenvalues for original problem, 1/eigenvalues for reformulated one
    SingleVector* loadFactors_;

    UInt numEigenValues_;

    // order of the sorted modes
    StdVector<int> modeOrder_;

    //! input parameter, set method of Mode sizing
    BaseEigenSolver::ModeNormalization modeNormalization_;
};

}

#endif
