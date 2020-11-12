#ifndef FILE_EIGEN_VALUE_DRIVER_HH
#define FILE_EIGEN_VALUE_DRIVER_HH

#include "SingleDriver.hh"
#include <fstream>
#include <iterator>

#include "Domain/Domain.hh"
#include "OLAS/solver/BaseEigenSolver.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"



namespace CoupledField {

template<class TYPE> class Vector;
class SingleVector;
//class SimState;

//! Driver class for calculating a general eigenvalue problem
class EigenValueDriver: public virtual SingleDriver {

  public:

    //! constructor
    //! \param sequenceStep current step in multisequence simulation
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    EigenValueDriver(UInt sequenceStep,
                  bool isPartOfSequence,
                  shared_ptr<SimState> state,
                  Domain *domain,
                  PtrParamNode paramNode,
                  PtrParamNode infoNode);

    //! Destructor
    ~EigenValueDriver();

    //! Initialization method
    void Init(bool restart);

    //! This method constitutes the actual driving method which controls the
    //! solution process for the problem.
    void SolveProblem();

    //! Return current time / frequency step of simulation
    UInt GetActStep(const std::string &pdename) {
      return 1;
    }

    // for eigenSolver always true, otherwise eigenmodes return with bad cast.
    virtual bool IsComplex() {
      return true;
    }

    /** Return the number of eigenmodes to be calculated.
     * @see BaseDriver::GetNumSteps() */
    unsigned int GetNumSteps() {
      return numSteps_;
    }

    // TODO implement optimization (StoreResults() is currently functionless)
    /** @see BaseDriver::StoreResults()
     * stepNum and step_val are ignored!! */
    void StoreResults(UInt stepNum, double step_val);

    // return proportionality factor (currently unused)
    double GetPropFactor(unsigned int idx) const;

    virtual void SetToStepValue(UInt stepNum, Double stepVal);

    Vector<Double> eigenValues_;
    Vector<Double> errBounds_;
    Vector<Complex> eigenValuesComplex_;
    Vector<Complex> errBoundsComplex_;

    StdVector<int> modeOrder_; // for sorting the obtained modes
    void SortModes();

  private:

    // print eigenValues to the console
    void PrintResult();

    // actually calculate eigenValues
    void CalcEigenValues();

    // calculate eigenmodes for a given eigenvalue
    void CalcModes();

    //! input parameter for input method 2, number of modes to be calculated
    unsigned int numValue_;

    //! input parameter for input method 2, shift for eigenValues
    Double shiftPoint_Real_;
    Double shiftPoint_Imag_;
    Complex shiftPoint_;

    //! input parameter for input method 1, minimum eigenvalue
    Double minVal_;

    //! input parameter for input method 1, maximum eigenvalue
    Double maxVal_;

    //! set input methode
    UInt inputMethod_;

    //! set EVP type: 1 - StandardEVP; 2-GeneralizedEVP; 3-QuadraticEVP
    UInt evpType_;

    //! Storage for Matrix A (quadratic matrix of QEVP or generalized or standard Matrix)
    FEMatrixType matrixA_;

    //! Storage for Matrix B (linear matrix of QEVP or generalized)
    FEMatrixType matrixB_;

    //! Storage for Matrix C (constant matrix of QEVP)
    FEMatrixType matrixC_;

    //! Matrix storage type
    bool isStoredSymmetric_;

    //! input parameter, true if eigenmodes shall be calculated
    bool calcModes_;

    //! needed parameter see calcMode()
    unsigned int save_step_;

    //! Number of converged EV
    UInt numSteps_ = 1;
    //! input parameter, set method of Mode sizing
    BaseEigenSolver::ModeNormalization modeNormalization_;


};

}

#endif
