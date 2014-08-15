#ifndef FILE_EIGENVALUE_DRIVER_HH
#define FILE_EIGENVALUE_DRIVER_HH

#include "SingleDriver.hh"
#include <fstream>

namespace CoupledField {

template <class TYPE> class Vector;
class SingleVector;
  
  //! Driver class for calculating a general eigenvalue problem
  class EigenFrequencyDriver : public SingleDriver {

  public:

    //! constructor
    //! \param sequenceStep current step in multisequence simulation
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    EigenFrequencyDriver(  UInt sequenceStep,
                           bool isPartOfSequence,
                           shared_ptr<SimState> state,
                           Domain* domain,
                           PtrParamNode paramNode, PtrParamNode infoNode );
    
    //! Destructor 
    ~EigenFrequencyDriver();

    //! Initialization method
    void Init(bool restart);
  
    //! Main method solution method

    //! This method constitutes the actual driving method which controls the
    //! solution process for the problem.
    void SolveProblem();
    
    //! Return current time / frequency step of simulation
    UInt GetActStep( const std::string& pdename ) { return 1;}

    /** Helper method which determines if an AnalyisType is complex.
     * This means that the solution vector is complex to allow to define velocity and acceleration.
     * The matrix itself shall only be complex for isQudratic or isBloc_ */
    virtual bool IsComplex() { return true; };

    /** @see BaseDriver::DoBlochModeEigenfrequency() */
    bool DoBlochModeEigenfrequency() const { return isBloch_; };

    /** @see current_wave_vector_ */
    Vector<double>& GetCurrentWaveVector() { assert(current_wave_vector_.GetSize() > 0); return current_wave_vector_; }

  private:

    /** fill wave_vectors_ */
    void FillWaveVectors(PtrParamNode bloch_pn);

    /** This is the templated form to handle the general and quadratic case */
    template <class T>
    void PrintResult(SingleVector* frequencies, Vector<Double>& bounds,
                     ResultHandler* resHandler, UInt numConverged, int wave_vector_step = -1);
    
    //! Flag indicating, if a quadratic eigenvalue problem is to
    //! be solved
    bool isQuadratic_;

    /** perform complex bloch analysis to show band-gaps in the irreducible Brillouin zone */
    bool isBloch_;
    
    //! Number of eigenfrequencies to be calculated
    UInt numFreq_;

    // In case we do Bloch mode analysis, the number of steps (to be multiplied by numFreq_)
    UInt blochSteps_;

    //! Shift for eigenvalues
    Double freqShift_;

    //! Flag for writing the eigenmods into the file
    bool writeModes_;

    /** This is the current wave vector, a copy form the an wave_vectors_ entry */
    Vector<double> current_wave_vector_;

    /** this is the list of wave vectors we have to process.
     * Obtained arbitrary */
    StdVector<Vector<double> > wave_vectors_;

    /** here we output the bloch mode data for direct plotting */
    std::ofstream bloch_plot_;

    /** do we do ibz? Then bloch.plot will repeat the first step as final step. */
    bool ibz_;

    /** store the first plot.dat line to be repeated in the ibz_ case as last step */
    std::string first_plot_line_;
  };

}

#endif
