#ifndef FILE_EIGENVALUE_DRIVER_HH
#define FILE_EIGENVALUE_DRIVER_HH

#include "SingleDriver.hh"
#include <fstream>
#include <iterator>

#include "Domain/Domain.hh"
#include "OLAS/solver/BaseEigenSolver.hh"

namespace CoupledField {

template <class TYPE> class Vector;
class SingleVector;
//class SimState;
  
  //! Driver class for calculating a general eigenvalue problem
  class EigenFrequencyDriver : public virtual SingleDriver {

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
    
    /** part of bloch mode analysis which needs to be excitation by excitation in the optimization case.
      Otherwise called by SolveProblem() */
    void ComputeBlochWaveVector(int wave_vector_step);

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

    /** We need to set it for optimization to get the proper stiffness matrices during function gradient evaluation.
     * Note that this is done after all systems are solved and stored. */
    void SetCurrentWaveVector(unsigned int index) { current_wave_vector_ = wave_vectors[index]; }

    /** we need to store current_wave_vector, find the index :(. Not very fast! */
    unsigned int GetCurrentWaveVectorIndex() const;

    void Eig2Freq(Vector<Double>& eigRe, Vector<Double> & freq ) {
        freq.Resize( eigRe.GetSize() );
        Double twoPi = 8.0*atan(1.0);
        for (UInt i=0; i < eigRe.GetSize(); i++) {
          freq[i] = sqrt(std::abs(eigRe[i]))/twoPi; // take the absolute value only to avoid problems with very small but negative EVs
        }
    }
    void QuadEig2FreqDamp(Vector<Double>& eigRe, Vector<Double>& eigIm, Vector<Double> & freq_undamped, Vector<Double> &freq_damped, Vector<Double> & damp ) {
        assert( eigRe.GetSize() == eigIm.GetSize() );
        freq_undamped.Resize( eigRe.GetSize() );
        freq_damped.Resize( eigRe.GetSize() );
        damp.Resize( eigRe.GetSize() );
        Double twoPi = 8.0*atan(1.0);
        for (UInt i=0; i < eigRe.GetSize(); i++) {
            double eigRatio = eigIm[i]/eigRe[i];
            damp[i] = sqrt( 1.0/(1.0+eigRatio*eigRatio) );
            freq_undamped[i] = std::abs(eigRe[i])/damp[i]/twoPi;
            freq_damped[i] = std::abs(eigIm[i])/twoPi;
        }
    }
    void Eig2FreqDamp(Vector<Double>& eigRe, Vector<Double>& eigIm, Vector<Double> & freq, Vector<Double> & damp ) {
        assert( eigRe.GetSize() == eigIm.GetSize() );
        freq.Resize( eigRe.GetSize() );
        damp.Resize( eigRe.GetSize() );
        Double twoPi = 8.0*atan(1.0);
        for (UInt i=0; i < eigRe.GetSize(); i++) {
            double eigRatio = eigIm[i]/eigRe[i];
            damp[i] = sqrt( 1.0/(1.0-eigRatio*eigRatio) );
            freq[i] = std::abs(eigRe[i])/damp[i]/twoPi;
        }
    }

    void Eig2FreqDamp(Vector<Complex>& eig, Vector<Double> & freq, Vector<Double> & damp ) {
        //assert( eigRe.GetSize() == eigIm.GetSize() );
        freq.Resize( eig.GetSize() );
        damp.Resize( eig.GetSize() );
        Double twoPi = 8.0*atan(1.0);
        for (UInt i=0; i < eig.GetSize(); i++) {
            // no idea if this is correct in general ...
          Complex lam = sqrt(eig[i]);
          double omega = lam.real();
          damp[i] = lam.imag();
          freq[i] = omega/twoPi;
          /*double omega = sqrt(eig[i].real()); # alternative implementation from sharedopt
          damp[i] = eig[i].imag()/omega;
          freq[i] = omega/twoPi;*/
        }
    }
    void SortModes(){
        modeOrder_.Resize(frequency_.GetSize());
        std::size_t n(0);
        // allocate modeOrder_
        std::generate(std::begin(modeOrder_), std::end(modeOrder_), [&]{ return n++; });
        // sort it by requency
        std::sort(std::begin(modeOrder_),std::end(modeOrder_),[&](int i1, int i2) { return frequency_[i1] < frequency_[i2]; } );
    }

    /** Return the number of eigenfrequencies to be calculated. Not the number of wave_vectors!!
     * @see BaseDriver::GetNumSteps() */
    unsigned int GetNumSteps() { return numFreq_; }

    static unsigned int GetNumModes(PtrParamNode node);

    /** Helper for ContextManager in multi sequence optimzation */
    static unsigned int GetNumBlochWave(PtrParamNode node);

    /** @see BaseDriver::StoreResults()
     * stepNum and step_val are ignored!! */
    void StoreResults(UInt stepNum, double step_val);


    /** eigenFreqs might be complex in the quadaratic, then we need to extract the real frequency by from the imaginary part
     * @param mode index within eigenFreq (0-based)
     * @return might be negative! */
    double GetFrequency(unsigned int idx) const;

    /** is the real part of the quadratic eigenFrequency
     * @return 0.0 if not quadratic */
    double GetDamping(unsigned int idx) const;

    /** create header for .bloch.dat file. For cfs -d the iteration is added to the filename */
    void SetupBlochPlot();

    /** for multi sequence optimization we need some information before driver instantiation */
    static bool DoBloch(PtrParamNode node) { return node->Has("bloch"); }

    void SetToStepValue(UInt stepNum, Double stepVal )  {
      // ensure that this method is only called if simState has input
      if( false ) { // ! simState_->HasInput()
        EXCEPTION( "Can only set external time step, if simulation state "
                << "is read from external file" );
      }

      //actFreqStep_ = stepNum;;
      //actFreq_ = stepVal;

      // Set current frequency value in the mathParser
      domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER, "f", stepVal );
      domain_->GetMathParser()->SetValue( MathParser::GLOB_HANDLER, "step", stepNum );

    }

    /** the resent calculated eigenvalues. Might be complex, @see GetFrequency(). Corresponds with errBounds_ */
    SingleVector* eigenFreqs; //ToDo: remove due to new structure -> frequency_

    //! the recent calculated eigen values, might be complex
    SingleVector* eigenVals_;

    Vector<Double> frequency_; // frequency
    Vector<Double> dampingRatio_; // damping ratio
    Vector<Double> eigsRe_; // real part
    Vector<Double> eigsIm_; // imag part

    /** this is the list of wave vectors we have to process.
     * Obtained arbitrary */
    StdVector<Vector<double> > wave_vectors;

    StdVector<int> modeOrder_; // for sorting the obtained modes

  private:

    /** fill wave_vectors_ */
    void FillWaveVectors(PtrParamNode bloch_pn);

    /** Prints info.xml, console and bloch.dat output. Handles if we are in the optimization case.
     * Does NOT write stuff to output files, This is done via StoreResults() */
    void PrintResult(int wave_vector_step = -1);

    /** the different bloch ibz types. SYMMETRIC is the triangle/tetrahedron - only valid vor symmetric designs! */
    typedef enum { NO_IBZ, SYMMETRIC, QUADRANT, HORIZONZAL, FULL } Boundary;

    Enum<Boundary> boundary;

    /** the bloch IBZ type we have */
    Boundary boundary_;

    /** corresponds with eigenFreqs */
    Vector<Double> errBounds_;

    //! Flag indicating, if a quadratic eigenvalue problem is to
    //! be solved
    bool isQuadratic_;

    /** perform complex bloch analysis to show band-gaps in the irreducible Brillouin zone */
    bool isBloch_;
    
    bool sort_;

    //! Number of eigenfrequencies to be calculated
    unsigned int numFreq_;

    // In case we do Bloch mode analysis, the number of steps (to be multiplied by numFreq_)
    unsigned int blochSteps_;

    //! Shift for eigenvalues
    Double freqShift_;

    //! minimum Eigenvalue (to define search interval for feast)
    Double minVal_;

    //! maximum Eigenvalue (to define search interval for feast)
    Double maxVal_;

    //! Flag for writing the eigenmods into the file
    bool writeModes_;

    /** This is the current wave vector index, a copy form the an wave_vectors_ entry.
     * We may not store only the index as StrainOperatorBloch2D stores the pointer. */
    Vector<double> current_wave_vector_;

    /** here we output the bloch mode data for direct plotting */
    std::ofstream bloch_plot_;

    /** the filename for bloch_plot_ */
    std::string bloch_name_;

    /** do we do ibz? Then bloch.plot will repeat the first step as final step. */
    bool ibz_;

    /** store the first plot.dat line to be repeated in the ibz_ case as last step */
    std::string first_plot_line_;

    /** the step number is complicated with bloch and or optimization. Count by ourselves here! */
    unsigned int save_step_;

    bool eigenValuesAreReal_;

    BaseEigenSolver::ModeNormalization modeNormalization_;
  };

}

#endif
