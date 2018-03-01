#ifndef FILE_EIGENVALUE_DRIVER_HH
#define FILE_EIGENVALUE_DRIVER_HH

#include "SingleDriver.hh"
#include <fstream>
#include <iterator>

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
            freq[i] = sqrt(eigRe[i])/twoPi;
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
            freq[i] = eigRe[i]/damp[i]/twoPi;
        }
    }
    void Eig2FreqDamp(Vector<Complex>& eig, Vector<Double> & freq, Vector<Double> & damp ) {
        //assert( eigRe.GetSize() == eigIm.GetSize() );
        freq.Resize( eig.GetSize() );
        damp.Resize( eig.GetSize() );
        Double twoPi = 8.0*atan(1.0);
        for (UInt i=0; i < eig.GetSize(); i++) {
            // no idea if this is correct in general ...
            double omega = sqrt(eig[i].real());
            damp[i] = eig[i].imag()/omega;
            freq[i] = omega/twoPi;
        }
    }
    void SortModes(){
        ModeOrder.Resize(Frequency.GetSize());
        std::size_t n(0);
        // allocate ModeOrder
        std::generate(std::begin(ModeOrder), std::end(ModeOrder), [&]{ return n++; });
        // sort it by requency
        std::sort(std::begin(ModeOrder),std::end(ModeOrder),[&](int i1, int i2) { return Frequency[i1] < Frequency[i2]; } );
        std::cout << " ModeOrder : "<<ModeOrder.ToString()<<"\n";
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

    /** the resent calculated eigenvalues. Might be complex, @see GetFrequency(). Corresponds with errBounds_ */
    SingleVector* eigenFreqs;

    //! the recent calculated eigen values, might be complex
    SingleVector* eigenVals;

    Vector<Double> Frequency; // frequency
    Vector<Double> DampingRatio; // damping ratio
    Vector<Double> eigsRe; // real part
    Vector<Double> eigsIm; // imag part

    /** this is the list of wave vectors we have to process.
     * Obtained arbitrary */
    StdVector<Vector<double> > wave_vectors;

    StdVector<int> ModeOrder; // for sorting the obtained modes

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

    /** do we do ibz? Then bloch.plot will repeat the first step as final step. */
    bool ibz_;

    /** store the first plot.dat line to be repeated in the ibz_ case as last step */
    std::string first_plot_line_;

    /** the step number is complicated with bloch and or optimization. Count by ourselves here! */
    unsigned int save_step_;

    bool eigenValuesAreReal_;

  };

}

#endif
