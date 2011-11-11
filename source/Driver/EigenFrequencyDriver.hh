#ifndef FILE_EIGENVALUE_DRIVER_HH
#define FILE_EIGENVALUE_DRIVER_HH

#include "SingleDriver.hh"

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
                           bool isPartOfSequence = false );
    
    //! Destructor 
    ~EigenFrequencyDriver();

    //! Initialization method
    void Init();
  
    //! Main method solution method

    //! This method constitutes the actual driving method which controls the
    //! solution process for the problem.
    void SolveProblem(bool write_results = true, PtrParamNode given_analysis_id = PtrParamNode(), AdjointParameters* adjointParams = NULL);
    
    //! Return current time / frequency step of simulation
    UInt GetActStep( const std::string& pdename ) { return 1;}

  private:
    /** This is the templated form to handle the general and quadratic case */
    template <class T>
    void PrintResult(SingleVector* frequencies, Vector<Double> bounds, 
                     ResultHandler* resHandler, UInt numConverged);
    
    //! Flag indicating, if a quadratic eigenvalue problem is to
    //! be solved
    bool isQuadratic_;
    
    //! Number of eigenfrequencies to be calculated
    UInt numFreq_;

    //! Shift for eigenvalues
    Double freqShift_;

    //! Flag for writing the eigenmods into the file
    bool writeModes_;
  };

}

#endif
