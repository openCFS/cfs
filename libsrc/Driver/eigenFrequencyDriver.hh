#ifndef FILE_EIGENVALUE_DRIVER_HH
#define FILE_EIGENVALUE_DRIVER_HH

#include "singleDriver.hh"

namespace CoupledField {


  //! Driver class for calculating a general eigenvalue problem
  class EigenFrequencyDriver : public SingleDriver {

  public:

    //! constructor
    EigenFrequencyDriver( std::string driverTag = "anyTag",
                          bool isPartOfSequence = false );
    
    //! Destructor 
    ~EigenFrequencyDriver();

    //! Initialization method
    void Init();
  
    //! Main method solution method

    //! This method constitutes the actual driving method which controls the
    //! solution process for the problem.
    void SolveProblem();
    
    //! Return current time / frequency step of simulation
    UInt GetActStep( const std::string& pdename ) { return 1;}

  private:

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
