#ifndef FILE_EIGENVALUE_DRIVER_HH
#define FILE_EIGENVALUE_DRIVER_HH

#include "singleDriver.hh"

namespace CoupledField {


  //! Driver class for calculating a general eigenvalue problem
  class EigenFrequencyDriver : public SingleDriver {

  public:

    //! constructor
    EigenFrequencyDriver(Domain * adomain,
                 UInt stepOffset = 0,
                 Double timeOffset = 0.0,
                 std::string driverTag = "anyTag",
                 Boolean isPartOfSequence = FALSE);
    
    //! Destructor 
    ~EigenFrequencyDriver();
  
    //! Main method solution method

    //! This method constitutes the actual driving method which controls the
    //! solution process for the problem.
    void SolveProblem();
    
  private:

    //! Number of eigenfrequencies to be calculated
    UInt numFreq_;

    //! Shift for eigenvalues
    Double freqShift_;

    //! Flag for using shift-mode of the eigenvalue solver
    Boolean shiftMode_;

    //! Flag for writing the eigenmods into the file
    Boolean writeModes_;
  };

}

#endif
