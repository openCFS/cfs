// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_STATICDRIVER_2001
#define FILE_STATICDRIVER_2001

#include "singleDriver.hh"

namespace CoupledField {

  //! driver for static problems. it is derived from BaseDriver
  class StaticDriver : public SingleDriver {

  public:
    //! constructor
    //! \param sequenceStep current step in multisequence simulation
    //! \param true, if driver is part of  multiSequence
    StaticDriver( UInt sequenceStep,
                  bool isPartOfSequence = false );

    //! Destructor 
    ~StaticDriver();
  
    //! Return current time / frequency step of simulation
    UInt GetActStep( const std::string& pdename ) { return 1;}

    //! Initialization method
    void Init();

    /** Main method solution method. 
     * This method constitutes the actual driving method which controls the
     * solution process for the problem. The results are not implicitley written.
     * @param optimizationIteration to be called frequently by the optimizer */
    void SolveProblem(int optimizationIteration);
    
    void SolveProblem()
    {
        SolveProblem(1);
    }
    
    /** SolveProblem() does not write the results - used for optimization */ 
    void StoreResults(double step_val = -1.0);
    
  private:
    
    
  };

}

#endif // FILE_STATICDRIVER
