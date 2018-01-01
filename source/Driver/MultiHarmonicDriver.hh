// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     MultiHarmonicDriver.hh
 *       \brief    Driver for nonlinear frequency domain analysis
 *
 *       \date     Jan 2, 2017
 *       \author   kroppert
 */
//================================================================================================

#ifndef FILE_MULTIHARMONICDRIVER
#define FILE_MULTIHARMONICDRIVER

#include "SingleDriver.hh"

#include <boost/shared_ptr.hpp>


namespace CoupledField
{


class Timer;

//! driver for multiharmonic problems (nonlinear in frequency-domain). it is derived from BaseDriver
class MultiHarmonicDriver : public virtual SingleDriver
{

public:

  //! Constructor
  //! \param sequenceStep current step in multisequence simulation
  //! \param isPartOfSequence true, if driver is part of  multiSequence
  MultiHarmonicDriver(UInt sequenceStep, bool isPartOfSequence,
                 shared_ptr<SimState> state,
                 Domain* domain, PtrParamNode paramNode, PtrParamNode infoNode );


  //! Destructor
  virtual ~MultiHarmonicDriver();

  //! Return current time / frequency step of simulation
  UInt GetActStep( const std::string& pdename ) { return 1;}

  //! Initialization method
  void Init(bool restart);

  //! Main method, where harmonic analysis is implemented.
  void SolveProblem();

  /** Helper method which determines if an AnalyisType is complex. */
  virtual bool IsComplex() { return true; };


protected:

  //! Base frequency for which a simulation is performed
  Double baseFreq_;
  
  //! Number of considered harmonics for the solution quantity
  UInt numHarmonics_N_;
  
  //! Number of considered harmonics for the nonlinearity, e.g. nonlinear BH-curve in electromagnetics
  UInt numHarmonics_M_;

  //! Timer for estimating remaining runtime
  boost::shared_ptr<Timer> timer_;

};

}

#endif // FILE_MULTIHARMONICDRIVER
