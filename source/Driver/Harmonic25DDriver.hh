// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
/*!
 *       \file     Harmonic25DDriver.hh
 *       \brief    Driver for 2.5D Harmonic analysis
 *
 *       \date     November 2024
 *       \author   Likun Luo
 */
//================================================================================================
#ifndef FILE_HARMONIC_25D_HH
#define FILE_HARMONIC_25D_HH

#include "SingleDriver.hh"
#include <memory>
#include "Utils/Timer.hh"
#include <fstream>

namespace CoupledField {

class Timer;

class Harmonic25DDriver : public virtual SingleDriver {

  public:

    //! Constructor
    //! \param sequenceStep current step in multisequence simulation
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    Harmonic25DDriver(UInt sequenceStep,
                  bool isPartOfSequence,
                  shared_ptr<SimState> state,
                  Domain *domain,
                  PtrParamNode paramNode,
                  PtrParamNode infoNode);
    
    //! Destructor
    virtual ~Harmonic25DDriver();

    //! Return current frequency step of simulation
    UInt GetActStep(const std::string& pdename) {
      return actFreqStep_;
    }

    /** Helper method which determines if an AnalyisType is complex. */
    virtual bool IsComplex() { return true; };

    //! Initialization method
    void Init(bool restart);

    //! Main method, which computes the wavenumber spectrum and performs fourier inverse transform
    void SolveProblem();

    //! \copydoc SingleDriver::SetToStepValue
    virtual void SetToStepValue(UInt stepNum, Double stepVal );

    /** This is the list of all frequencies. As long as we have no adaptive
    * frequeceny steps this makes no problem. */
    struct Frequency
    {
      /** the frequencys are 1-based such freqs[i].step should be i+1 */
      unsigned int step;

      double freq;

      friend std::ostream &
      operator<<(std::ostream & os, Frequency const & freqStp)
      {
        return os << "Frequency: " << freqStp.freq << ", Step: " 
          << freqStp.step;
      }
    };
    
    //! Inclusive Index range [lo,hi] into waveNum_ that is still unsolved
    struct Interval {UInt lo, hi;};

    // Structure for IFT evaluation position
    // {
    //   /** the frequencys are 1-based such freqs[i].step should be i+1 */
    //   unsigned int step;

    //   double pos;

    //   friend std::ostream &
    //   operator<<(std::ostream & os, Position const & posStp)
    //   {
    //     return os << "Position: " << posStp.pos << ", Step: " 
    //       << posStp.step;
    //   }
    // };


    //! Vector containing the wavenumbers
    StdVector<Frequency> waveNum_;

    //! Vector containing the iFT evalution positions
    // StdVector<Position> iftPos_;

    // compute single frequency step for one wavenumber
    // Double ComputeFrequencyStep(Frequency const& freqStp);

    // Solve one step by index; forceRefactor=true for an interval anchor.
    // Writes result under waveNum_[idx].step, sets solved_[idx], returns t_solve.
    Double ComputeFrequencyStep(UInt idx, bool forceRefactor);

    //! Center-out scheduler: solve one interval (anchor + left/right sweeps).
    void SolveInterval(UInt lo, UInt hi);
    //! Pick the anchor index for [lo,hi] (index midpoint for now).
    UInt SelectCenter(UInt lo, UInt hi) const { return (lo + hi) / 2; }


    //! Compute the whole wavenumber spectrum
    void CalcWavenumberSpectrum();

    //! Perform fourier inverse transform on the wavenumber spectrum to get the spatial solution
    // void CalcInverseTransform();
  
  protected:

    void ReadWaveNumbers();

    // bool ReadEvalPosList();

    //! Base frequency for which a simulation is performed
    Double baseFreq_;

    //! cutoff frequency of the wavenumber spectrum
    Double freqCutoff_;

    //! frequency resolution of the wavenumber spectrum
    Double freqResolution_;

    //! First frequency for which a simulation is performed
    Double startFreq_ = 0.0;

    //! Last frequency for which a simulation is performed
    Double stopFreq_;
  
    //! Last frequency step
    UInt stopFreqStep_;

    //! Number of frequencies for which a simulation is performed
    UInt numFreq_;

    //! Current frequency value
    Double actFreq_;
    
    //! Current frequency step
    UInt actFreqStep_;

    bool warmStart_ = true;
    bool reuseFactorization_ = false;
    UInt maxIterBeforeRefactorize_ = 10;
    bool adaptiveTime_ = false;
    Double refactorizeFraction_ = 0.8;
    Double lastSolveTime_ = 0.0;
    std::ofstream sweepLog_;                //!< per-step CSV
    Timer solveTimer_;
    //! Center-out sweep scheduler state.
    bool  centerOut_ = false;               //!< false = legacy monotone sweep
    bool reuseEnabled_ = false;
    bool firstFactorization_ = true;
    Double lastNumFact_ = 0.0;              //!< T_numfact of the current interval's anchor
    std::vector<Interval> pending_;         //!< LIFO of unsolved index ranges
    std::vector<bool>     solved_;          //!< coverage guard, sized numFreq_
    std::vector<UInt>     solveOrder_;      //!< visiting order, logging only

    //! Number of evaluation positions for IFT
    // UInt numPos_;

    //! Current position value
    // Double actPos_;

    //! Current position step
    // UInt actPosStep_;

    //! Flag for writing the wavenumber spectrum to results or not
    // bool storeSpectrum_;

};

}


#endif //FILE_HARMONIC_25D