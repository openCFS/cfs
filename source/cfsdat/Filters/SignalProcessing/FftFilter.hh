/*
 * FftFilter.hh
 *
 *  Created on: Dec 5, 2016
 *      Author: jens
 */

#ifndef SOURCE_CAATOOL_FILTERS_SIGNALPROCESSING_FFTFILTER_HH_
#define SOURCE_CAATOOL_FILTERS_SIGNALPROCESSING_FFTFILTER_HH_

#include <cfsdat/Filters/BaseFilter.hh>
#include <General/Enum.hh>
#include "MatVec/Matrix.hh"

namespace CFSDat {

class FftFilter : public BaseFilter {
    
  public:
    
    //! Enum of available window functions
    typedef enum {
      WIN_RECT = 0,
      WIN_BLACKMAN,
      WIN_BLACKMAN_HARRIS,
      WIN_HAMMING,
      WIN_HANN,
      WIN_NUTTAL
    } WindowType;
    
    //! Auxiliary object for enum WindowType 
    static Enum<WindowType> WindowTypeEnum;
    
    //! Compute a window function at discrete points
    static void GetWindow(Vector<Double> & winVec, WindowType winType,
                          UInt numSamples, Double param = 0.0);
    
    //! Initialize the static Enum objects
    static CF::Enum<FftFilter::WindowType> InitWindowTypeEnum();
    
    //! Constructor
    FftFilter(UInt numWorkers, PtrParamNode config, PtrResultManager resMan);
    
    //! Destructor
    virtual ~FftFilter() {};
    
    //! Run the filter operation
    virtual bool Run();
    
  protected:
    
    //! Set up upstream results
    virtual ResultIdList SetUpstreamResults();
    
    //! Set up this filter's own results
    virtual void AdaptFilterResults();
    
    
    //! Number of samples per FFT
    UInt nFFT_;
    
    //! Total number of samples in time domain
    UInt numSamples_;
    
    //! Number of frequencies in result
    UInt numFreqs_;
    
    //! Type of window function
    WindowType windowType_;
    
    //! Window parameter (some window functions are parametric)
    Double windowParam_;
    
    //! Length of window (excluding zero padding)
    UInt windowLength_;
    
    //! Overlap of signal blocks
    Double overlap_;
    
    //! Time step size
    Double deltaT_;
    
    //! Frequency step size
    Double deltaF_;
    
    //! Buffer for all input data of the FFT (one time series per row)
    //! Will be overwritten with the result of the FFT.
    Matrix<Double> data_;
    
    //! Iterator pointing to the current frequency step
    std::map<UInt, Double>::iterator myStepIter_;
    
};

} /* namespace CAATool */

#endif /* SOURCE_CAATOOL_FILTERS_SIGNALPROCESSING_FFTFILTER_HH_ */
