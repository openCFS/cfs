// =====================================================================================
//
//       Filename:  fft.h
//
//    Description:  This is the header File for the implementation of FFT and IFFT
//                  based on the FFTW library. Desired features include:
//                  - perform FFT of transient hdf5 simulation results for arbitrary quantities
//                  - optionally apply windowing functions to improve the FFT output
//                  - optionally apply filter funtions (Low-Pass,High-Pass,Band-Pass)
//                  - perform IFFT
//                  - read in data in chunks to prevent memory overflow
//                  Code is to be created by ragula Srikanth Reddy for master thesis
//                  SrikanthReddy.Ragula@aau.at
//
//        Version:  1.0
//        Created:  05/25/2011 11:03:53 AM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
//
// =====================================================================================

#ifndef FILE_FFT_HH
#define FILE_FFT_HH

#define SIZE 1000

#include "fftw3.h"
#include "General/environment.hh"
#include "MatVec/vector.hh"
#include "MatVec/matrix.hh"
#include "Utils/result.hh"
#include "Domain/resultInfo.hh"
#include "General/defs.hh"
#include<fstream>
#include <vector>
#include "DataInOut/simInput.hh"
#include "DataInOut/simOutput.hh"
#include <complex>
#include <cmath>


namespace CoupledField {
  class FFT{
    public:
      //constructor

    typedef enum{

      NOWINDOW,
      HANNING,
      RECTANGULAR,
      BLACKMAN,
      BLACKMANHARRIS,
      NUTTALL,
      FLATTOP

    }WindowType;

    void FT(const Matrix<Double> & timeValues, Double Fs, Double fmin, Double fmax,
             std::vector< std::vector<Complex> > & freqValues, StdVector<Double>& freqSteps);

    void WIN(WindowType window,const Matrix<Double> & timeValues, Matrix<Double> & wtimeValues);

    void WINFFT(const Matrix<Double> & timeValues, const Matrix<Double> & wtimeValues, Double fmin, Double fmax,
            std::vector< std::vector<Complex> > & wfValues, const StdVector<Double>& freqSteps );

    void IFT(const Matrix<Double> & timeValues, const std::vector< std::vector<Complex> > & wfValues,
         std::vector< std::vector<Double> > & ifftValues, Double tstep, const StdVector<Double>& freqSteps, StdVector<Double>& timeSteps);




    //destructor
    ~FFT();

    protected:

    private:
      fftw_complex  *data, *wdata, *fft_result, *wfft_result, *fdata;
      double *ifft_result;
      fftw_plan plan_forward, plan_backward,wplan_forward;

      double fs;
      UInt xi;
      Vector<Double>winVal;

      StdVector<Double>fSteps;
      StdVector<Double>tSteps;

      Complex temp;


  };
}


#endif
