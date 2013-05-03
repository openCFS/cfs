// =====================================================================================
//
//       Filename:  fft.cc
//
//    Description:  Implementation file for the FFT feature
//
//        Version:  1.0
//        Created:  05/25/2011 11:11:57 AM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
//
// =====================================================================================

#include "fft.hh"
#include "fftw3.h"
#include "General/environment.hh"
#include "MatVec/vector.hh"
#include "MatVec/matrix.hh"
#include "Utils/result.hh"
#include "Domain/resultInfo.hh"
#include "General/defs.hh"
#include "DataInOut/simInput.hh"
#include "DataInOut/simOutput.hh"

#define PI 3.14159265

#include<iostream>

#include <fstream>
#include<cmath>
#include<stdlib.h>
#include<new>
#include<iomanip>
#include <vector>
#include<complex>
//namespace boost{

namespace CoupledField {


 void FFT::FT(const Matrix<Double> & timeValues, Double Fs, Double fmin, Double fmax,
              std::vector< std::vector<Complex> > & freqValues, StdVector<Double>& freqSteps){

    UInt numNodes = timeValues.GetNumCols();
    UInt numTsteps = timeValues.GetNumRows();

    Double freqInterval = Fs/(numTsteps);
    UInt numFreqs = std::ceil((fmax +1)/freqInterval);
    UInt cutMinNumFreqs = std::ceil(fmin/freqInterval);
    Double fmin_first = cutMinNumFreqs * freqInterval;
    numFreqs -= cutMinNumFreqs;

    freqSteps.Resize(numFreqs);
    freqSteps.Init();
    fSteps.Resize(numFreqs);
    fSteps.Init();
    fSteps[0]=0.0;

    for(UInt actF=0;actF<numFreqs;actF++){

      fSteps[actF] = fmin_first + actF * freqInterval;
      freqSteps[actF]=fSteps[actF];
    }

    freqValues.resize(numNodes,std::vector<Complex>(numFreqs));


    data = (fftw_complex*) fftw_malloc (sizeof (fftw_complex)*(numTsteps));


    fft_result = (fftw_complex*) fftw_malloc (sizeof (fftw_complex)*(numTsteps));
    plan_forward = fftw_plan_dft_1d(numTsteps,data, fft_result, FFTW_FORWARD, FFTW_MEASURE );

 //ifft_result = (fftw_complex*) fftw_malloc (sizeof (fftw_complex)*(k));

 ///plan_backward = fftw_plan_dft_1d(numTsteps, fft_result, ifft_result, FFTW_BACKWARD, FFTW_MEASURE);

    for(UInt i=0; i<numNodes; i++){
      for(UInt j=0; j<numTsteps; j++){
        data[j][0] = timeValues[j][i];
        data[j][1] = 0;
      }

      fftw_execute(plan_forward );
      // fftw_execute(plan_backward);

      for(UInt j=0; j<numFreqs; j++){
      if(fSteps[j]>=fmin && fSteps[j]<=fmax){

         freqValues[i][j]= Complex(fft_result[j][0]*2/numFreqs,fft_result[j][1]*2/numFreqs);
         //std::cout<< freqValues[i][j]<<std::endl;
        }
      }
    }
    fftw_destroy_plan(plan_forward);
    fftw_free(data);
    fftw_free(fft_result);
  }




     void FFT::WIN(WindowType window,const Matrix<Double> & timeValues, Matrix<Double> & wtimeValues){

    UInt numNodes = timeValues.GetNumCols();
      UInt numTsteps = timeValues.GetNumRows();
      Double a0,a1,a2,alpha=0.16;
      a0=(1-alpha)/2;
      a1=0.5;
      a2=alpha/2;

      Double b0,b1,b2,b3;
      b0 = 0.35875;
      b1 = 0.48829;
      b2 = 0.14128;
      b3 = 0.01168;

      Double c0,c1,c2,c3;

      c0 = 0.355768;
      c1 = 0.487396;
      c2 = 0.141232;
      c3 = 0.012604;

      Double d0,d1,d2,d3,d4;

      d0 = 1;
      d1 = 1.93;
      d2 = 1.29;
      d3 = 0.388;
      d4 = 0.032;

      winVal.Resize(numTsteps);
      wtimeValues.Resize(numNodes,numTsteps);

      for(UInt i=0; i<numTsteps; i++){

        //if(i>=0 && i<= numTsteps){
         if(window == FFT::HANNING){
        xi=i;
        //if(i>=0 && i< std::ceil(numTsteps/4)){
          winVal[i]=(0.5*(1-std::cos(2*PI*xi/numTsteps)));


         }
         if(window==FFT::RECTANGULAR){
        winVal[i]=1;
         }
         if(window==FFT::BLACKMAN){
        xi=i;
        winVal[i]=(a0-a1*(std::cos(2*PI*xi/numTsteps))+a2*(std::cos(4*PI*xi/numTsteps)));
         }
         if(window==FFT::BLACKMANHARRIS){
          xi=i;
        winVal[i]= (b0-b1*(std::cos(2*PI*xi/numTsteps))+b2*(std::cos(4*PI*xi/numTsteps))-b3*(std::cos(6*PI*xi/numTsteps)));
         }
         if(window==FFT::NUTTALL){
          xi=i;
        winVal[i]= (c0-c1*(std::cos(2*PI*xi/numTsteps))+c2*(std::cos(4*PI*xi/numTsteps))-c3*(std::cos(6*PI*xi/numTsteps)));
         }
         if(window==FFT::FLATTOP){
          xi=i;
        winVal[i]= (d0-d1*(std::cos(2*PI*xi/numTsteps))+d2*(std::cos(4*PI*xi/numTsteps))-d3*(std::cos(6*PI*xi/numTsteps))+d4*(std::cos(8*PI*xi/numTsteps)));
         }
      }

        for(UInt i=0; i<numNodes; i++){
          for(UInt j=0; j<numTsteps; j++){
            wtimeValues[i][j]= timeValues[j][i]*winVal[j];

            //std::cout<< wdata[j][0]<<std::endl;
          }
        }

 }
 void FFT::WINFFT(const Matrix<Double> & timeValues, const Matrix<Double> & wtimeValues, Double fmin, Double fmax,
        std::vector< std::vector<Complex> > & wfValues, const StdVector<Double>& freqSteps ){

    UInt numNodes = timeValues.GetNumCols();
      UInt numTsteps = timeValues.GetNumRows();
      UInt numFreqs = freqSteps.GetSize();

      wfValues.resize(numNodes,std::vector<Complex>(numFreqs));

      wdata = (fftw_complex*) fftw_malloc (sizeof (fftw_complex)*(numTsteps));
      wfft_result = (fftw_complex*) fftw_malloc (sizeof (fftw_complex)*(numFreqs));
      wplan_forward = fftw_plan_dft_1d(numFreqs,wdata, wfft_result, FFTW_FORWARD, FFTW_MEASURE );


        for(UInt i=0; i<numNodes; i++){
          for(UInt j=0; j<numTsteps; j++){
            wdata[j][0] = wtimeValues[i][j];
            wdata[j][1] = 0;
            //std::cout<< wdata[j][0]<<std::endl;
          }


          fftw_execute(wplan_forward );

          for(UInt j=0; j<numFreqs; j++){
            if(freqSteps[j]>=fmin && freqSteps[j]<=fmax){
            wfValues[i][j]= Complex(wfft_result[j][0]*2/numFreqs,wfft_result[j][1]*2/numFreqs);
            }
          }
        }
        fftw_destroy_plan(wplan_forward);
        fftw_free(wdata);
        fftw_free(wfft_result);

 }





 void FFT::IFT(const Matrix<Double> & timeValues, const std::vector< std::vector<Complex> > & wfValues,
       std::vector< std::vector<Double> > & ifftValues, Double tstep, const StdVector<Double>& freqSteps , StdVector<Double>& timeSteps){

      UInt numNodes = timeValues.GetNumCols();
      UInt numTsteps = timeValues.GetNumRows();
     // std::cout<< numTsteps<<std::endl;
      UInt numFsteps = std::ceil(numTsteps);
      Double timeInterval = tstep;

      timeSteps.Resize(numTsteps);
      timeSteps.Init();

      tSteps.Resize(numTsteps);
      tSteps.Init();

      for(UInt actF=0;actF<numTsteps;actF++){
       tSteps[actF] = actF * timeInterval;
       timeSteps[actF]=tSteps[actF];

      }

      ifftValues.resize(numNodes,std::vector<Double>(numTsteps));
      fdata = (fftw_complex*) fftw_malloc (sizeof (fftw_complex)*(numTsteps));


      ifft_result = (double*) fftw_malloc (sizeof (double)*(numFsteps));


      plan_backward = fftw_plan_dft_c2r_1d(numTsteps, fdata, ifft_result, FFTW_MEASURE);
      for(UInt i=0; i<numNodes; i++){
        for(UInt j=0; j<freqSteps.GetSize(); j++){
        //if(tSteps[j]>0 && tSteps[j]<50){

         fdata[j][0] = (wfValues[i][j].real()*numFsteps/2);
         fdata[j][1] = (wfValues[i][j].imag()*numFsteps/2);

        }


        fftw_execute(plan_backward);

        for(UInt j=0; j<numTsteps; j++){

            timeSteps[j]=tSteps[j];

         ifftValues[i][j] = ifft_result[j]*1/numTsteps;

        }
      }

 }

  FFT::~FFT(){



  }
}

