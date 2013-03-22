% ==============================================================
%
%    Time2FreqHDF5
%
%
%    GENERAL
%      Fourier transform for CFS++ HDF5 data.
%
%      This script reads transient results from HDF5 files and performs a
%      FFT on them. Afterwards the harmonic data is written to a (different)
%      HDF5 file.
%
%      ATTENTION: PATH must include h5tool for this script to work correctly.
%
%
%    INPUT/S
%      infile    - path of input HDF5 file
%      outfile   - path of output HDF5 file
%      quantity  - which quantity to convert
%      region    - region the quantity is defined on
%      lowfreq   - lowest frequency to be stored (0 for unlimited)
%      highfreq  - highest frequency to be stored (0 for unlimited)
%      bufsize   - maximum memory consumption (in megabytes)
%        
%    OUTPUT/S
%
%      
%    ABOUT
%
%      -Created:     Jan 2006
%      -Last update: 3 Jun 2008
%      -Revision:    0.6
%      -Authors:     Max Escobar, Simon Triebenbacher, Jens Grabinger
%
% =======================================================================


function [] =  Time2FreqHDF5(infile, outfile, quantity, region, lowfreq, highfreq, bufsize)

FftHdf5Core(1, infile, outfile, quantity, region, lowfreq, highfreq, bufsize);
