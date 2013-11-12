%    Time2FreqHDF5
%
%      Fourier transform for CFS++ HDF5 data.
%
%      This script reads transient results from HDF5 files and performs a
%      FFT on them. Afterwards the harmonic data is written to a (different)
%      HDF5 file.
%
% Input Parameters
%   * infile   - path of input HDF5 file
%   * outfile  - path of output HDF5 file
%   * quantity - which quantity to convert
%   * region   - region the quantity is defined on
%   * lowfreq  - lowest frequency to be stored (0 for unlimited)
%   * highfreq - highest frequency to be stored (0 for unlimited)
%   * bufsize  - maximum memory consumption (in megabytes)
%
% Return Value
%   None
%
% About
%   * Created:  Jan 2006
%   * Authors:  Max Escobar, Simon Triebenbacher, Jens Grabinger
%   * Revision: $Id$


function [] =  Time2FreqHDF5(infile, outfile, quantity, region, lowfreq, highfreq, bufsize)

FftHdf5Core(1, infile, outfile, quantity, region, lowfreq, highfreq, bufsize);
