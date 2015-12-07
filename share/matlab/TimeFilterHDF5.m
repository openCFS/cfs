%    TimeFilterHDF5
%
%      Frequency filter tool for CFS++ HDF5 data.
%
%      This script reads transient results from HDF5 files and performs a FFT
%      on them. Then unwanted frequencies are set to zero and an inverse FFT
%      is applied to obtain the filtered signal in the time domain. Afterwards
%      the data is written to a (different) HDF5 file.
%
% Input Parameters
%   * infile   - path of input HDF5 file
%   * outfile  - path of output HDF5 file
%   * quantity - which quantity to filter
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


function [] =  TimeFilterHDF5(infile, outfile, quantity, region, lowfreq, highfreq, bufsize)

FftHdf5Core(2, infile, outfile, quantity, region, lowfreq, highfreq, bufsize);
