% ==============================================================
%
%    GENERAL
%    Writes the FileInfo group of an HDF5 file.
%
%    INPUT/S
%      filename  - full path of a new or existing HDF5 file
%      content   - integer vector indicating the content of the file
%
%    OUTPUT/S
%
%      
%    ABOUT
%
%      -Created:     25 Oct 2007
%      -Revision:    0.1
%      -Authors:     Jens Grabinger
% ==============================================================


function [] = WriteFileInfoHDF5(filename, content)

user = getenv('USER');
host = getenv('HOSTNAME');
os = getenv('OSTYPE');
if isempty(os)
  os = getenv('OS');
end
proc = getenv('HOSTTYPE');

creator = sprintf('CFS++ HDF5 tools for MATLAB run by %s@%s (MATLAB %s, %s %s)', user, host, version, os, proc);

if exist(filename, 'file') == 2
  hdf5write(filename, '/FileInfo/Creator', creator, 'WriteMode', 'append');
else
  hdf5write(filename, '/FileInfo/Creator', creator);
end

hdf5write(filename, '/FileInfo/Date', datestr(now, 0), 'WriteMode', 'append');

hdf5write(filename, '/FileInfo/Version', '0.9', 'WriteMode', 'append');

hdf5write(filename, '/FileInfo/Content', uint32(content), 'WriteMode', 'append');
