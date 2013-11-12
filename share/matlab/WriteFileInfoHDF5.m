%WRITEFILEINFOHDF5 writes the FileInfo group of a CFS++ HDF5 file.
%
% WRITEFILEINFOHDF5(FILENAME,CONTENT) opens FILENAME and creates the group
% /FileInfo, in which four new dataset are created:
%
% - "Creator" contains the signature of this program (including user,
%   hostname, OS version, processor architecture and MATLAB version).
%
% - "Date" contains the current date and time.
%
% - "Version" is the version of the CFS++ HDF5 specification (= 0.9).
%      
% - "Content" is a 32-bit integer array filled with the CONTENT parameter.
%   It describes the contents of the file by a combination the following flags:
%     1: file contains a mesh
%     2: file contains volume results
%     4: file contains history results

% $Id$

function [] = WriteFileInfoHDF5(filename, content)

user = getenv('USER');
host = getenv('HOSTNAME');
os = getenv('OSTYPE');
if isempty(os)
  os = getenv('OS');
end
proc = getenv('HOSTTYPE');

creator = sprintf('CFS++ HDF5 tools for MATLAB run by %s@%s (MATLAB %s, %s %s)', ...
                  user, host, version, os, proc);
try
  h5datacreate(filename, '/FileInfo/Content', 'type', 'uint32', ...
               'size', length(content));
  h5varput(filename, '/FileInfo/Content', uint32(content));
catch %#ok<CTCH>
end

h5WriteVLStrDset(filename, '/FileInfo', 'Creator', creator, true);

h5WriteVLStrDset(filename, '/FileInfo', 'Date', datestr(now, 0), true);

h5WriteVLStrDset(filename, '/FileInfo', 'Version', '0.9', true);
