function h5id=h5filecreate(h5name,varargin)
%H5FILECREATE  Create HDF5 file.
%
%   h5id = H5CREATEFILE(HDFFILE,parameter, value,...) creates and file with 
%   the specified parameters.  The return value is the HDF5 file id.  The
%   file is closed if no return argument is specified.
%
%     'buffered'    - sets file access to use the standard I/O driver.
%     'buffer_size' - sets size of data buffer.  Can only be used with
%                     'buffered' option.  
%     'truncate'    - whether or not to create a few file (default is true)
%     'userblock'   - size of the userblock
%
%   This function requires R2009b or higher.
%
%   Example:
%       fid = h5filecreate('myfile.h5','userblock',512);
%       H5F.close(fid);
%
%   Credit where credit is due:  Philip Top at LLNL
%
%   See also  H5P.set_fapl_stdio, H5P.set_sieve_buf_size, 
%   H5P.set_userblock, h5datacreate.

%   Copyright 2010 The MathWorks, Inc.

p = inputParser;
p.addParamValue('buffer_size',0,@isnumeric);
p.addParamValue('buffered',false,@islogical);
p.addParamValue('truncate',true,@islogical);
p.addParamValue('userblock',[],@isnumeric);
p.parse(varargin{:});
params = p.Results;

if ~params.buffered && (params.buffer_size ~= 0)
    error('HDF5TOOLS:bufferSizeWithoutBuffered', ...
        'Cannot use ''buffer_size'' without also specifying ''buffered''.' );
end

fcpl=H5P.create('H5P_FILE_CREATE');
if (params.userblock>0)
    ubsize=512;
    while(ubsize<params.userblock)
        ubsize=ubsize*2;
    end
    H5P.set_userblock(fcpl,ubsize);
end


fapl=H5P.create('H5P_FILE_ACCESS');
%check for buffered mode
if (params.buffered)
    H5P.set_fapl_stdio(fapl);
    if (params.buffer_size>0)
        H5P.set_sieve_buf_size(fapl,params.buffer_size);
    end
end


if (params.truncate)
    mode='H5F_ACC_TRUNC';
else
    mode='H5F_ACC_EXCL';
end

h5id=H5F.create(h5name,mode,fcpl,fapl);
if nargout == 0
    H5F.close(h5id);
end
