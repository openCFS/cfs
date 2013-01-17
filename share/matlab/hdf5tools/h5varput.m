function h5varput ( varargin )
%H5VARPUT  Write HDF5 dataset.
%   H5VARPUT(HDFFILE,VARNAME,DATA) writes an entire dataset to the variable 
%   given by VARNAME.
%
%   H5VARPUT(HDF5FILE,VARNAME,START,COUNT,DATA) writes a contiguous portion 
%   of a dataset.
%
%   H5VARPUT(HDF5FILE,VARNAME,START,COUNT,STRIDE,DATA) writes a 
%   non-contiguous portion of a dataset.
%
%   HDF5FILE may also be an identifier for an already open HDF5 file.
%
%   If the dataset is extendible, and if START, COUNT, and STRIDE are 
%   specified, and if the resulting dataspace selection makes sense, then
%   the dataset will be properly extended.
%
%   See also h5varget.

%   Copyright 2009-2010 The MathWorks, Inc.

error(nargchk(3,6,nargin,'struct'));

[h5file,varname,offset,count,stride,data] = parse_h5_varput_options ( varargin{:} );

flags = 'H5F_ACC_RDWR';
plist_id = 'H5P_DEFAULT';

if isa(h5file,'H5ML.id')
	file_id = h5file;
	close_file = false;
else
	file_id     = H5F.open(h5file,flags,plist_id);
	close_file = true;
end
if isa(varname,'H5ML.id')
    datasetId = varname;
    close_dset = false;
else
    datasetId  = H5D.open(file_id,varname);
    close_dset = true;
end

write_dataset(datasetId,data,offset,stride,count);

if close_dset
    H5D.close(datasetId);
end
if close_file
	H5F.close(file_id);
end


%--------------------------------------------------------------------------
function write_complex_dataset(datasetId,data,offset,stride,count)

% ok, read the references
dref = H5D.read(datasetId,'H5T_STD_REF_OBJ','H5S_ALL','H5S_ALL','H5P_DEFAULT');

% The first dereference should be to the real part.
dsetR = H5R.dereference(datasetId,'H5R_OBJECT',dref(:,1));
dsetI = H5R.dereference(datasetId,'H5R_OBJECT',dref(:,2));

% Write the real and complex parts.
write_dataset(dsetR,real(data),offset,stride,count);
write_dataset(dsetI,imag(data),offset,stride,count);




%--------------------------------------------------------------------------
function write_dataset(datasetId,data,offset,stride,block)

datatype_id = H5D.get_type(datasetId);

dcpl = H5D.get_create_plist(datasetId);
dxpl = 'H5P_DEFAULT';

% Try to honor a fill value if the datatype is floating point.
if H5T.equal(datatype_id,'H5T_NATIVE_DOUBLE') || H5T.equal(datatype_id,'H5T_NATIVE_FLOAT')
	% Turn any fill value into NaN
	if (H5P.fill_value_defined(dcpl) == H5ML.get_constant_value('H5D_FILL_VALUE_USER_DEFINED'))
		fillvalue    = H5P.get_fill_value(dcpl,datatype_id);
		data(isnan(data)) = fillvalue;
	end
end

fileSpaceId = getFileSpace(datasetId,offset,stride,block);

% Create the appropriate output dataspace.
if isempty(offset) && isempty(block) && isempty(stride)
	mem_space_id = 'H5S_ALL';
else
	mem_space_id = H5S.create_simple(length(offset),fliplr(block),fliplr(block));
end

H5D.write(datasetId,'H5ML_DEFAULT',mem_space_id,fileSpaceId,dxpl,data);

H5T.close(datatype_id);






%--------------------------------------------------------------------------
function fileSpaceId = getFileSpace(datasetId,offset,stride,count)

if isempty(offset) && isempty(count) && isempty(stride)
	% We write to the entire file space as it exists.  We don't allow
	% the user to extend in this case.
	fileSpaceId = 'H5S_ALL';
	return;
end

% Extents were specified.  Make the selections on the hyperslab.
fileSpaceId = H5D.get_space(datasetId);
H5S.select_hyperslab(fileSpaceId, 'H5S_SELECT_SET', ...
	                     fliplr(offset), fliplr(stride), ...
						 ones(1,length(offset)), ...
						 fliplr(count) );

% If any extents have been specified, then check to see if we extend the
% dataset.
[nSpaceDims, spaceDims] = H5S.get_simple_extent_dims(fileSpaceId); %#ok<ASGLU>
spaceDims = spaceDims(:);
[boundsStart, boundsEnd] = H5S.get_select_bounds(fileSpaceId); %#ok<ASGLU>
boundsEnd = boundsEnd(:);

% Do any of the bounding box dimensions exceed the current dimensions
% of the dataset?
%
% R2009a and earlier return boundsEnd as a column vector, but also returned
% spaceDims as a row vector.  Have to normalize for that.
idx = find(boundsEnd > (spaceDims-1), 1);
if isempty(idx)
	% Nope, no need to extend.  
	return;
end

new_dims = max(spaceDims,boundsEnd+1);

% Ok, at least one specified extent is beyond that of the current extents.
% Extend the dataset and return.
H5S.close(fileSpaceId);

v = version('-release');
switch(v)
	case { '2006b', '2007a' , '2007b', '2008a', '2008b' }
		% Hopefully the dataset won't be shrinking :-)
		H5D.extend(datasetId,new_dims);

	otherwise
		% In 9a and beyond, this function should be used.
		H5D.set_extent(datasetId,new_dims);

end
fileSpaceId = H5D.get_space(datasetId);
H5S.select_hyperslab(fileSpaceId, 'H5S_SELECT_SET', ...
	                     fliplr(offset), fliplr(stride), fliplr(count), ...
						 ones(1,length(offset)));


return




%===============================================================================
function [hfile,varname,start,count,stride,data] = parse_h5_varput_options ( varargin )
%
% Have to be able to check the following signatures.

% H5_VARGET(HFILE,VARIABLE,DATA) 
% H5_VARGET(HFILE,VARIABLE,START,COUNT,DATA)
% H5_VARGET(HFILE,VARIABLE,START,COUNT,STRIDE,DATA)

% First argument should be the filename.
if ~(ischar(varargin{1}) || isa(varargin{1},'H5ML.id'))
	error ( 'MATLAB:H5TOOLS:badInput', 'File argument must be a filename or an id.')
end

hfile = varargin{1};

%
% 2nd argument should be the variable name.
if ~(ischar(varargin{2}) || isa(varargin{2},'H5ML.id'))
	error ( 'MATLAB:H5VARPUT:badInput', ...
        'Variable name argument must be character or an id.')
end

varname = varargin{2};


switch nargin
case 3

	start = [];
	stride = [];
	count = [];
	data = varargin{3};

case 4
	error ( 'MATLAB:H5VARPUT:badInput', 'Cannot have 4 input arguments.')

case 5
	%
	% Given the start, stride, and count.
	if ~isnumeric(varargin{3}) 
		error ( 'MATLAB:H5VARPUT:badInput', 'Start argument must be numeric')
	end
	start = varargin{3};

	if ~isnumeric(varargin{4}) 
		error ( 'MATLAB:H5VARPUT:badInput', 'Count argument must be numeric')
	end
	count = varargin{4};

	stride = ones(size(count));
	data = varargin{5};

	
case 6

	%
	% Given the start, stride, and count.
	if ~isnumeric(varargin{3}) 
		error ( 'MATLAB:H5VARPUT:badInput', 'Start argument must be numeric')
	end
	start = varargin{3};

	if ~isnumeric(varargin{4}) 
		error ( 'MATLAB:H5VARPUT:badInput', 'Count argument must be numeric')
	end
	count = varargin{4};

	if ~isnumeric(varargin{5}) 
		error ( 'MATLAB:H5VARPUT:badInput', 'Stride argument must be numeric')
	end
	stride = varargin{5};
	data = varargin{6};

otherwise
	error ( 'MATLAB:H5VARPUT:badInput', 'Bad number of input arguments.')

end

return






