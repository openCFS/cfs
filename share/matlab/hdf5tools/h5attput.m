function h5attput(h5file,varname,attname,attvalue)
%H5ATTPUT  Write HDF5 attribute.
%
%   H5ATTPUT(H5FILE,VARNAME,ATTNAME,ATTVALUE) creates/overwrites the 
%   attribute named ATTNAME with the value ATTVALUE.  The parent object 
%   VARNAME can be either an HDF5 variable or group.  VARNAME must be a 
%   complete pathname.
%
%   Simple strings will be created in a scalar dataspace.
%
%   If H5FILE does not exist, it will be created.  In this case, the parent 
%   object must be '/' , (the root group)
%
%   See also h5attget.

%   Copyright 2008-2010 The MathWorks, Inc.

error(nargchk(4,4,nargin,'struct'));

if ~ischar(h5file)
	error('MATLAB:H5ATTPUT:badDatatype', ...
	      'Filename input argument must have datatype char.' );
end

if ~ischar(varname)
	error('MATLAB:H5ATTPUT:badDatatype', ...
	      'VARNAME input argument must have datatype char.' );
end

if ~ischar(attname)
	error('MATLAB:H5ATTPUT:badDatatype', ...
	      'ATTNAME input argument must have datatype char.' );
end


flags = 'H5F_ACC_RDWR';
plist_id = 'H5P_DEFAULT';

if ~exist(h5file,'file')
    create_plist = H5P.create('H5P_FILE_CREATE');
    file_id = H5F.create(h5file,'H5F_ACC_TRUNC', create_plist, 'H5P_DEFAULT');
    H5P.close(create_plist);
else
    file_id = H5F.open ( h5file, flags, plist_id );
end

[parent_id, parent_obj_is_group] = set_parent_id ( file_id, varname );
dataspace_id                     = set_dataspace_id ( attvalue );
datatype_id                      = set_datatype_id ( attvalue );
att_id                           = set_attribute_id ( parent_id, attname, datatype_id, dataspace_id );
H5A.write(att_id,datatype_id,attvalue);

H5T.close(datatype_id);
H5S.close(dataspace_id);
H5A.close(att_id);

if parent_obj_is_group
	H5G.close(parent_id);
else
	H5D.close(parent_id);
end

H5F.close(file_id);





%===============================================================================
% SET_ATTRIBUTE_ID
%
% Setup the attribute ID.  We need to check as to whether or not the attribute
% already exists.
function att_id = set_attribute_id ( parent_id, attname, datatype_id, dataspace_id )

try
	att_id = H5A.open_name ( parent_id, attname );
catch %#ok<CTCH>
	att_id = H5A.create ( parent_id, attname, datatype_id, dataspace_id,'H5P_DEFAULT' );
end

%===============================================================================
% SET_DATASPACE_ID
%
% Setup the dataspace ID.  This just depends on how many elements the 
% attribute actually has.
function dataspace_id = set_dataspace_id ( attvalue )

if ischar(attvalue)
    if ( ndims(attvalue) == 2 ) && ( any(size(attvalue) ==1) )
        dataspace_id = H5S.create('H5S_SCALAR');
        return
    elseif ndims(attvalue < 3)
        rank = 1;
        dims = size(attvalue,2);
    else
        error('HDF5TOOLS:h5attput:badStringSize', ...
            'Cannot accept a string input with ndims > 2.');
    end
else
    if ( ndims(attvalue) == 2 ) && ( any(size(attvalue) ==1) )
        rank = 1;
        dims = numel(attvalue);
    else
        % attribute is a "real" 2D value.		
        rank = ndims(attvalue);
	    dims = fliplr(size(attvalue));
    end
end
dataspace_id = H5S.create_simple ( rank, dims, dims );




%===============================================================================
% SET_PARENT_ID
%
% If the given variable is "/", then we know we are creating a group attribute.
% Otherwise try to open the variable as a dataset.  If that fails, then it
% must be a subgroup.
function [parent_id, parent_obj_is_group] = set_parent_id ( file_id, varname )
if strcmp(varname,'/')
    parent_obj_is_group = true;
    parent_id = H5G.open ( file_id, varname );
else
    try
        parent_id=H5D.open(file_id,varname);
    	parent_obj_is_group = false;
    catch %#ok<CTCH>
        parent_id = H5G.open ( file_id, varname );
    	parent_obj_is_group = true;
    end
end

%===============================================================================
% SET_DATATYPE_ID
%
% We need to choose an appropriate HDF5 datatype based upon the attribute
% data.
function datatype_id = set_datatype_id ( attvalue )
switch class(attvalue)
	case 'double'
	    datatype_id = H5T.copy('H5T_NATIVE_DOUBLE');
	case 'single'
	    datatype_id = H5T.copy('H5T_NATIVE_FLOAT');
	case 'int64'
	    datatype_id = H5T.copy('H5T_NATIVE_LLONG');
	case 'uint64'
	    datatype_id = H5T.copy('H5T_NATIVE_ULLONG');
	case 'int32'
	    datatype_id = H5T.copy('H5T_NATIVE_INT');
	case 'uint32'
	    datatype_id = H5T.copy('H5T_NATIVE_UINT');
	case 'int16'
	    datatype_id = H5T.copy('H5T_NATIVE_SHORT');
	case 'uint16'
	    datatype_id = H5T.copy('H5T_NATIVE_USHORT');
	case 'int8'
	    datatype_id = H5T.copy('H5T_NATIVE_SCHAR');
	case 'uint8'
	    datatype_id = H5T.copy('H5T_NATIVE_UCHAR');
	case 'char'
	    datatype_id = H5T.copy('H5T_C_S1');
        if any(size(attvalue) == 1)
            H5T.set_size(datatype_id,numel(attvalue));
        else
            H5T.set_size(datatype_id,size(attvalue,1));
        end
		H5T.set_strpad(datatype_id,'H5T_STR_NULLTERM');
	otherwise 
		error('MATLAB:H5ATTPUT:unsupportedDatatype', ...
		       '''%s'' is not a supported H5ATTPUT datatype.\n', class(attvalue) );
end
return

