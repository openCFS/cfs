function attval = h5attget(hfile,varname,attname)
%H5ATTGET  Retrieve attribute from HDF5 file.
% 
%   ATTVAL = H5ATTGET(HFILE,LOCATION,ATTNAME) retrieves an attribute named
%   ATTNAME from the HDF5 file HFILE.  LOCATION can be either a group
%   or a traditional dataset.  Use '/' for the root group.
%
%   Scalar string attributes are reported as a row string, rather than a
%   column string.  Other string attributes are returned as-is.
% 
%   1D attributes are returned as a single row, rather than as a column.
% 
%   Example:  
%       d = h5attget('example.h5','/g1/g1.1/dset1.1.1','att1');
%
%   See also h5attput.

%   Copyright 2009-2010 The MathWorks, Inc.

parent_is_dataset = true;

%
% Just use the defaults for now?
flags = 'H5F_ACC_RDONLY';
plist_id = 'H5P_DEFAULT';

file_id = H5F.open ( hfile, flags, plist_id );

if strcmp(varname,'/')
    parent_is_dataset = false;
    parent_id = H5G.open ( file_id, varname );
else
    try
        parent_id=H5D.open(file_id,varname);
    catch %#ok<CTCH>
        parent_is_dataset = false;
        parent_id = H5G.open ( file_id, varname );
    end
end

attr_id = H5A.open_name ( parent_id, attname );

attval = H5A.read ( attr_id, 'H5ML_DEFAULT' );

% Read the datatype information and use that to possibly post-process
% the attribute data.
attr_datatype = H5A.get_type(attr_id);
attr_class = H5T.get_class(attr_datatype);

attr_space_id = H5A.get_space(attr_id);
ndims = H5S.get_simple_extent_dims(attr_space_id);

% Do we alter the size, shape of the attribute?
switch (attr_class)

	case H5ML.get_constant_value('H5T_STRING')
		% Is it a scalar attribute?
		if ( ndims == 0 )
			% Yep, it's a scalar.  Transpose it.
			attval = attval';
		end

	case { H5ML.get_constant_value('H5T_INTEGER'), ...
	       H5ML.get_constant_value('H5T_FLOAT') }
		if ( ndims == 1 )
			% If 1D, make it a row vector.
			attval = reshape(attval,1,numel(attval));
		end
end

H5T.close(attr_datatype);


H5A.close(attr_id);

if parent_is_dataset
    H5D.close(parent_id);    
else
    H5G.close(parent_id);    
end
H5F.close(file_id);
