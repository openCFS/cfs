function h5dump(varargin)
%H5DUMP  print hdf5 metadata.
%
%   H5DUMP(FILENAME) prints the entire file's metadata.
%
%   H5DUMP(FILENAME,OBJ) prints the metadata for the specified
%   object.  If obj is a group, all objects below the group will
%   be described.
%
%   H5DUMP(...,'attr_mode',MODE) dumps the attributes according to
%   MODE, which can have two values.
%
%       'medium' - dump the attribute name, datatype, and dataspace
%       'lite'   - only print the attribute name
%     
%   This function requires R2007b or higher.

%   Copyright 2007-2010 The MathWorks, Inc.


global indent_char;
global dump_attr_mode;

indent_char = '   ';
dump_attr_mode = 'medium';

hargs = parse_args(varargin{:});

v = version('-release');
switch(v)
	case { '2006b', '2007a' }
		error('MATLAB:h5dump:requiresR2007bOrHigher', ...
		    'This function requires R2007b or higher to run.' );
	case { '2007b', '2008a', '2008b', '2009a' }
		h5dump_2009a(hargs{:});
	otherwise
		h5dump_2009b(hargs{:});
end


%--------------------------------------------------------------------------
function hargs = parse_args(varargin)
	
	global dump_attr_mode;
    
    switch nargin
        case {1, 3}
            hargs{1} = varargin{1};
            
        case 2
            hargs{1} = varargin{1};
            hargs{2} = varargin{2};
            
        otherwise
            hargs = varargin(1:2);
    end

    if nargin < 3
        return
    end
    
	j = 1;
	while true
		key = varargin{j};
		if ischar(key)

			if nargin == j
				error('HDF5TOOLS:h5dump:tooFewInputs', ...
				      'Must supply an argument for ''%s'' parameter.', key );
			end

			switch lower(key)
				case 'attr_mode'
					if ischar(varargin{j+1})
						value = lower(varargin{j+1});
						switch(value)
							case {'medium','lite'}
								dump_attr_mode = value;
							otherwise
								error('HDF5TOOLS:h5dump:unrecognizedParamValue', ...
								      'Unrecognized parameter value ''%s'' for ''%s''.', ...
									  value, key );
						end									  
					else
						error('HDF5TOOLS:h5dump:tooFewInputs', ...
						      'Must supply an argument for ''%s'' parameter.', key );
					end
					j = j+2;
                    
                otherwise
                    j = j+1;
			end

		else
			j = j+1;
		end

		if j > numel(varargin)
			break;
		end

	end
return
%--------------------------------------------------------------------------
function h5dump_2009a(fileName,obj)
% Use the group interface to iterate.
file = H5F.open (fileName, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

if nargin == 1
	obj = '/';
end

op_func(file,obj);

H5F.close (file);

return


%--------------------------------------------------------------------------
function status=op_func (loc_id, name)
% This is the master iteration function.  It just basically serves as a
% switchyard.  If the current object is a group, the iteration recurses
% there.

statbuf=H5G.get_objinfo (loc_id, name, 0);

switch (statbuf.type)
    case H5ML.get_constant_value('H5G_GROUP')
        dump_group(loc_id,name);
        
    case H5ML.get_constant_value('H5G_DATASET')
        dump_dataset(loc_id,name,1);
        
    case H5ML.get_constant_value('H5G_TYPE')
        dump_datatype(loc_id,name,1);
        
    case H5ML.get_constant_value('H5G_LINK')
        dump_link(loc_id,name,1);
        
    otherwise
        fprintf ( '  Unknown: %s\n', name);
end

status=0;

return


%--------------------------------------------------------------------------
function dump_link(loc_id,name, indent_level)

global indent_char;
indent = repmat(indent_char, 1,indent_level);

%statbuf=H5G.get_objinfo (loc_id, name, 0)
buf = H5G.get_linkval(loc_id,name);

fprintf(1,'%sSoft Link:  ''%s''==>''%s''\n',indent, name,buf);

%--------------------------------------------------------------------------
function dump_datatype(dtype,name, indent_level) 

global indent_char;
indent = repmat(indent_char, 1,indent_level);

%dclass = H5T.get_class(dtype);
%sz = H5T.get_size(dtype);

dtypeStr = get_datatype_string(dtype,indent_level);

if numel(name) > 0
    fprintf ('%sDatatype:  ''%s''  %s\n', indent, name, dtypeStr);
else
    fprintf ('%sDatatype:  %s\n', indent, dtypeStr);
end

return

%--------------------------------------------------------------------------
function dtypeStr = get_datatype_string(dtype,indent_level)


% Check vlen strings first.
if H5T.is_variable_str(dtype)
    dtypeStr = get_vlen_string(dtype,indent_level);
	return
end

class_id = H5T.get_class(dtype);

switch ( class_id )
    case H5ML.get_constant_value('H5T_BITFIELD')
        if H5T.equal(dtype,'H5T_STD_B8BE')
            dtypeStr = 'H5T_STD_B8BE';
        elseif H5T.equal(dtype,'H5T_STD_B8LE')
            dtypeStr = 'H5T_STD_B8LE';
        elseif H5T.equal(dtype,'H5T_STD_B16LE')
            dtypeStr = 'H5T_STD_B16LE';
        elseif H5T.equal(dtype,'H5T_STD_B16BE')
            dtypeStr = 'H5T_STD_B16BE';
        elseif H5T.equal(dtype,'H5T_STD_B32LE')
            dtypeStr = 'H5T_STD_B32LE';
        elseif H5T.equal(dtype,'H5T_STD_B32BE')
            dtypeStr = 'H5T_STD_B32BE';
        elseif H5T.equal(dtype,'H5T_STD_B64LE')
            dtypeStr = 'H5T_STD_B64LE';
        elseif H5T.equal(dtype,'H5T_STD_B64BE')
            dtypeStr = 'H5T_STD_B64BE';
        else
            error('HDF5TOOLS:h5dump:unrecognizedBitFieldType', ...
                  'Encountered an unrecognized bitfield type.' );
        end
    case H5ML.get_constant_value('H5T_INTEGER')
        if H5T.equal(dtype,'H5T_STD_I8BE')
            dtypeStr = 'H5T_STD_I8BE';
        elseif H5T.equal(dtype,'H5T_STD_I8LE')
            dtypeStr = 'H5T_STD_I8LE';
        elseif H5T.equal(dtype,'H5T_STD_U8BE')
            dtypeStr = 'H5T_STD_U8BE';
        elseif H5T.equal(dtype,'H5T_STD_U8LE')
            dtypeStr = 'H5T_STD_U8LE';
        elseif H5T.equal(dtype,'H5T_STD_I16BE')
            dtypeStr = 'H5T_STD_I16BE';
        elseif H5T.equal(dtype,'H5T_STD_I16LE')
            dtypeStr = 'H5T_STD_I16LE';
        elseif H5T.equal(dtype,'H5T_STD_U16BE')
            dtypeStr = 'H5T_STD_U16BE';
        elseif H5T.equal(dtype,'H5T_STD_U16LE')
            dtypeStr = 'H5T_STD_U16LE';
        elseif H5T.equal(dtype,'H5T_STD_I32BE')
            dtypeStr = 'H5T_STD_I32BE';
        elseif H5T.equal(dtype,'H5T_STD_I32LE')
            dtypeStr = 'H5T_STD_I32LE';
        elseif H5T.equal(dtype,'H5T_STD_U32BE')
            dtypeStr = 'H5T_STD_U32BE';
        elseif H5T.equal(dtype,'H5T_STD_U32LE')
            dtypeStr = 'H5T_STD_U32LE';
        elseif H5T.equal(dtype,'H5T_STD_I64BE')
            dtypeStr = 'H5T_STD_I64BE';
        elseif H5T.equal(dtype,'H5T_STD_I64LE')
            dtypeStr = 'H5T_STD_I64LE';
        elseif H5T.equal(dtype,'H5T_STD_U64BE')
            dtypeStr = 'H5T_STD_U64BE';
        elseif H5T.equal(dtype,'H5T_STD_U64LE')
            dtypeStr = 'H5T_STD_U64LE';
        else
            dtypeStr = sprintf('Unusual %d-bit integer', H5T.get_size(dtype)*8);
        end

    case H5ML.get_constant_value('H5T_FLOAT')
        if H5T.equal(dtype,'H5T_IEEE_F32BE')
            dtypeStr = 'H5T_IEEE_F32BE';
        elseif H5T.equal(dtype,'H5T_IEEE_F32LE')
            dtypeStr = 'H5T_IEEE_F32LE';
        elseif H5T.equal(dtype,'H5T_IEEE_F64BE')
            dtypeStr = 'H5T_IEEE_F64BE';
        elseif H5T.equal(dtype,'H5T_IEEE_F64LE')
            dtypeStr = 'H5T_IEEE_F64LE';
        else
            dtypeStr = sprintf('Unusual %d-bit floating point type', H5T.get_size(dtype)*8);
        end

    case H5ML.get_constant_value('H5T_OPAQUE')
        dtypeStr = get_opaque_string(dtype,indent_level);

    case H5ML.get_constant_value('H5T_REFERENCE')
        dtypeStr = 'H5T_REFERENCE';

    case H5ML.get_constant_value('H5T_COMPOUND')
        dtypeStr = get_compound_string(dtype,indent_level);


    case H5ML.get_constant_value('H5T_ENUM')
        dtypeStr = get_enum_string(dtype,indent_level);


    case H5ML.get_constant_value('H5T_ARRAY')
        dtypeStr = get_array_string(dtype,indent_level);

    case H5ML.get_constant_value('H5T_VLEN')
        dtypeStr = dump_vlen_datatype(dtype,indent_level+1);

    case H5ML.get_constant_value('H5T_STRING')
        dtypeStr = 'H5T_STRING';

    otherwise
        error('HDF5TOOLS:h5dump:unrecognizedClass', ...
              'Encountered an unrecognized datatype class %d.', ...
              class_id);
end


return

%--------------------------------------------------------------------------
function dtypeStr = get_opaque_string(dtype,indent_level)

global indent_char;
indent = repmat(indent_char, 1,indent_level);


tag = H5T.get_tag(dtype);
dtypeStr = sprintf('H5T_OPAQUE\n%sOpaque Tag: %s', indent, tag );

return


%--------------------------------------------------------------------------
function dtypeStr = get_compound_string(dtype,indent_level)

global indent_char;
indent = repmat(indent_char, 1,indent_level);
indent_p1 = repmat(indent_char, 1,indent_level+1);

dtypeStr = 'H5T_COMPOUND';

n = H5T.get_nmembers(dtype);
if n>0
    dtypeStr = sprintf('%s {\n', dtypeStr);
    for j = 0:n-1
        membName = H5T.get_member_name(dtype,j);
        membType = H5T.get_member_type(dtype,j);
        membDtypeStr = get_datatype_string(membType,indent_level+1);
        dtypeStr = sprintf('%s%s%s "%s"\n',dtypeStr,indent_p1,membDtypeStr,membName);
        H5T.close(membType);
    end
    dtypeStr = sprintf('%s%s}', dtypeStr, indent);
end

%--------------------------------------------------------------------------
function dtypeStr = get_array_string(dtype,indent_level)

global indent_char;
indent = repmat(indent_char, 1,indent_level);
indent_p1 = repmat(indent_char, 1,indent_level+1);

dtypeStr = 'H5T_ARRAY';

% What are the dimensions of the array?
[n,dims] = H5T.get_array_dims(dtype);
dims = fliplr(dims);
if (n>0)
    dimStr = num2str(dims(1));
    for j = 2:n
        dimStr = sprintf('%s x %d', dimStr, dims(j));
    end
    dtypeStr = sprintf('%s [%s]', dtypeStr,dimStr);
end

dtypeStr = sprintf('%s {\n', dtypeStr );

% What is the array datatype?
super_dtype = H5T.get_super(dtype);
super_dtype_str = get_datatype_string(super_dtype,indent_level+1);
H5T.close(super_dtype);
dtypeStr = sprintf('%s%s%s\n', dtypeStr, indent_p1, super_dtype_str);

dtypeStr = sprintf('%s%s}', dtypeStr, indent);

%--------------------------------------------------------------------------
function dtypeStr = get_enum_string(dtype,indent_level)

global indent_char;
indent = repmat(indent_char, 1,indent_level);
indent_p1 = repmat(indent_char, 1,indent_level+1);

dtypeStr = sprintf('%s {\n', 'H5T_ENUM' );
super_dtype = H5T.get_super(dtype);
super_dtype_str = get_datatype_string(super_dtype,indent_level+1);
H5T.close(super_dtype);
dtypeStr = sprintf('%s%s%s;\n', dtypeStr, indent_p1, super_dtype_str);

n = H5T.get_nmembers(dtype);
for j = 0:n-1
    name = H5T.get_member_name(dtype,j);       
    dtypeStr = sprintf('%s%s%-12s\t%4d\n', dtypeStr,indent_p1,name,j);
end

dtypeStr = sprintf('%s%s}', dtypeStr, indent);

return

%--------------------------------------------------------------------------
function dtypeStr = get_vlen_string(dtype,indent_level)

global indent_char;
indent = repmat(indent_char, 1,indent_level);
indent_p1 = repmat(indent_char, 1,indent_level+1);

    dtypeStr = sprintf('%s {\n', 'H5T_STRING' );
    dtypeStr = sprintf('%s%sStrSize H5T_VARIABLE;\n',dtypeStr,indent_p1);

    strPad = H5T.get_strpad(dtype);
	switch ( strPad )
		case H5ML.get_constant_value('H5T_STR_NULLTERM')
    		dtypeStr = sprintf('%s%sStrPad H5T_STR_NULLTERM;\n', dtypeStr,indent_p1);

		case H5ML.get_constant_value('H5T_STR_NULLPAD')
    		dtypeStr = sprintf('%s%sStrPad H5T_STR_NULLPAD;\n', dtypeStr,indent_p1);

		case H5ML.get_constant_value('H5T_STR_SPACEPAD')
    		dtypeStr = sprintf('%s%sStrPad H5T_STR_SPACEPAD;\n', dtypeStr,indent_p1);

		otherwise
			error('HDF5TOOLS:h5dump:unrecognizedStrPad', ...
			      'Unrecognized character string storage mechanism %d.', ...
				  strPad );
	end

    cset = H5T.get_cset(dtype);
    switch (cset)
        case H5ML.get_constant_value('H5T_CSET_ASCII')
            dtypeStr = sprintf('%s%sCset H5T_CSET_ASCII;\n', dtypeStr,indent_p1);
            
        otherwise
            error('HDF5TOOLS:h5dump:unrecognizedCset', ...
                'Unrecognized character set type %d.', ...
                cset);
    end
    

    dtypeStr = sprintf('%s%s}', dtypeStr, indent);

return


%--------------------------------------------------------------------------
function str = dump_vlen_datatype(dtype_id,indent_level)

super_dtype = H5T.get_super(dtype_id);
super_str = get_datatype_string(super_dtype, indent_level+1);
H5T.close(super_dtype);

str = sprintf ( 'H5T_VLEN { %s }', super_str );
return

%--------------------------------------------------------------------------
function dump_attribute(attr_id,indent_level)
global indent_char;
global dump_attr_mode;

indent = repmat(indent_char, 1,indent_level);
aname = H5A.get_name(attr_id);

fprintf('%sAttribute: ''%s'' ', indent, aname );
if strcmp(dump_attr_mode,'lite')
	fprintf('\n' );
	return
end

fprintf('{\n' );

dtype = H5A.get_type(attr_id);
dump_datatype(dtype,'',indent_level+1);
H5T.close(dtype);

space_id = H5A.get_space(attr_id);
dump_dataspace(space_id,indent_level+1);
H5S.close(space_id);

fprintf('%s}\n',indent);
return

%--------------------------------------------------------------------------
function dump_attributes(loc_id,indent_level)
% iterate over the group attributes
%H5A.iterate(loc_id,[],@attr_iterator);
total_attrs = H5A.get_num_attrs(loc_id);
for j = 0:total_attrs-1
    attr_id = H5A.open_idx(loc_id,j);
    dump_attribute(attr_id,indent_level);
    H5A.close(attr_id);
end
return

%--------------------------------------------------------------------------
function dump_group(loc_id,name)

if (name ~= '/')
    fullname = sprintf('%s/%s', H5I.get_name(loc_id), name);
else
	fullname = name;
end


if (numel(fullname)>1) && strcmp(fullname(1:2),'//')
    fullname = fullname(2:end);
end

fprintf ('\nGroup ''%s'' {\n', fullname);

gid = H5G.open(loc_id,name);
dump_attributes(gid,1);
H5G.iterate (loc_id, name,[] , @op_func);
H5G.close(gid);
fprintf ('} End Group ''%s''\n', name);
return

%--------------------------------------------------------------------------
function dump_dataspace(space_id,indent_level)

global indent_char;
indent = repmat(indent_char, 1,indent_level);

[ndims,dims,maxdims] = H5S.get_simple_extent_dims(space_id);
if ndims == 0
    fprintf ('%sDataspace scalar\n', indent);
else
    str = sprintf('%.0f ', flipud(dims(:)));
    str = deblank(str);

    strmax = '';
    for j = ndims:-1:1
        if (maxdims(j) == 18446744073709551616) ||(maxdims(j) == -1)  
            strmax = [strmax 'H5S_UNLIMITED ']; %#ok<AGROW>
        else
            strmax = [strmax num2str(maxdims(j)) ' ']; %#ok<AGROW>
        end
    end
    strmax = deblank(strmax);

    fprintf ('%sDataspace simple {[%s] / [%s]}\n', indent, str, strmax);
end

return

%--------------------------------------------------------------------------
function dump_layout(dset_id,indent_level)

global indent_char;
indent = repmat(indent_char, 1,indent_level);

dcpl = H5D.get_create_plist(dset_id);

layout = H5P.get_layout(dcpl);
switch ( layout )
    case H5ML.get_constant_value('H5D_COMPACT')
		layout_string = sprintf('Layout:  Compact');

    case H5ML.get_constant_value('H5D_CONTIGUOUS')
		layout_string = sprintf('Layout:  Contiguous');

    case H5ML.get_constant_value('H5D_CHUNKED')
		layout_string = sprintf('Layout:  Chunked');
		[rank,dims] = H5P.get_chunk(dcpl);
		str = '';
    	for j = rank:-1:1
	        str = [str num2str(dims(j)) ' ']; %#ok<AGROW>
	    end
		str = sprintf('{[%s]}', deblank(str));
		layout_string = sprintf ( '%s %s', layout_string, str);

end

fprintf ('%s%s\n', indent, layout_string);
H5P.close(dcpl);

return

%--------------------------------------------------------------------------
function dump_dataset(loc_id,name,indent_level)

global indent_char;
indent = repmat(indent_char, 1,indent_level);

fprintf ('%sDataset ''%s'' {\n', indent, name);

dset_id = H5D.open(loc_id,name);

dtype = H5D.get_type(dset_id);
dump_datatype(dtype,'',indent_level+1);
H5T.close(dtype);

space_id = H5D.get_space(dset_id);
dump_dataspace(space_id,indent_level+1);
H5S.close(space_id);

dump_layout(dset_id,indent_level+1);

dump_attributes(dset_id,2);


H5D.close(dset_id);

fprintf ('%s}\n', indent);
return




%--------------------------------------------------------------------------
function h5dump_2009b(filename,objName)
% Use the H5L interface to iterate.


file_id = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
if nargin == 1
	gid = H5G.open(file_id,'/');
	dump_group_2009b(gid);
	H5G.close(gid);
	H5F.close(file_id);
else
	oid = H5O.open(file_id,objName,'H5P_DEFAULT');
	info = H5O.get_info(oid);
	switch ( info.type )
		case H5ML.get_constant_value('H5O_TYPE_GROUP')
			dump_attributes(oid,1);
			dump_group_2009b(oid);
		case H5ML.get_constant_value('H5O_TYPE_DATASET')
			dump_dataset_2009b(oid);           
        case H5ML.get_constant_value('H5O_TYPE_NAMED_DATATYPE')
            dump_datatype(obj_id,name,1);
		otherwise
			fprintf('unknown object type %d', info.type );
	end
	H5O.close(oid);
end





%--------------------------------------------------------------------------
function dump_link_2009b(link_id,name,indent_level)
global indent_char;
indent = repmat(indent_char, 1,indent_level);

info = H5L.get_info(link_id,name,'H5P_DEFAULT');

switch(info.type)
    case H5ML.get_constant_value('H5L_TYPE_SOFT')
        v = H5L.get_val(link_id,name,'H5P_DEFAULT');
        fprintf('%sSoft Link:  ''%s'' ==> ''%s''\n', indent, name, v{1} );
            
    case H5ML.get_constant_value('H5L_TYPE_HARD')
        fprintf('%sHard Link:  %s\n', indent, name );
        
    case H5ML.get_constant_value('H5L_TYPE_EXTERNAL')
        fprintf('%sExternal Link:  %s\n', indent, name );
end

return

%--------------------------------------------------------------------------
function [status opdata_out]= H5Literator_func(gid,name,opdata)
opdata_out = opdata;
status = 0;
try
	obj_id = H5O.open(gid,name,'H5P_DEFAULT');
catch me
    switch (me.identifier)
            
        case {'MATLAB:hdf5lib2:H5G_traverse_real:failure', ...
                'MATLAB:hdf5lib2:H5Oopen:failure'}
            % Am I a link?  Should there not be a better way to catch this?
            dump_link_2009b(gid,name,2);
            return
            
        otherwise
            rethrow(me);
    end
end


info = H5O.get_info(obj_id);

switch ( info.type )
	case H5ML.get_constant_value('H5O_TYPE_GROUP')
		dump_group_2009b(obj_id);
	case H5ML.get_constant_value('H5O_TYPE_DATASET')
		dump_dataset_2009b(obj_id);
    case H5ML.get_constant_value('H5O_TYPE_NAMED_DATATYPE')
        dump_datatype(obj_id,name,1);
	otherwise
		fprintf('unknown object type %d', info.type );
end
H5O.close(obj_id);

status = 0;
return


%-----------------------------------------------------------------------------
function dump_group_2009b(gid)

gname = H5I.get_name(gid);
fprintf ('\nGroup ''%s'' {\n', gname);
dump_attributes(gid,1);
H5L.iterate(gid, 'H5_INDEX_NAME', 'H5_ITER_INC', 0, @H5Literator_func, []);
fprintf ('} End Group ''%s''\n', gname);

%-----------------------------------------------------------------------------
function dump_dataset_2009b(dset_id)

global indent_char;

dset_name = H5I.get_name(dset_id);
fprintf ('%sDataset ''%s'' {\n', indent_char, dset_name);

if H5DS.is_scale(dset_id)
    fprintf('%s%sDimension Scale\n', indent_char, indent_char);
end


dtype = H5D.get_type(dset_id);
dump_datatype(dtype,'',2);
H5T.close(dtype);

space_id = H5D.get_space(dset_id);
dump_dataspace(space_id,2);
H5S.close(space_id);

dump_layout(dset_id,2);
dump_attributes(dset_id,2);

fprintf ('%s} End Dataset ''%s''\n\n', indent_char, dset_name);



