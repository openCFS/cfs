function pstruct = h5propinfo(property_list_id)
%H5PROPINFO:  return structure of information about a property list.
%   info = h5propinfo(plist_id) returns a structure of information about 
%   the property list identified by plist_id.
%
%   In most cases, the fieldnames are actual properties.  There is one
%   exception, though.  info also provides a CLASS field that specifies
%   the class name of the property list.  
%
%   EXAMPLE:
%       fapl = H5P.create('H5P_FILE_ACCESS');
%       fapl_info = h5propinfo(fapl);
%
%   Not all property names are valid MATLAB fieldnames.  These property 
%   names are changed by passing them thru GENVARNAMES.  The usual effect 
%   of this is to make the names camelCase with no spaces.
%
%   Property list types currently handled (to varying degrees) includes (in 
%   MATLAB versions R2008b and below):
%
%       dataset creation
%       dataset transfer
%       file access
%       file creation
%
%   MATLAB versions R2009a and higher also include support for
%
%       dataset access
%       object copy
%       object create
%       file mount
%       group access
%       group create
%       datatype access
%       datatype create
%       string create
%       attribute create
%       link access
%       link create
%

%   Copyright 2009-2010 The MathWorks, Inc.


propnames = cell(0);
H5P.iterate(property_list_id, 0, @iter_func );

struct_args = cell(2*numel(propnames),1);

for j = 1:numel(propnames)

    newpropname = genvarname(propnames{j});
	if ~strcmp(newpropname,propnames{j})
		warning('HDF5TOOLS:h5propinfo:invalidPropName', ...
		        'The property ''%s'' name is changed to ''%s'' to make it a legal MATLAB struct field name.',  ...
				propnames{j}, newpropname);
	end
    struct_args{2*(j-1)+1} = newpropname;

    % Some of the properties should be retrieved via specialized
    % functions rather than the generic "get".
    switch ( propnames{j} )

        case 'alloc_time'
            struct_args{2*j} = interpret_alloc_time(property_list_id);
            continue

        case 'btree_rank'
            struct_args{2*j} = interpret_btree_rank(property_list_id);
            continue

        case 'chunk_size'
            struct_args{2*j} = interpret_chunked_layout(property_list_id);
            continue

        case 'close_degree'
            struct_args{2*j} = interpret_close_degree(property_list_id);
            continue

        case 'fill_time'
            struct_args{2*j} = interpret_fill_time(property_list_id);
            continue

        case 'layout'
            struct_args{2*j} = interpret_layout(property_list_id);
            continue

        case 'symbol_leaf'
            struct_args{2*j} = interpret_symbol_leaf(property_list_id);
            continue

    end

    val = H5P.get(property_list_id,propnames{j});
    sz = H5P.get_size(property_list_id,propnames{j});


    % Handle any type-specific cases by name
    switch ( propnames{j} )
        case 'rdcc_w0'
            % We know that this one is double precision
            struct_args{2*j} = typecast(val','double');
            continue

    end

    % Try to intelligently typecast, guessing on total number of bytes.
    switch sz
        case { 4 }
            struct_args{2*j} = typecast(val','uint32');
        case { 8 }
            struct_args{2*j} = typecast(val','uint64');
        otherwise
            struct_args{2*j} = val';
    end
end


pstruct = struct(struct_args{:});

% Add the class
c = H5P.get_class(property_list_id);
pstruct.CLASS = H5P.get_class_name(c);


    function status = iter_func(id,name) %#ok<INUSL>
        % iteration callback function, just rolls up the property names.
        propnames{end+1} = name;
        status = 0;
    end

end

%-------------------------------------------------------------------------------
function alloc_time = interpret_alloc_time(property_list_id)
n = H5P.get_alloc_time(property_list_id);
switch (n )
    case H5ML.get_constant_value('H5D_ALLOC_TIME_DEFAULT')
        alloc_time = 'H5D_ALLOC_TIME_DEFAULT';
        
    case H5ML.get_constant_value('H5D_ALLOC_TIME_EARLY')
        alloc_time = 'H5D_ALLOC_TIME_EARLY';
        
    case H5ML.get_constant_value('H5D_ALLOC_TIME_INCR')
        alloc_time = 'H5D_ALLOC_TIME_INCR';
        
    case H5ML.get_constant_value('H5D_ALLOC_TIME_LATE')
        alloc_time = 'H5D_ALLOC_TIME_LATE';
        
    otherwise
        error('HDF5TOOLS:h5propinfo:badAllocTime', 'Unable to recognize alloc time property.');
end
end


%-------------------------------------------------------------------------------
function fill_time = interpret_fill_time(property_list_id)
n = H5P.get_fill_time(property_list_id);
switch (n )
    case H5ML.get_constant_value('H5D_FILL_TIME_IFSET')
        fill_time = 'H5D_FILL_TIME_IFSET';
        
    case H5ML.get_constant_value('H5D_FILL_TIME_ALLOC')
        fill_time = 'H5D_FILL_TIME_ALLOC';
        
    case H5ML.get_constant_value('H5D_FILL_TIME_NEVER')
        fill_time = 'H5D_FILL_TIME_NEVER';
        
    otherwise
        error('HDF5TOOLS:h5propinfo:badFillTime', 'Unable to recognize fill time property.');
end
end

%-------------------------------------------------------------------------------
function layout = interpret_layout(property_list_id)
n = H5P.get_layout(property_list_id);
switch (n )
    case H5ML.get_constant_value('H5D_COMPACT')
        layout = 'H5D_COMPACT';
        
    case H5ML.get_constant_value('H5D_CONTIGUOUS')
        layout = 'H5D_CONTIGUOUS';
        
    case H5ML.get_constant_value('H5D_CHUNKED')
        layout = 'H5D_CHUNKED';
        
    otherwise
        error('HDF5TOOLS:h5propinfo:badLayout', 'Unable to recognize layout property value.');
end
end



%-------------------------------------------------------------------------------
function close_degree = interpret_close_degree(property_list_id)
n = H5P.get_fclose_degree(property_list_id);
switch (n )
    case H5ML.get_constant_value('H5F_CLOSE_WEAK')
        close_degree = 'H5F_CLOSE_WEAK';
        
    case H5ML.get_constant_value('H5F_CLOSE_SEMI')
        close_degree = 'H5F_CLOSE_SEMI';
        
    case H5ML.get_constant_value('H5F_CLOSE_STRONG')
        close_degree = 'H5F_CLOSE_STRONG';
        
    case H5ML.get_constant_value('H5F_CLOSE_DEFAULT')
        close_degree = 'H5F_CLOSE_DEFAULT';
        
    otherwise
        error('HDF5TOOLS:h5propinfo:badCloseDegree', ...
              'Unable to recognize close degree property value.');
end
end



%-------------------------------------------------------------------------------
function btree_rank = interpret_btree_rank(plist_id)
btree_rank = H5P.get_sym_k(plist_id);
end

%-------------------------------------------------------------------------------
function symtab = interpret_symbol_leaf(plist_id)
[btree_rank, symtab] = H5P.get_sym_k(plist_id); %#ok<ASGLU>
end

%-------------------------------------------------------------------------------
function chunksize = interpret_chunked_layout(plist)

% Must check if we have a chunked layout first.
n = H5P.get_layout(plist);
if ( n == H5ML.get_constant_value('H5D_CHUNKED') )
	[rank,chunksize] = H5P.get_chunk(plist); %#ok<ASGLU>
else
	chunksize = [];
end

end
