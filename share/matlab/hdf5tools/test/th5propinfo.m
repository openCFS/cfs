function th5propinfo()

%   Copyright 2008-2009 The MathWorks, Inc.
warning ('off', 'HDF5TOOLS:h5propinfo:invalidPropName' );
%
test_datasetTransfer;    
%test_datasetCreate;  
%test_fileAccess; 
%test_fileCreate; 

switch ( version('-release') )
case { '2008b', '2008a', '2007b', '2007a', '2006b', '2006a' }
    % do nothing

otherwise
	test_v18;
end


fprintf ( 'All H5PROPINFO tests succeeded.\n\n' );

%--------------------------------------------------------------------------
function test_datasetTransfer()

switch(computer)
    case { 'MAC', 'MACI', 'GLNX86', 'PCWIN' }
        platform = 'thirtytwo_bit';
    case { 'MACI64', 'GLNXA64', 'PCWIN64', 'SOL64' }
        platform = 'sixtyfour_bit';
end


dxpl = H5P.create('H5P_DATASET_XFER');
actInfo = h5propinfo(dxpl);
H5P.close(dxpl);

load('th5propinfo.mat');

[majnum,minnum,relnum] = H5.get_libversion();
switch(minnum)
    case 6
        expInfo = expInfo.v16x.(platform).dxpl;
    case 8
        expInfo = expInfo.v18x.(platform).dxpl;
    otherwise
        error('unhandled hdf5 version');
end


if ~isequal(expInfo,actInfo)
    error('DXPL info struct did not match');
end


fprintf ( 1, 'Retrieval of H5P_DATASET_XFER property values succeeded .\n' );

%--------------------------------------------------------------------------
function test_datasetCreate()

switch(computer)
    case { 'MAC', 'MACI', 'GLNX86', 'PCWIN' }
        platform = 'thirtytwo_bit';
    case { 'MACI64', 'GLNXA64', 'PCWIN64', 'SOL64' }
        platform = 'sixtyfour_bit';
end


dcpl = H5P.create('H5P_DATASET_CREATE');
actInfo = h5propinfo(dcpl);
H5P.close(dcpl);

load('th5propinfo.mat');

[majnum,minnum,relnum] = H5.get_libversion();
switch(minnum)
    case 6
        expInfo = expInfo.v16x.(platform).dcpl;
    case 8
        expInfo = expInfo.v18x.(platform).dcpl;
    otherwise
        error('unhandled hdf5 version');
end

% pline and efl seem somewhat volatile, let's remove them
actInfo = rmfield(actInfo,{'efl','pline','fill_value'});
expInfo = rmfield(expInfo,{'efl','pline','fill_value'});


if ~isequal(expInfo,actInfo)
    error('DCPL info struct did not match');
end


fprintf ( 1, 'Retrieval of H5P_DATASET_CREATE property values succeeded .\n' );

%--------------------------------------------------------------------------
function test_v18()

load('th5propinfo.mat');
tmp = expInfo.v18x;

switch(computer)
    case { 'MAC', 'MACI', 'GLNX86', 'PCWIN' }
        platform = 'thirtytwo_bit';
    case { 'MACI64', 'GLNXA64', 'PCWIN64', 'SOL64' }
        platform = 'sixtyfour_bit';
end

h5_plist_name = {};                          ml_plist_name = {};
h5_plist_name{end+1} = 'H5P_DATASET_ACCESS';   ml_plist_name{end+1} = 'dapl';
h5_plist_name{end+1} = 'H5P_OBJECT_COPY';      ml_plist_name{end+1} = 'ocppl';
h5_plist_name{end+1} = 'H5P_OBJECT_CREATE';    ml_plist_name{end+1} = 'ocpl';
h5_plist_name{end+1} = 'H5P_FILE_MOUNT';       ml_plist_name{end+1} = 'fmpl';
h5_plist_name{end+1} = 'H5P_GROUP_CREATE';     ml_plist_name{end+1} = 'gcpl';
h5_plist_name{end+1} = 'H5P_GROUP_ACCESS';     ml_plist_name{end+1} = 'gapl';
h5_plist_name{end+1} = 'H5P_DATATYPE_CREATE';  ml_plist_name{end+1} = 'dtcpl';
h5_plist_name{end+1} = 'H5P_DATATYPE_ACCESS';  ml_plist_name{end+1} = 'dtapl';
h5_plist_name{end+1} = 'H5P_STRING_CREATE';    ml_plist_name{end+1} = 'scpl';
h5_plist_name{end+1} = 'H5P_ATTRIBUTE_CREATE'; ml_plist_name{end+1} = 'acpl';
h5_plist_name{end+1} = 'H5P_LINK_CREATE';      ml_plist_name{end+1} = 'lcpl';
h5_plist_name{end+1} = 'H5P_LINK_ACCESS';      ml_plist_name{end+1} = 'lapl';


for j = 1:numel(h5_plist_name)
	plist = H5P.create(h5_plist_name{j});
	actInfo = h5propinfo(plist);
	expInfo = tmp.(platform).(ml_plist_name{j});
	
	if ~isequal(expInfo,actInfo)
	    error('%s info struct did not match', h5_plist_name{j});
	end

	fprintf ( 1, 'Retrieval of %s property values succeeded .\n', h5_plist_name{j} );
end



%--------------------------------------------------------------------------
function test_fileAccess()

switch(computer)
    case { 'MAC', 'MACI', 'GLNX86', 'PCWIN' }
        platform = 'thirtytwo_bit';
    case { 'MACI64', 'GLNXA64', 'PCWIN64', 'SOL64' }
        platform = 'sixtyfour_bit';
end


fapl = H5P.create('H5P_FILE_ACCESS');
actInfo = h5propinfo(fapl);
H5P.close(fapl);

load('th5propinfo.mat');

[majnum,minnum,relnum] = H5.get_libversion();
switch(minnum)
    case 6
        expInfo = expInfo.v16x.(platform).fapl;
    case 8
        expInfo = expInfo.v18x.(platform).fapl;
    otherwise
        error('unhandled hdf5 version');
end

% mdc_initCacheCfg is volatile, remove it
actInfo = rmfield(actInfo,{'mdc_initCacheCfg'});
expInfo = rmfield(expInfo,{'mdc_initCacheCfg'});

if ~isequal(expInfo,actInfo)
    error('FAPL info struct did not match');
end

fprintf ( 1, 'Retrieval of H5P_FILE_ACCESS property values succeeded .\n' );
%--------------------------------------------------------------------------
function test_fileCreate()

switch(computer)
    case { 'MAC', 'MACI', 'GLNX86', 'PCWIN' }
        platform = 'thirtytwo_bit';
    case { 'MACI64', 'GLNXA64', 'PCWIN64', 'SOL64' }
        platform = 'sixtyfour_bit';
end


fcpl = H5P.create('H5P_FILE_CREATE');
actInfo = h5propinfo(fcpl);
H5P.close(fcpl);

load('th5propinfo.mat');

[majnum,minnum,relnum] = H5.get_libversion();
switch(minnum)
    case 6
        expInfo = expInfo.v16x.(platform).fcpl;
    case 8
        expInfo = expInfo.v18x.(platform).fcpl;
    otherwise
        error('unhandled hdf5 version');
end


if ~isequal(expInfo,actInfo)
    error('FCPL info struct did not match');
end

fprintf ( 1, 'Retrieval of H5P_FILE_CREATE property values succeeded .\n' );
