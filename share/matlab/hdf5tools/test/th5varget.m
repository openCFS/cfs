function th5varget()

%   Copyright 2008-2010 The MathWorks, Inc.
fprintf('\tTesting H5VARGET...  ');
setup_th5varget;

% Retrieve full datasets.
%
%    single precision
%    vlen(int32)
%    vlen(string)
%
test_retrieveFullSinglePrecisionDataset; 
test_retrieveFullVlenInt32Dataset; 
test_retrieveFullVlenStringDataset; 

% Retrieve contiguous dataset.
%
%    compound
%    single precision
%    vlen(int32)
%    vlen(string)
%
test_retrievePartialContiguousCompoundDataset; 
test_retrievePartialContiguousSinglePrecisionDataset; 
test_retrievePartialVlenInt32Dataset; 
test_retrievePartialVlenStringDataset; 

% Retrieve strided dataset.
%
%    single precision
%    vlen(string)
%
test_retrievePartialStridedSinglePrecisionDataset; 
test_retrieveStridedVlenStringDataset; 


test_openFileHandle; 

%test_retrieveFullEnumDataset; 
%test_retrievePartialCharEnumDataset; 

test_fillValueWithInt32Dataset; 
test_fillValueWithDoubleDataset; 

test_retrieveTooMuchData;

fprintf('OK\n');

%--------------------------------------------------------------------------
function setup_th5varget()
setpref('HDF5TOOLS','ENUMERATE',false);

%--------------------------------------------------------------------------
%function test_retrievePartialCharEnumDataset()
%
% In the example file, we have 
%
% HDF5 "data/h5ex_t_enum.h5" {
%  GROUP "/" {
%     DATASET "DS1" {
%        DATATYPE  H5T_ENUM {
%              H5T_STD_I16BE;
%              "SOLID"            0;
%              "LIQUID"           1;
%              "GAS"              2;
%              "PLASMA"           3;
%           }
%        DATASPACE  SIMPLE { ( 4, 7 ) / ( 4, 7 ) }
%        DATA {
%        (0,0): SOLID, SOLID, SOLID, SOLID, SOLID, SOLID, SOLID,
%        (1,0): SOLID, LIQUID, GAS, PLASMA, SOLID, LIQUID, GAS,
%        (2,0): SOLID, GAS, SOLID, GAS, SOLID, GAS, SOLID,
%        (3,0): SOLID, PLASMA, GAS, LIQUID, SOLID, PLASMA, GAS
%        }
%     }
%  }
%  }


%actData = h5varget('data/h5ex_t_enum.h5','/DS1',[0 3], [4 1],'enumerate',true);
%
%load('th5varget.mat');
%if ~isequal(actData,expData.partialCharEnum)
%	error ( 'failed to verify full enum dataset retrieval');
%end

%fprintf ( 1, 'Retrieval of a partial char enum dataset succeeded.\n' );

%--------------------------------------------------------------------------
%function test_retrieveFullEnumDataset()
%
% In the example file, we have 
%
% HDF5 "data/h5ex_t_enum.h5" {
%  GROUP "/" {
%     DATASET "DS1" {
%        DATATYPE  H5T_ENUM {
%              H5T_STD_I16BE;
%              "SOLID"            0;
%              "LIQUID"           1;
%              "GAS"              2;
%              "PLASMA"           3;
%           }
%        DATASPACE  SIMPLE { ( 4, 7 ) / ( 4, 7 ) }
%        DATA {
%        (0,0): SOLID, SOLID, SOLID, SOLID, SOLID, SOLID, SOLID,
%        (1,0): SOLID, LIQUID, GAS, PLASMA, SOLID, LIQUID, GAS,
%        (2,0): SOLID, GAS, SOLID, GAS, SOLID, GAS, SOLID,
%        (3,0): SOLID, PLASMA, GAS, LIQUID, SOLID, PLASMA, GAS
%        }
%     }
%  }
%  }


%actData = h5varget('data/h5ex_t_enum.h5','/DS1','enumerate',true);

%load('th5varget.mat');
%if ~isequal(actData,expData.fullEnum)
%	error ( 'failed to verify full enum dataset retrieval');
%end

%fprintf ( 1, 'Retrieval of a full enum dataset succeeded.\n' );

%--------------------------------------------------------------------------
function test_openFileHandle()
%
% In the example file, we have 
%
% 
% HDF5 "h5ex_t_vlen.h5" {
% GROUP "/" {
%    DATASET "DS1" {
%       DATATYPE  H5T_VLEN { H5T_STD_I32LE}
%       DATASPACE  SIMPLE { ( 2 ) / ( 2 ) }
%       DATA {
%       (0): (3, 2, 1), (1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144)
%       }
%    }
% }
% }

fid = H5F.open('data/h5ex_t_vlen.h5','H5F_ACC_RDONLY','H5P_DEFAULT');
actData = h5varget(fid,'/DS1',0,1);
H5F.close(fid);

expData = cell(1);
expData{1} = int32([3 2 1]');
if ~isequal(actData,expData)
	error ( 'failed to verify full dataset retrieval');
end


%--------------------------------------------------------------------------
function test_retrievePartialVlenInt32Dataset()
%
% In the example file, we have 
%
% 
% HDF5 "h5ex_t_vlen.h5" {
% GROUP "/" {
%    DATASET "DS1" {
%       DATATYPE  H5T_VLEN { H5T_STD_I32LE}
%       DATASPACE  SIMPLE { ( 2 ) / ( 2 ) }
%       DATA {
%       (0): (3, 2, 1), (1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144)
%       }
%    }
% }
% }

actData = h5varget('data/h5ex_t_vlen.h5','/DS1',0,1);

expData = cell(1);
expData{1} = int32([3 2 1]');
if ~isequal(actData,expData)
	error ( 'failed to verify full dataset retrieval');
end


%--------------------------------------------------------------------------
function test_retrieveStridedVlenStringDataset()
%
% In the example file, we have 
%
% 
%  
%  HDF5 "h5ex_t_vlstring.h5" {
%  GROUP "/" {
%     DATASET "DS1" {
%        DATATYPE  H5T_STRING {
%              STRSIZE H5T_VARIABLE;
%              STRPAD H5T_STR_SPACEPAD;
%              CSET H5T_CSET_ASCII;
%              CTYPE H5T_C_S1;
%           }
%        DATASPACE  SIMPLE { ( 4 ) / ( 4 ) }
%        DATA {
%        (0): "Parting", "is such", "sweet", "sorrow."
%        }
%     }
%  }
%  }

actData = h5varget('data/h5ex_t_vlstring.h5','/DS1',1,2,2);

expData = cell(2,1);
expData{1} = 'is such';
expData{2} = 'sorrow.';
if ~isequal(actData,expData)
	error ( 'failed to verify strided vlen string retrieval');
end


%--------------------------------------------------------------------------
function test_retrievePartialVlenStringDataset()
%
% In the example file, we have 
%
% 
%  
%  HDF5 "h5ex_t_vlstring.h5" {
%  GROUP "/" {
%     DATASET "DS1" {
%        DATATYPE  H5T_STRING {
%              STRSIZE H5T_VARIABLE;
%              STRPAD H5T_STR_SPACEPAD;
%              CSET H5T_CSET_ASCII;
%              CTYPE H5T_C_S1;
%           }
%        DATASPACE  SIMPLE { ( 4 ) / ( 4 ) }
%        DATA {
%        (0): "Parting", "is such", "sweet", "sorrow."
%        }
%     }
%  }
%  }

actData = h5varget('data/h5ex_t_vlstring.h5','/DS1',1,2);

expData = cell(2,1);
expData{1} = 'is such';
expData{2} = 'sweet';
if ~isequal(actData,expData)
	error ( 'failed to partial vlen string retrieval');
end


%--------------------------------------------------------------------------
function test_retrieveFullVlenStringDataset()
%
% In the example file, we have 
%
% 
%  
%  HDF5 "h5ex_t_vlstring.h5" {
%  GROUP "/" {
%     DATASET "DS1" {
%        DATATYPE  H5T_STRING {
%              STRSIZE H5T_VARIABLE;
%              STRPAD H5T_STR_SPACEPAD;
%              CSET H5T_CSET_ASCII;
%              CTYPE H5T_C_S1;
%           }
%        DATASPACE  SIMPLE { ( 4 ) / ( 4 ) }
%        DATA {
%        (0): "Parting", "is such", "sweet", "sorrow."
%        }
%     }
%  }
%  }

actData = h5varget('data/h5ex_t_vlstring.h5','/DS1');

expData = cell(4,1);
expData{1} = 'Parting';
expData{2} = 'is such';
expData{3} = 'sweet';
expData{4} = 'sorrow.';
if ~isequal(actData,expData)
	error ( 'failed to verify full vlen string retrieval');
end


%--------------------------------------------------------------------------
function test_retrieveFullVlenInt32Dataset()
%
% In the example file, we have 
%
% 
% HDF5 "h5ex_t_vlen.h5" {
% GROUP "/" {
%    DATASET "DS1" {
%       DATATYPE  H5T_VLEN { H5T_STD_I32LE}
%       DATASPACE  SIMPLE { ( 2 ) / ( 2 ) }
%       DATA {
%       (0): (3, 2, 1), (1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144)
%       }
%    }
% }
% }

actData = h5varget('data/h5ex_t_vlen.h5','/DS1');

expData = cell(2,1);
expData{1} = int32([3 2 1]');
expData{2} = int32([1 1 2 3 5 8 13 21 34 55 89 144]');
if ~isequal(actData,expData)
	error ( 'failed to verify full dataset retrieval');
end


%--------------------------------------------------------------------------
function test_retrieveFullSinglePrecisionDataset()
% Retrieve an entire dataset.  
%
% In the example file, we have 
%
%         Group:  /g2
%             Dataset /g2/dset2.1
%                 Datatype:  
%                     Class:  H5T_IEEE_F32BE
%                 Dims:  10 
%                 MaxDims:  10 
%                 Layout:  contiguous

actData = h5varget('example.h5','/g2/dset2.1');
expData = single([1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9]');
if ~isequal(actData,expData)
	error ( 'failed to verify full dataset retrieval');
end



%--------------------------------------------------------------------------
function test_retrievePartialContiguousSinglePrecisionDataset()
% Retrieve partial contiguous dataset.  
%
% In the example file, we have 
%
%         Group:  /g2
%             Dataset /g2/dset2.1
%                 Datatype:  
%                     Class:  H5T_IEEE_F32BE
%                 Dims:  10 
%                 MaxDims:  10 
%                 Layout:  contiguous

actData = h5varget('example.h5','/g2/dset2.1',1, 3);
expData = single([1.1 1.2 1.3]');
if ~isequal(actData,expData)
	error ( 'failed to verify full dataset retrieval');
end


%--------------------------------------------------------------------------
function test_retrievePartialContiguousCompoundDataset()
% Retrieve one element
%

switch( version('-release') )
	case { '2008b', '2008a', '2007b', '2007a', '2006b' }
		fprintf ( 'Filtering this t_retrievePartialContiguousCompoundDataset test for now, bug in variable length compound fields...\n' );
		return
end

actData = h5varget('data/h5ex_t_cmpd.h5','/DS1',2,2);
load('th5varget.mat');


if ~isequal(actData,expData.contiguousCompound)
	error ( 'failed to verify partial compound dataset retrieval');
end


%--------------------------------------------------------------------------
function test_retrievePartialStridedSinglePrecisionDataset()
% Retrieve strided dataset.  
%
% In the example file, we have 
%
%         Group:  /g2
%             Dataset /g2/dset2.1
%                 Datatype:  
%                     Class:  H5T_IEEE_F32BE
%                 Dims:  10 
%                 MaxDims:  10 
%                 Layout:  contiguous

actData = h5varget('example.h5','/g2/dset2.1',1, 3, 2);
expData = single([1.1 1.3 1.5]');
if ~isequal(actData,expData)
	error ( 'failed to verify full dataset retrieval');
end

%--------------------------------------------------------------------------
function test_fillValueWithInt32Dataset()
% 
% Only double precision and single precision datasets see any effect from
% a fill value.  The retrieved value should be equal to the fill value, 
% which is 99.

actData = h5varget('data/h5ex_d_fillval.h5','/DS1',[9 5], [1 1]);
expData = int32(99);
if ~isequal(actData,expData)
	error ( 'failed to verify integer dataset with fill value');
end


%--------------------------------------------------------------------------
function test_fillValueWithDoubleDataset()
% 
% Only double precision and single precision datasets see any effect from
% a fill value.  The retrieved value should be equal to NaN

actData = h5varget('data/h5ex_d_dfillval.h5','/DS1',[9 5], [1 1]);
if ~isnan(actData)
	error ( 'failed to verify double dataset with fill value');
end


%--------------------------------------------------------------------------
function test_retrieveTooMuchData()
% Try to retrieve too much data.
%
% In the example file, we have 
%
%         Group:  /g2
%             Dataset /g2/dset2.1
%                 Datatype:  
%                     Class:  H5T_IEEE_F32BE
%                 Dims:  10 
%                 MaxDims:  10 
%                 Layout:  contiguous

try
    h5varget('example.h5','/g2/dset2.1',5, 10);
catch me
    if ~strcmp(me.identifier,'HDF5TOOLS:H5VARGET:invalidFilespace')
        error ( 'failed to verify trying to retrieve too much of a 1D dataset.' );
    end
    return
end
error('failed to verify trying to retrieve too much of a 1D dataset.');

