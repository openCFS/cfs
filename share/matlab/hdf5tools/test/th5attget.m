function th5attget()

%   Copyright 2008-2009 The MathWorks, Inc.

fprintf('\tTesting H5ATTGET...  ');

test_3string;  % retrieve a non-scalar string attribute
test_arrayatt; % retrieve scalar nullterm string attribute
test_scalar_nullterm_string; % retrieve scalar nullterm string attribute
test_attrNotThere;           % try to retrieve an attribute that is not there

fprintf ('OK\n');

%--------------------------------------------------------------------------
function test_3string()
% Retrieve a non scalar string attribute
%
%
% GROUP "/" {
%    GROUP "g1" {
%       DATASET "dset1" {
%          ATTRIBUTE "attr1" {
%             DATATYPE  H5T_STRING {
%                   STRSIZE 11;
%                   STRPAD H5T_STR_NULLTERM;
%                   CSET H5T_CSET_ASCII;
%                   CTYPE H5T_C_S1;
%                }
%             DATASPACE  SIMPLE { ( 3 ) / ( 3 ) }
%             DATA {
%             (0): "0123456789", "abcdefghij", "ABCDEFGHIJ"
%             }
%          }
% 
actData = h5attget('data/hstr2.h5','/g1/dset1', 'attr1');
expData = ['0123456789' char(0); ...
    'abcdefghij' char(0); ...
    'ABCDEFGHIJ', char(0)]';
if ~isequal(actData,expData) && isa(actData,'char') 
	error ( 'failed to verify non scalar string attribute retrieval');
end


%--------------------------------------------------------------------------
function test_arrayatt()
actData = h5attget('data/h5ex_t_arrayatt.h5','/DS1', 'A1');

if ~isa(actData,'int64')
	error ( 'failed to verify int64 array attribute datatype');
end

if ~isequal(size(actData),[5 3 4])
	error ( 'failed to verify size of int64 array attribute datatype');
end

load('th5attget.mat');

if ~isequal(actData,expData.arrayAtt)
	error ( 'failed to verify int64 array attribute');
end


%--------------------------------------------------------------------------
function test_scalar_nullterm_string()
% Retrieve a nullterm string attribute
%
% In the example file, we have 
%
%    ATTRIBUTE "attr5" {
%       DATATYPE  H5T_STRING {
%             STRSIZE 17;
%             STRPAD H5T_STR_NULLTERM;
%             CSET H5T_CSET_ASCII;
%             CTYPE H5T_C_S1;
%          }
%       DATASPACE  SCALAR
%       DATA {
%       (0): "string attribute"
%       }
%    }
% 
actData = h5attget('data/hattr.h5','/', 'attr5');
expData = { 'string attribute' };
if ~isequal(actData,expData) && isa(actData,'cell') 
	error ( 'failed to verify single precision attribute retrieval');
end



%--------------------------------------------------------------------------
function test_attrNotThere()
% Try to retrieve an attribute that is not there.
%
% 

try
	h5attget('example.h5','/','attr3');
catch me
    switch(version('-release'))
        case { '2008b', '2008a', '2007b', '2007a', '2006b' }
            expEid = 'MATLAB:H5ML_hdf5:invalidID';
        case '2009a'
            expEid = 'MATLAB:hdf5lib2:H5Aopen_name:attributeOpenFailed';
        otherwise
            expEid = 'MATLAB:hdf5lib2:H5Aopen_name:failure';
    end
    if strcmp(me.identifier,expEid)
        return
    end
end


