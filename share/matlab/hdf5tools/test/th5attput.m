function th5attput()

%   Copyright 2008-2009 The MathWorks, Inc.
fprintf('\tTesting H5ATTPUT...' );

test_writeFileDoesNotExist('notThere.h5');
test_write1DCharAttToGroup('foo.h5');
test_write2DCharAttToVar('foo.h5');
test_writeInt8ToGroup('foo.h5');   
test_writeUint8ToGroup('foo.h5');    
test_writeInt16ToGroup('foo.h5');    
test_writeUint16ToGroup('foo.h5');   
test_writeInt32ToGroup('foo.h5');    
test_writeUint32ToGroup('foo.h5');    
test_writeInt64ToGroup('foo.h5');    
test_writeUint64ToGroup('foo.h5');    
test_009('foo.h5');    % write char attribute to a group

fprintf('OK\n');
delete('foo.h5');


%--------------------------------------------------------------------------
function test_write1DCharAttToGroup(hfile)

create_testfile(hfile);

expData = 'abc';
h5attput(hfile,'/', 'attr1', expData);
actData = h5attget(hfile,'/','attr1');

if ~isequal(actData,expData)
	error ( 'failed to verify');
end



%--------------------------------------------------------------------------
function test_write2DCharAttToVar(hfile)

create_testfile(hfile);

expData = strvcat('abc','def'); %#ok<VCAT>
h5attput(hfile,'/DS1', 'attr1', expData);
actData = h5attget(hfile,'/DS1','attr1');

if ~isequal(actData,expData)
	error ( 'failed to verify');
end


%--------------------------------------------------------------------------
function test_writeFileDoesNotExist(hfile)

expData = 0;
h5attput(hfile,'/', 'test', expData);
actData = h5attget(hfile,'/','test');

if ~isequal(actData,expData)
    error('failed to verify');
end

delete(hfile);



%--------------------------------------------------------------------------
function create_testfile(hfile)

create_plist = H5P.create('H5P_FILE_CREATE');
file_id = H5F.create(hfile, 'H5F_ACC_TRUNC',create_plist , 'H5P_DEFAULT');

space_id = H5S.create_simple(2,[10 8],[10 8]);

dataset_id = H5D.create(file_id,'/DS1','H5T_NATIVE_DOUBLE', space_id, 'H5P_DEFAULT');
H5D.close(dataset_id);

H5F.close(file_id);
H5P.close(create_plist);

%--------------------------------------------------------------------------
function test_writeInt8ToGroup(hfile)
% Write a signed byte attribute
%

create_testfile(hfile);

expData = int8([3 5]);
h5attput(hfile,'/', 'attr1', expData);
actData = h5attget(hfile,'/','attr1');

if ~isequal(actData,expData) && isa(actData,'int8')
	error ( 'failed to verify int8 precision attribute creation');
end




%--------------------------------------------------------------------------
function test_writeUint8ToGroup(hfile)
% Write an unsigned byte attribute
%

create_testfile(hfile);

expData = uint8([3 5]);
h5attput(hfile,'/', 'attr1', expData);
actData = h5attget(hfile,'/','attr1');

if ~isequal(actData,expData) && isa(actData,'uint8')
	error ( 'failed to verify int8 precision attribute creation');
end




%--------------------------------------------------------------------------
function test_writeInt16ToGroup(hfile)
% Write int16 attribute
%

create_testfile(hfile);

expData = int16([3 5]);
h5attput(hfile,'/', 'attr1', expData);
actData = h5attget(hfile,'/','attr1');

if ~isequal(actData,expData) && isa(actData,'int16')
	error ( 'failed to verify int16 precision attribute creation');
end




%--------------------------------------------------------------------------
function test_writeUint16ToGroup(hfile)
% Write uint16 attribute
%

create_testfile(hfile);

expData = uint16([3 5]);
h5attput(hfile,'/', 'attr1', expData);
actData = h5attget(hfile,'/','attr1');

if ~isequal(actData,expData) && isa(actData,'uint16')
	error ( 'failed to verify int16 precision attribute creation');
end




%--------------------------------------------------------------------------
function test_writeInt32ToGroup(hfile)
% Write int32 attribute
%

create_testfile(hfile);

expData = int32([3 5]);
h5attput(hfile,'/', 'attr1', expData);
actData = h5attget(hfile,'/','attr1');

if ~isequal(actData,expData) && isa(actData,'int32')
	error ( 'failed to verify int32 precision attribute creation');
end




%--------------------------------------------------------------------------
function test_writeUint32ToGroup(hfile)
% Write uint32 attribute
%

create_testfile(hfile);

expData = uint32([3 5]);
h5attput(hfile,'/', 'attr1', expData);
actData = h5attget(hfile,'/','attr1');

if ~isequal(actData,expData) && isa(actData,'uint32')
	error ( 'failed to verify int32 precision attribute creation');
end




%--------------------------------------------------------------------------
function test_writeInt64ToGroup(hfile)
% Write int64 attribute
%

create_testfile(hfile);

expData = int64([3 5]);
h5attput(hfile,'/', 'attr1', expData);
actData = h5attget(hfile,'/','attr1');

if ~isequal(actData,expData) && isa(actData,'int64')
	error ( 'failed to verify int64 precision attribute creation');
end



%--------------------------------------------------------------------------
function test_writeUint64ToGroup(hfile)
% Write uint64 attribute
%

create_testfile(hfile);

expData = uint64([3 5]);
h5attput(hfile,'/', 'attr1', expData);
actData = h5attget(hfile,'/','attr1');

if ~isequal(actData,expData) && isa(actData,'uint64')
	error ( 'failed to verify int64 precision attribute creation');
end




%--------------------------------------------------------------------------
function test_009(hfile)
% Write char attribute
%

create_testfile(hfile);

expData = 'abcd';
h5attput(hfile,'/', 'attr1', expData);
actData = h5attget(hfile,'/','attr1');

if ~isequal(actData,expData) && isa(actData,'char')
	error ( 'failed to verify char attribute creation');
end





