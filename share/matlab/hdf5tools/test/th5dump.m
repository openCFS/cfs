function th5dump()

%   Copyright 2008-2009 The MathWorks, Inc.

fprintf ('\tTesting H5DUMP...  ' );
test_compound;  
test_dataset_reference;  
test_bitfield;  
fprintf ('OK\n' );
return


%--------------------------------------------------------------------------
function test_compound()
% Dump example.h5

load('th5dump.mat');
output_text = evalc('h5dump(''data/h5ex_t_cmpd.h5'');');

v = version('-release');
switch(v)
	case { '2007b', '2008a', '2008b','2009a'}
		if ~strcmp(d.R2007b.cmpd,output_text)
		    error('failed');
		end
	otherwise
		if ~strcmp(d.cmpd,output_text)
		    error('failed');
		end

end





%--------------------------------------------------------------------------
function test_dataset_reference()
% Test a reference dataset dump.

load('th5dump.mat');
cmd = 'h5dump(''data/h5ex_t_objref.h5'',''/DS1'');';
output_text = evalc(cmd);

v = version('-release');
switch(v)
	case { '2007b', '2008a', '2008b','2009a'}
		if ~strcmp(d.R2007b.reference_dataset_dump,output_text)
		    error('failed');
		end
	otherwise
		if ~strcmp(d.reference_dataset_dump,output_text)
		    error('failed');
		end

end







%--------------------------------------------------------------------------
function test_bitfield()
% Test a bitfield dump.

load('th5dump.mat');

cmd = 'h5dump(''data/h5ex_t_bit.h5'',''DS1'');';
output_text = evalc(cmd);

v = version('-release');
switch(v)
	case { '2007b', '2008a', '2008b','2009a'}
		if ~strcmp(d.R2007b.bitfield_dataset_dump,output_text)
		    error('failed');
		end
	otherwise
		if ~strcmp(d.bitfield_dataset_dump,output_text)
		    error('failed');
		end

end



