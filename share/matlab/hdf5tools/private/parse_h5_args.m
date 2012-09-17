function args = parse_h5_args ( varargin )
%
% Have to be able to check the following signatures.

% DATA = H5_VARGET(HFILE,VARIABLE) 
% DATA = H5_VARGET(HFILE,VARIABLE,START,COUNT)
% DATA = H5_VARGET(HFILE,VARIABLE,START,COUNT,STRIDE)
%
% DATA = H5_VARGET(...,'enumerate',true);
%
% DATA = H5IMREAD(HFILE) 
% DATA = H5IMREAD(HFILE,VARIABLE) 
% DATA = H5IMREAD(HFILE,VARIABLE,START,COUNT)
% DATA = H5IMREAD(HFILE,VARIABLE,START,COUNT,STRIDE)
% 

% First argument should be the filename.
if ~(ischar(varargin{1}) || isa(varargin{1},'H5ML.id'))
	error ( 'MATLAB:H5TOOLS:badInput', 'File argument must be a filename or an id.')
end

args = struct('filename','', ...
              'dataset', '', ...
              'offset',   [], ...
              'count',   [], ...
              'stride',  [], ...
              'enumerate',  false  );
args.filename = varargin{1};

%
% 2nd argument should be the variable name.
if ~ischar(varargin{2})
	error ( 'MATLAB:H5TOOLS:badInput', 'Variable name argument must be character')
end

args.dataset = varargin{2};


% Was there an enumerate option?
remove_idx = [];
if nargin > 2
	for j = 3:nargin
		if ischar(varargin{j})
			if strcmpi(varargin{j},'enumerate')
				if (nargin ==j)
					error ( 'H5TOOLS:noEnumerateArgumentGiven', ...
					        'No enumerate argument was given.' );
				end
				args.enumerate = varargin{j+1};	
				
				remove_idx = [j j+1];
			end
		end
	end
end
% Remove the pair
varargin(remove_idx) = [];

% Now get the start, count, and stride args.
switch numel(varargin)

case 2
	;
case 3
	error ( 'MATLAB:H5TOOLS:badInput', 'Cannot have 3 input arguments.')

case 4
	%
	% Given the start, stride, and count.
	if ~isnumeric(varargin{3}) 
		error ( 'MATLAB:H5TOOLS:badInput', 'Start argument must be numeric')
	end
	args.offset = varargin{3};

	if ~isnumeric(varargin{4}) 
		error ( 'MATLAB:H5TOOLS:badInput', 'Count argument must be numeric')
	end
	args.count = varargin{4};


	
case 5
	%
	% Given the start, stride, and count.
	if ~isnumeric(varargin{3}) 
		error ( 'MATLAB:H5TOOLS:badInput', 'Start argument must be numeric')
	end
	args.offset = varargin{3};

	if ~isnumeric(varargin{4}) 
		error ( 'MATLAB:H5TOOLS:badInput', 'Count argument must be numeric')
	end
	args.count = varargin{4};

	if ~isnumeric(varargin{5}) 
		error ( 'MATLAB:H5TOOLS:badInput', 'Stride argument must be numeric')
	end
	args.stride = varargin{5};

otherwise
	error ( 'MATLAB:H5TOOLS:badInput', 'Bad number of input arguments.')

end

return






