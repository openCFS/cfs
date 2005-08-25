function mat = sbmread( basename, nrows, ncols )

% Import Super-Block-Matrix 
%
% -----------------------------------------------------------------------------
% mat = sbmread( basename, nrows, ncols )
%
% The function imports a Super-Block-Matrix. It assumes that the sub-matrices
% of the SBM matrix are stored in separate files following the MatrixMarket
% format specification. It first imports the sub-matricies and then assembles
% them into a large sparse matrix.
% 
% Note: The function relies on mminfo.m and mmread.m for importing the
%       sub-matrices
%
% The parameters have the following meaning
% 
% basename .... basename of the files containing the sub-matrices; we expect
%               the sub-matric at position (i,j) to be contained in a file
%               named <basename>_i_j.mtx
% nrows ....... number of block rows
% ncols ....... number of block columns
% -----------------------------------------------------------------------------


% startup message
disp( ' ------------------------------------------------------------------' );
disp( '  Importing SBM Matrix' );
disp( ' ------------------------------------------------------------------' );

% determine size information
disp( ' ' );
disp( ' STEP 1:' );
disp( ' ' );

rowStart(1) = 1;
colStart(1) = 1;
sumEntries = 0;

for i = 1:nrows

  j = 1;
  fname = sprintf( '%s_%d_1.mtx', basename, i );
  mesg = sprintf( ' -> Reading info from %s', fname );
  disp( mesg );
  [rows, cols, entries, rep, field, symmetry] = mminfo( fname );
  rowStart(i+1) = rowStart(i) + rows;
  colStart(j+1) = colStart(j) + cols;

  for j = 2:ncols
    fname = sprintf( '%s_%d_%d.mtx', basename, i, j );
    mesg = sprintf( ' -> Reading info from %s', fname );
    disp( mesg );
    [rows, cols, entries, rep, field, symmetry] = mminfo( fname );
    if i == 1
      colStart(j+1) = colStart(j) + cols;
    end
    
    sumEntries = sumEntries + entries;
  end
end

myRows = rowStart( nrows ) - 1;
myCols = colStart( ncols ) - 1;
disp( ' ' );
mesg = sprintf( ' Block matrix is of dimension %d x %d with nnz = %d', ...
    myRows, myCols, sumEntries );
disp(mesg);

% Generate sparse SBM matrix
mat = spalloc( myRows, myCols, sumEntries );

% importing sub-matrices
disp( ' ' );
disp( ' STEP 2:' );
disp( ' ' );
for i = 1:nrows
  for j = 1:ncols
    fname = sprintf( '%s_%d_%d.mtx', basename, i, j );
    mesg = sprintf( ' -> Importing sub-matrix (%d,%d) from %s', ...
	i, j, fname );
    disp(mesg);
    mesg = sprintf( '    to sub-block [%d:%d,%d:%d]\n', ...
	rowStart(i), rowStart(i+1) - 1, ...
	colStart(j), colStart(j+1) - 1 );
    disp(mesg);
    mat(rowStart(i):rowStart(i+1)-1,colStart(j):colStart(j+1)-1) = ...
	mmread( fname );
  end
end

% finished
disp( ' ------------------------------------------------------------------' );
