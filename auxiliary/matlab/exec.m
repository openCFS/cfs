% ==============================================================
%
%    exec
%
%
%    GENERAL
%    Executes a shell command like MATLAB's system function, but uses a
%    clean LD_LIBRARY_PATH.
%
%    INPUT/S
%      cmd       - shell command to be executed
%        
%    OUTPUT/S
%      status    - the shell's exit status
%      result    - the command's output (stdout)
%      
%    ABOUT
%
%      -Created:     30 Oct 2007
%      -Last update: 07 Nov 2007
%      -Revision:    0.1
%      -Author:      Jens Grabinger
%
% ==============================================================

function [status, result] = exec(cmd)

% contruct temporary filename
tmpfile = sprintf('.exec%d.sh', ceil(666*rand));

% write shell script to temp file
fid = fopen(tmpfile, 'w');
fprintf(fid, '#!/bin/sh\nunset LD_LIBRARY_PATH\n');
fprintf(fid, 'source /home/data/libraries/CFSDEPS/import_cfsdeps.sh\n');
fprintf(fid, '%s\n', cmd);

% execte script
[status,result] = system(sprintf('%s %s', getenv('SHELL'), tmpfile));

% delete temp file
delete(tmpfile);
