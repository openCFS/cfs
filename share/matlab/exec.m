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
%      -Last update: 15 Nov 2007
%      -Revision:    0.1
%      -Author:      Jens Grabinger
%
% ==============================================================

function [status, result] = exec(cmd)

% contruct temporary filename
tmpfile = sprintf('.exec%d.sh', ceil(666*rand));

% write shell script to temp file
fid = fopen(tmpfile, 'w');
fprintf(fid, '#!/bin/sh\n\n');

fprintf(fid, '# Split LD_LIBRARY_PATH which has been augmented by Matlab\n');
fprintf(fid, 'PATHS=`echo $LD_LIBRARY_PATH | sed "s/:/ /g"`\n\n');

fprintf(fid, '# Unset old lib path\n');
fprintf(fid, 'unset LD_LIBRARY_PATH\n\n');

fprintf(fid, '# Reverse the ordering of the old library path so that Matlab paths are the last ones\n');
fprintf(fid, 'for p in $PATHS\n');
fprintf(fid, 'do\n');
fprintf(fid, '  LD_LIBRARY_PATH="$p:$LD_LIBRARY_PATH"\n');
fprintf(fid, 'done\n\n');

fprintf(fid, '# Determine machine type and add standard paths in front of lib path\n');
fprintf(fid, 'MACH=`uname -m`\n');
fprintf(fid, 'case $MACH in\n');
fprintf(fid, '           i[3-6]86)\n');
fprintf(fid, '              LD_LIBRARY_PATH="/lib:/usr/lib:$LD_LIBRARY_PATH"\n');
fprintf(fid, '              ;;\n');
fprintf(fid, '           x86_64)\n');
fprintf(fid, '              LD_LIBRARY_PATH="/lib64:/usr/lib64:$LD_LIBRARY_PATH"\n');
fprintf(fid, '              ;;\n');
fprintf(fid, 'esac\n\n');

fprintf(fid, '# Set new lib path\n');
fprintf(fid, 'export LD_LIBRARY_PATH\n\n');

fprintf(fid, '# Execute command\n');
fprintf(fid, '%s\n', cmd);
fclose(fid);

% execute script
[status,result] = system(sprintf('%s %s', getenv('SHELL'), tmpfile));

% delete temp file
delete(tmpfile);
