% exec
%
% Executes a shell command like MATLAB's system function, but uses a
% clean LD_LIBRARY_PATH.
%
% Input Parameters
%   * cmd       - shell command to be executed
%     
% Return Values
%   * status    - the shell's exit status
%   * result    - the command's output (stdout)
%   
% About
%   * Created:  30 Oct 2007
%   * Authors:  Simon Triebenbacher, Jens Grabinger
%   * Revision: $Id$


function [status, result] = exec(cmd)

% contruct temporary filename
rng('shuffle')
tmpfile = sprintf('.exec%d.sh', randi(999999));

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
fprintf(fid, 'case $HOSTTYPE in\n');
fprintf(fid, '           i[3-6]86)\n');
fprintf(fid, '              LD_LIBRARY_PATH="/lib:/usr/lib:/lib/$HOSTTYPE-linux-gnu:/usr/lib/$HOSTTYPE-linux-gnu:$LD_LIBRARY_PATH"\n');
fprintf(fid, '              ;;\n');
fprintf(fid, '           x86_64)\n');
fprintf(fid, '              LD_LIBRARY_PATH="/lib64:/usr/lib64:/lib/$HOSTTYPE-linux-gnu:/usr/lib/$HOSTTYPE-linux-gnu:$LD_LIBRARY_PATH"\n');
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
