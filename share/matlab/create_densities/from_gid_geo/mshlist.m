function [ list ] = mshlist( filter, dirname )
%mshlist Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    dirname = 'msh/';
end

if nargin < 1
    msh_name = strcat(dirname, '*.msh');
else
    msh_name = strcat(dirname, filter);
end

list=get_sort_dir(msh_name);

for i=1:length(list)
    list{i} = strcat(dirname, list{i});
end

list = list';

return;

function [filelist]=get_sort_dir(what)
% gets sorted list of files specified with suffix 'what', checks their size
% and if they are not directories

  list=dir(what);
  n=length(list);
  reallist=[];
  pos=1;
  for i=1:n
    if (list(i).bytes>0 && ~list(i).isdir)  % && exist(list(i).name,'file'))
      reallist{pos}=list(i).name;
      pos = pos+1;
    end
  end
  filelist=sort(reallist);
  
return;
