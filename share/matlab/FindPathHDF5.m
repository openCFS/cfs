%    FindPathHDF5
%
%    Searches a HDF5 file for the path to a given result.
%
% Input Parameters
%   * toplevel  - Top level node (/) of a HDF5 file
%   * multistep - which multistep to search for
%   * step      - which step to search for
%   * quantity  - what quantity
%   * region    - on what region
%      
% Return Values
%   * found    - integer that indicates at which node level the search
%                  was stopped.
%
%                  toplevel                             found=0
%                  |
%                  `-Results                            found=1
%                    |
%                    `-Mesh                             found=2
%                      |
%                      `-MultiStep_<n>                  found=3
%                        |
%                        `-Step_<t>                     found=4
%                          |
%                          `-<quantity>                 found=5
%                            |
%                            `-<region>                 found=6
%                              |
%                              `-[ Nodes | Elements ]   found=7
%                                |
%                                `-Real                 found=8
%
%   * resgroup - group node where the search ended (according to output
%                  parameter 'found', if found >= 3). Does NOT point to
%                  dataset 'Real' even if found == 8.
%   * restype  - if found >= 7, restype is the name of the group containing
%                  the Real dataset (currently either 'Nodes' or 'Elements')
%   * msgroup  - node of the multistep group (valid if found >= 3)
%   * datafile - the data file's top level group (used with external files).
%                  Equals the input parameter toplevel, if no external files
%                  are present.
%
% About
%   * Created:  Jun 2007
%   * Authors:  Simon Triebenbacher, Jens Grabinger
%   * Revision: $Id$


function [found resgroup restype msgroup datafile] =  FindPathHDF5(toplevel, multistep, step, quantity, region)

found = 0;
restype = 0;
ext_files = 0;
datafile = toplevel;

% look for MultiStep_x group
number_of_groups = length(toplevel.Groups);
for group=1:number_of_groups
  actgroup = toplevel.Groups(group);
  group_name = actgroup.Name;

  if strcmp(group_name, '/Results')
%    sprintf('We have found the data group')
    ndatagroups = length(actgroup.Groups);

    found = 1;

    for datagroup=1:ndatagroups
      actgroup2 = actgroup.Groups(datagroup);
      group_name = actgroup2.Name;
      if strcmp(group_name, '/Results/Mesh')
%        sprintf('We have found the volume data group')
        nvoldatagroups = length(actgroup2.Groups);

        found = 2;

        try
          ext_files = h5attget(toplevel.Filename, '/Results/Mesh', 'ExternalFiles');
        catch %#ok<CTCH>
        end
        
        for voldatagroup=1:nvoldatagroups
          actgroup3 = actgroup2.Groups(voldatagroup);
          group_name = actgroup3.Name;

          cmpstr = sprintf('/Results/Mesh/MultiStep_%d', multistep);
%         if anybody reads this line, he gets a beer from me (zoerner)
          if strcmp(group_name, cmpstr)
%            sprintf('We have found the multistep %d group', multistep)
            basepath = cmpstr;
            found = 3;
            break
          end
        end
        break
      end
    end
    break
  end
end

if found < 3
  return
end
msgroup = actgroup3;
resgroup = msgroup;

% look for time step
nmstepgroups = length(actgroup3.Groups);
for mstepgroup=1:nmstepgroups
  actgroup = actgroup3.Groups(mstepgroup);
  group_name = actgroup.Name;

  cmpstr = sprintf('%s/Step_%d', basepath, step);
  if strcmp(group_name, cmpstr)
%    sprintf('We have found the Step_%d group', step)
    found = 4;
    break
  end

  cmpstr = sprintf('%s/Step %d', basepath, step);
  if strcmp(group_name, cmpstr)
%    sprintf('We have found the Step %d group', step)
    found = 4;
    break
  end

end

if found < 4
  return
end
resgroup = actgroup;
curpath = resgroup.Name;

if ext_files

  fn_attr = h5attget(toplevel.Filename, curpath, 'ExtHDF5FileName');
  ext_filename = fn_attr{1};
  
  % extract path to input file, because external time step files are given
  % relative to this path
  pos = strfind(toplevel.Filename, '/');
  if size(pos, 2)
    prefix = toplevel.Filename(1:pos(size(pos, 2)));
    if prefix(length(prefix)) == '/'
      ext_filename = sprintf('%s%s', prefix, ext_filename);
    else
      ext_filename = sprintf('%s/%s', prefix, ext_filename);
    end
  end

  if exist(ext_filename,'file') == 2
    try
      df_info = hdf5info(ext_filename);
    catch %#ok<CTCH>
      return
    end
    datafile = df_info.GroupHierarchy;
    actgroup = datafile;
    curpath = '';
  else
    return
  end
end

% look for the quantity group
cmpstr = sprintf('%s/%s', curpath, quantity);
nstepgroups = length(actgroup.Groups);
for stepgroup=1:nstepgroups
  actgroup2 = actgroup.Groups(stepgroup);
  group_name = actgroup2.Name;

  if strcmp(group_name, cmpstr)
%    sprintf('We have found the quantity %s group', quantity)
    found = 5;
    break
  end

end

if found < 5
  return
end
resgroup = actgroup2;
curpath = resgroup.Name;

% look for the region group
cmpstr = sprintf('%s/%s', curpath, region);
nqgroups = length(actgroup2.Groups);
for qgroup=1:nqgroups
  actgroup = actgroup2.Groups(qgroup);
  group_name = actgroup.Name;

  if strcmp(group_name, cmpstr)
%    sprintf('We have found the region group %s', region)
    found = 6;
    break
  end
end

if found < 6
  return
end
resgroup = actgroup;
curpath = resgroup.Name;

% search resultDescription for type of result
cmpstr = sprintf('%s/ResultDescription', basepath);
nmstepgroups = length(msgroup.Groups);
for mstepgroup=1:nmstepgroups
  actgroup3 = msgroup.Groups(mstepgroup);
  group_name = actgroup3.Name;

  if strcmp(group_name, cmpstr)
%    sprintf('We have found the ResultDescription group')
    cmpstr = sprintf('%s/%s', group_name, quantity);

    nrdgroups = length(actgroup3.Groups);
    for rdgroup=1:nrdgroups
      actrdgroup = actgroup3.Groups(rdgroup);
      group_name = actrdgroup.Name;

      if strcmp(group_name, cmpstr)
%        sprintf('We have found the ResultDescription/%s group', quantity)
        cmpstr = sprintf('%s/DefinedOn', group_name);

        ndatasets = length(actrdgroup.Datasets);
        for dataset=1:ndatasets
          ds = actrdgroup.Datasets(dataset);
          ds_name = ds.Name;

          if strcmp(ds_name, cmpstr)
            restype = h5varget(ds.Filename, ds_name);
            break
          end
        end

        break
      end
    end

    break
  end
end

switch restype
case 0 % not found
%  errorstr = sprintf('result type of quantity %s not found.', quantity);
%  error(errorstr);
  return
case 1 % nodes
  path_suffix = 'Nodes';
case 4 % elements
  path_suffix = 'Elements';
otherwise
  error('result type %d not supported.', restype);
end

% look for Nodes/Elements group
cmpstr = sprintf('%s/%s', curpath, path_suffix);
nreggroups = length(actgroup.Groups);
for reggroup=1:nreggroups
  actgroup2 = actgroup.Groups(reggroup);
  group_name = actgroup2.Name;

  if strcmp(group_name, cmpstr)
%    sprintf('We have found the %s group', path_suffix)
    found = 7;
    break
  end
end

if found < 7
  return
end
resgroup = actgroup2;
curpath = resgroup.Name;

%look for the Real dataset
ndatasets = length(actgroup2.Datasets);
for dataset=1:ndatasets
  ds = actgroup2.Datasets(dataset);
  ds_name = ds.Name;

  cmpstr = sprintf('%s/%s', curpath, 'Real');
  if strcmp(ds_name, cmpstr)
%    sprintf('We have found the dataset')
    found = 8;
    break
  end

end
