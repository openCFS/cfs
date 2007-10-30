% ==============================================================
%
%    FindPathHDF5
%
%
%    GENERAL
%    Searches a HDF5 file for the path to a given result.
%
%    INPUT/S
%      toplevel  - Top level node (/) of a HDF5 file
%      multistep - which multistep to search for
%      step      - which step to search for
%      quantity  - what quantity
%      region    - on what region
%        
%    OUTPUT/S
%      found     - integer that indicates at which node level the search
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
%      resgroup  - group node where the search ended (according to output
%                  parameter 'found', if found >= 3). Does NOT point to
%                  dataset 'Real' even if found == 8.
%      msgroup   - node of the multistep group (valid if found >= 3)
%
%      
%    ABOUT
%
%      -Created:     Jun 2007
%      -Last update: 30 Oct 2007
%      -Revision:    0.2
%      -Authors:     Simon Triebenbacher, Jens Grabinger
%
% ==============================================================


function [found resgroup restype msgroup] =  FindPathHDF5(toplevel, multistep, step, quantity, region)

found = 0;

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

        for voldatagroup=1:nvoldatagroups
          actgroup3 = actgroup2.Groups(voldatagroup);
          group_name = actgroup3.Name;

          cmpstr = sprintf('/Results/Mesh/MultiStep_%d', multistep);
          if strcmp(group_name, cmpstr)
%            sprintf('We have found the multistep %d group', multistep)
            basepath = cmpstr;
            found = 3;
            break;
          end

          cmpstr = sprintf('/Results/Mesh/MultiStep %d', multistep);
          if strcmp(group_name, cmpstr)
%            sprintf('We have found the multistep %d group', multistep)
            basepath = cmpstr;
            found = 3;
            break;
          end
        end
        break;
      end
    end
    break;
  end
end

if found < 3
  return;
end
msgroup = actgroup3;
resgroup = msgroup;

%% commented out, because we look for the 'Real' dataset only
%cmpstr = sprintf('%s/AnalysisType', basepath);
%nmstepattrs = length(actgroup3.Attributes);
%for mstepattr=1:nmstepattrs
%  actattr = actgroup3.Attributes(mstepattr);
%  attr_name = actattr.Name;
%  if strcmp(attrname, cmpstr)
%    ana_type = hdf5read(actattr.Filename, attr_name);
%    break;
%  end
%end

% look for time step
nmstepgroups = length(actgroup3.Groups);
for mstepgroup=1:nmstepgroups
  actgroup = actgroup3.Groups(mstepgroup);
  group_name = actgroup.Name;

  cmpstr = sprintf('%s/Step_%d', basepath, step);
  if strcmp(group_name, cmpstr)
%    sprintf('We have found the step %d group', step)
    stepname = sprintf('Step_%d', step);
    found = 4;
    break;
  end

  cmpstr = sprintf('%s/Step %d', basepath, step);
  if strcmp(group_name, cmpstr)
%    sprintf('We have found the step %d group', step)
    stepname = sprintf('Step %d', step);
    found = 4;
    break;
  end

end

if found < 4
  return;
end
resgroup = actgroup;

% look for the quantity group
cmpstr = sprintf('%s/%s/%s', basepath, stepname, quantity);
nstepgroups = length(actgroup.Groups);
for stepgroup=1:nstepgroups
  actgroup2 = actgroup.Groups(stepgroup);
  group_name = actgroup2.Name;

  if strcmp(group_name, cmpstr)
%    sprintf('We have found the quantity %s group', quantity)
    found = 5;
    break;
  end

end

if found < 5
  return;
end
resgroup = actgroup2;

% look for the region group
cmpstr = sprintf('%s/%s/%s/%s', basepath, stepname, quantity, region);
nqgroups = length(actgroup2.Groups);
for qgroup=1:nqgroups
  actgroup = actgroup2.Groups(qgroup);
  group_name = actgroup.Name;

  if strcmp(group_name, cmpstr)
%    sprintf('We have found the region group %s', region)
    found = 6;
    break;
  end
end

if found < 6
  return;
end
resgroup = actgroup;

% search resultDescription for type of result
restype = 0;
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
            restype = hdf5read(ds.Filename, ds_name);
            break;
          end
        end

        break;
      end
    end

    break;
  end
end

switch restype
case 0 % not found
%  errorstr = sprintf('result type of quantity %s not found.', quantity);
%  error(errorstr);
  return;
case 1 % nodes
  path_suffix = 'Nodes';
case 4 % elements
  path_suffix = 'Elements';
otherwise
  errorstr = sprintf('result type %d not supported.', restype);
  error(errorstr);
  return;
end

% look for Nodes/Elements group
cmpstr = sprintf('%s/%s/%s/%s/%s', basepath, stepname, quantity, region, path_suffix);
nreggroups = length(actgroup.Groups);
for reggroup=1:nreggroups
  actgroup2 = actgroup.Groups(reggroup);
  group_name = actgroup2.Name;

  if strcmp(group_name, cmpstr)
%    sprintf('We have found the %s group', path_suffix)
    found = 7;
    break;
  end
end

if found < 7
  return;
end
resgroup = actgroup2;

%look for the Real dataset
ndatasets = length(actgroup2.Datasets);
for dataset=1:ndatasets
  ds = actgroup2.Datasets(dataset);
  ds_name = ds.Name;

  cmpstr = sprintf('%s/%s/%s/%s/%s/%s', basepath, stepname, quantity, region, path_suffix, 'Real');
  if strcmp(ds_name, cmpstr)
%    sprintf('We have found the dataset')
    found = 8;
    break;
  end

end
