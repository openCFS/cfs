% ==============================================================
%
%    GENERAL
%    Time data filtering tool for CFS++ HDF5 data.
%
%    INPUT/S
%      toplevel  - Top level node (/) of a HDF5 file
%      multistep - which multistep to search for
%      step      - which step to search for
%      quantity  - what quantity
%      region    - on what region
%        
%    OUTPUT/S
%      found     - found the dataset or any path between...
%      stepname  - 
%      msgroup   - node of the multistep group
%
%      
%    ABOUT
%
%      -Created:     Jun 2007
%      -Last update: Jun 2007
%      -Revision:    0.1
%      -Author:      Simon Triebenbacher
% ==============================================================


function [found stepname msgroup] =  FindPathHDF5(toplevel, multistep, step, quantity, region)

found = 0;

number_of_groups = length(toplevel.Groups);
for group=1:number_of_groups
  actgroup = toplevel.Groups(group);
  group_name = actgroup.Name;

  if strcmp(group_name, '/Data')
%    sprintf('We have found the data group')
    ndatagroups = length(actgroup.Groups);

    found = 1;

    for datagroup=1:ndatagroups
      actgroup2 = actgroup.Groups(datagroup);
      group_name = actgroup2.Name;
      if strcmp(group_name, '/Data/Volume')
%        sprintf('We have found the volume data group')
        nvoldatagroups = length(actgroup2.Groups);

        found = 2;

        for voldatagroup=1:nvoldatagroups
          actgroup3 = actgroup2.Groups(voldatagroup);
          group_name = actgroup3.Name;

          cmpstr = sprintf('/Data/Volume/Multistep_%d', multistep);
          if strcmp(group_name, cmpstr)
%            sprintf('We have found the multistep %d group', multistep)
            basepath = cmpstr;
            found = 3;
            break;
          end

          cmpstr = sprintf('/Data/Volume/Multistep %d', multistep);
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

msgroup = actgroup3;

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

nstepgroups = length(actgroup.Groups);
for stepgroup=1:nstepgroups
  actgroup2 = actgroup.Groups(stepgroup);
  group_name = actgroup2.Name;

  cmpstr = sprintf('%s/%s/%s', basepath, stepname, quantity);
  if strcmp(group_name, cmpstr)
%    sprintf('We have found the quantity %s group', quantity)
    found = 5;
    break;
  end

end

ndatasets = length(actgroup2.Datasets);
for dataset=1:ndatasets
  ds = actgroup2.Datasets(dataset);
  ds_name = ds.Name;

  cmpstr = sprintf('%s/%s/%s/%s', basepath, stepname, quantity, region);
  if strcmp(ds_name, cmpstr)
%    sprintf('We have found the dataset %s', region)
    found = 6;
    break;
  end

end

