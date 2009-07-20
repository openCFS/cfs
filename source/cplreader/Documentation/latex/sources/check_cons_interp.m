fileName = 'file.h5'
dsName = '/Results/Mesh/MultiStep_1/Step_X/acouRhsLoad/REGIONNAME/Nodes/Real'

vec = hdf5read(fileName, dsName);

sum = 0
for i=1:length(vec)
  sum = sum+vec(i);
end
