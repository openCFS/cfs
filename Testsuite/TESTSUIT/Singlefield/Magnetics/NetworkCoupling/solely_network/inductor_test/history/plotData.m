clear, clc;

e1 = extractvoltage('inductor_test-elecNetworkPotential-node-1-E1.hist');
e2 = extractvoltage('inductor_test-elecNetworkPotential-node-2-E1.hist');
e3 = extractvoltage('inductor_test-elecNetworkPotential-node-3-E4.hist');
e4 = extractvoltage('inductor_test-elecNetworkPotential-node-4-E2.hist');
t = e1(:,1);

figure(1), clf, hold on, grid minor
plot(t,e1(:,2)-e4(:,2))


figure(2), clf, hold on, grid minor
plot(t,e4(:,2)-e2(:,2))

function data = extractvoltage(fname)

fid = fopen(fname);
data = textscan(fid,'%f  %f','HeaderLines',3) ;
data = cell2mat(data);
fclose(fid);

end