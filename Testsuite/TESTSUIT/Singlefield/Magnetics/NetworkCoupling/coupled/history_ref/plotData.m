clear, clc;

e1 = extractvoltage('coupling-elecNetworkPotential-node-1-E3.hist');
e2 = extractvoltage('coupling-elecNetworkPotential-node-2-E1.hist');
e3 = extractvoltage('coupling-elecNetworkPotential-node-3-E2.hist');
e4 = extractvoltage('coupling-elecNetworkPotential-node-4-E3.hist');
t = e1(:,1);

figure(1), clf, hold on, grid minor
plot(t,(e1(:,2)-e4(:,2))/1)


figure(2), clf, hold on, grid minor
% plot(t,e4(:,2)-e2(:,2))
plot(t,e1(:,2))

function data = extractvoltage(fname)

fid = fopen(fname);
data = textscan(fid,'%f  %f','HeaderLines',3) ;
data = cell2mat(data);
fclose(fid);

end