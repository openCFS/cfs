clear, clc;

e1 = extractvoltage('band_stop-elecNetworkPotential-node-1-E1.hist');
e2 = extractvoltage('band_stop-elecNetworkPotential-node-2-E1.hist');
e3 = extractvoltage('band_stop-elecNetworkPotential-node-3-E4.hist');
e4 = extractvoltage('band_stop-elecNetworkPotential-node-4-E2.hist');
t = e1(:,1);

% resonance behaviour (T_res) via inductor voltage
figure(1), clf
semilogx(e1(:,1),e1(:,2)-e4(:,2)), hold on, grid minor

function data = extractvoltage(fname)

fid = fopen(fname);
data = textscan(fid,'%f  %f %f','HeaderLines',3) ;
data = cell2mat(data);
fclose(fid);

end

