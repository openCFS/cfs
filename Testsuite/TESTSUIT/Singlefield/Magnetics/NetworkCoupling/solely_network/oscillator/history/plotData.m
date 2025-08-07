clear, clc;


e1 = extractvoltage('oscillator-elecNetworkPotential-node-1-E1.hist');
e2 = extractvoltage('oscillator-elecNetworkPotential-node-2-E1.hist');
e3 = extractvoltage('oscillator-elecNetworkPotential-node-3-E4.hist');
e4 = extractvoltage('oscillator-elecNetworkPotential-node-4-E2.hist');
t = e1(:,1);

% resonance behaviour (T_res) via inductor voltage
figure(1), clf, hold on, grid minor
plot(t,e1(:,2)-e4(:,2))


% stationary voltage (e1 via I and Rint)
figure(2), clf, hold on, grid minor
plot(t,e1(:,2))


function data = extractvoltage(fname)

fid = fopen(fname);
data = textscan(fid,'%f  %f','HeaderLines',3) ;
data = cell2mat(data);
fclose(fid);

end

