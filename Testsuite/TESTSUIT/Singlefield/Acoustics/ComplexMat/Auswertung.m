%compute absorption coefficient form FE calculations
filename1 = 'history/impRohr1-acouPressure-node-12-mic1.hist';
filename2 = 'history/impRohr1-acouPressure-node-56-mic2.hist';

x2 = 20e-3; %distance of mic2 from sample
s = 40e-3;  %distance nnetween mic1 and mic2

c_0 = 340;

delimiterIn = ' ';
headerlinesIn = 3;
dataP1 = importdata(filename1,delimiterIn,headerlinesIn); 
dataP2 = importdata(filename2,delimiterIn,headerlinesIn); 

fsim=dataP1.data(:,1);
p1 = dataP1.data(:,2) + i*dataP1.data(:,3);;
p2 = dataP2.data(:,2) + i*dataP2.data(:,3);

H_21 = p2./p1;
k = 2*pi*fsim./c_0;

rSim = exp(i*2*k*(x2+s)).*( (H_21 -exp(-i*k*s))./( exp(i*k*s) - H_21));

%r2 = xp(i*2*k*
alphaSim= 1 -abs(rSim).^2;

plot(fsim,alphaSim),grid
