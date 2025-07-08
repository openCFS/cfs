clf;
load('Jexact_matlab.ref.mat');
plot(angles,Jexact,'DisplayName','matlab','Marker','o')

hold on
cfs = dlmread('Jexact.ref.txt',' ', 1,0);
plot(angles,cfs(:,2),'DisplayName','cfs');

hold on 
model = dlmread('Jmodel.ref.txt',' ', 1,0);
plot(angles,model(:,2),'DisplayName','model')

legend