clear all;
close all;

L = 10;
lambda = 13;
c_p = 40;
rho = 10;
c = sqrt(lambda/(rho*c_p));
C0 = 0; % Dirichlet left
W = 1/(rho*c_p); % source
u = 0.01;

% for finite domain
T = @(x) C0 + W*x/u + W*c^2/u^2*(exp(-u*L/c^2)-exp((x-L)*u/c^2));

% for infinite domain
%T = @(x) C0 + W*x/u;


xT = 0:0.01:L;

plot(xT,T(xT));
