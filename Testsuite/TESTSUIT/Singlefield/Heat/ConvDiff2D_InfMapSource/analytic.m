clear all;
close all;

L = 10;
c = 1;
C0 = 0; % Dirichlet left
W = 1; % source
u = 7;

% for finite domain
%T = @(x) C0 + W*x/u + W*c^2/u^2*(exp(-u*L/c^2)-exp((x-L)*u/c^2));

% for infinite domain
T = @(x) C0 + W*x/u;


xT = 0:0.01:L;

plot(xT,T(xT));
