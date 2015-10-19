function [c] = tricubic_coeff(E, dEda, dEdb, dEdc, dEdadb,dEdadc,dEdbdc,dEdadbdc, da, db,dc)

%Given arrays containing the function, gra-
%dients, and cross derivative at the eigth grid points of a cubic grid cell
%and given grid size da,db and dc
%Calculates the coefficients for the tricubic interpolation

dadb=1.;
dadc = 1.;
dbdc = 1.;
dadbdc = 1.;
%   dadb = da*db;
%   dadc = da*dc;
%   dbdc = db*dc;
%   dadbdc = da*db*dc;
x = zeros(64,1);
A = getA();
x(1:8) = E;
x(9:16) = dEda;
x(17:24) = dEdb;
x(25:32) = dEdc;
x(33:40) = dEdadb;
x(41:48) = dEdadc;
x(49:56) = dEdbdc;
x(57:64) = dEdadbdc;
c = A*x;
end

