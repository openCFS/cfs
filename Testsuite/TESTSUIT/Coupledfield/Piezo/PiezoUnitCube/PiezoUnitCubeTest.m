clear all
clc
% Unit Cube Analytical Approach
% CE is Elasticity Tensor
% S is Strain Tensor
% e is Coupling Tensor
% E is Electric Field Intensity
% T is Stress Tensor

CE=[ 1.376e+11 7.627e+10 7.313e+10 0.000e+00 0.000e+00 0.000e+00 
     7.627e+10 1.376e+11 7.313e+10 0.000e+00 0.000e+00 0.000e+00 
     7.313e+10 7.313e+10 1.145e+11 0.000e+00 0.000e+00 0.000e+00 
     0.000e+00 0.000e+00 0.000e+00 2.564e+10 0.000e+00 0.000e+00
     0.000e+00 0.000e+00 0.000e+00 0.000e+00 2.564e+10 0.000e+00 
     0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 3.058e+10];
e=[  0.0000e-003  0.0000e-003  0.0000e-003  0.0000e-003 12.7174e+000 0.0000e-003 
     0.0000e-003  0.0000e-003  0.0000e-003 12.7174e+000  0.0000e-003 0.0000e-003 
    -5.1714e+000 -5.1714e+000 15.1005e+000  0.0000e-003  0.0000e-003 0.0000e-003];
E=[0;0;1.0e+6];

% Step1: Stress Free
T1=[0;0;0;0;0;0];
S1=CE^-1*(e'*E+T1)               %Constitutive Equation%

% Step2: Uniaxial Stress
T2=[0;0;10e+6;0;0;0];
S2=CE^-1*(e'*E+T2)               %Constitutive Equation%

% Step3: Uniaxial Strain
S3=[0;0;1e-4;0;0;0];
T3=CE*S3-e'*E                      %Constitutive Equation%

