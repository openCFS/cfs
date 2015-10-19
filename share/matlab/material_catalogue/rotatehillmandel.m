function res = rotatehillmandel(Tensor, theta)
% ROTATEHILLMANDEL rotates a Tensor given in Hill-Mandel notation by theta
% Rotiert im Uhrzeigersinn!

assert(theta >= 0 && theta <= pi);

res = zeros(size(Tensor));

Q = zeros(3,3);

Q(1,1) = cos(theta)^2;
Q(1,2) = sin(theta)^2;
Q(1,3) = -1 * sqrt(2) / 2 * sin(2 * theta);
Q(2,1) = sin(theta)^2;
Q(2,2) = cos(theta)^2;
Q(2,3) = sqrt(2) / 2  * sin(2 * theta);
Q(3,1) = sqrt(2) / 2  * sin(2 * theta);
Q(3,2) = -1 * sqrt(2) / 2 * sin(2 * theta);
Q(3,3) = cos(2 * theta);

res = Q' * Tensor * Q;
