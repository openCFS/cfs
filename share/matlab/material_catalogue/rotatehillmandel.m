function res = rotatehillmandel(Tensor, theta)
% ROTATEHILLMANDEL rotates a Tensor given in Hill-Mandel notation by theta
% in mathematical positiv direction

assert(theta >= -pi && theta <= pi);

theta = -theta;

Q = zeros(3,3);

Q(1,1) = cos(theta)^2;
Q(1,2) = sin(theta)^2;
Q(1,3) = -1 * sqrt(2) / 2 * sin(2 * theta);
Q(2,1) = Q(1,2);
Q(2,2) = Q(1,1);
Q(2,3) = -Q(1,3);
Q(3,1) = Q(2,3);
Q(3,2) = Q(1,3);
Q(3,3) = cos(2 * theta);

res = Q' * Tensor * Q;
