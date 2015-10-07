function res = rotatevoigt(Tensor, theta)
% ROTATEHILLMANDEL rotates a Tensor given in Hill-Mandel notation by theta

assert(theta >= 0 && theta <= pi);

res = hillmandeltovoigt(rotatehillmandel(voigttohillmandel(Tensor),theta));