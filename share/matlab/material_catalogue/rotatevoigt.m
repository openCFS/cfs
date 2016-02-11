function res = rotatevoigt(Tensor, theta)
% ROTATEHILLMANDEL rotates a Tensor given in Hill-Mandel notation by theta

res = hillmandeltovoigt(rotatehillmandel(voigttohillmandel(Tensor),theta));