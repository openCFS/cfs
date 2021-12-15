function [c] = monocubic_coeff(E, dEda, da, atRightBoundary)
% Berechnung der Interpolationskoeffizienten für ein bestimmtes Intervall
% Polynom ist definiert auf dem Referenzintervall [0,1]
if nargin < 4 || ~atRightBoundary
    A=[1 0 0 0
       1 1 1 1
       0 1 0 0
       0 1 2 3];

    x = [E(1), E(2), dEda(1)*da, dEda(2)*da];
    c = A\x';
else
    A=[1 0 0
       1 1 1
       0 1 0];

    x = [E(1), E(2), dEda(1)*da];
    c = A\x';
    c(4) = 0;
end
end
