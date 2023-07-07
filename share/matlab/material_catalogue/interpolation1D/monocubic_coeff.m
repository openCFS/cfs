function [c] = monocubic_coeff(E, dEda, da, boundary)
% Berechnung der Interpolationskoeffizienten für ein bestimmtes Intervall
% Polynom ist definiert auf dem Referenzintervall [0,1]
if nargin < 4 || boundary == 0
    A=[1 0 0 0
       1 1 1 1
       0 1 0 0
       0 1 2 3];

    x = [E(1), E(2), dEda(1)*da, dEda(2)*da];
    c = A\x';
else
    if boundary == 1
        A=[1 0 0
           1 1 1
           0 1 0];
        x = [E(1), E(2), dEda(1)*da];
        c = A\x';
        c(4) = 0;
    else
        % Quadratische Interpolation kann zu Unterschwingern fuehren
%         A=[1 0 0
%            1 1 1
%            0 1 2];
%         x = [E(1), E(2), dEda(2)*da];
        A=[1 0 0 0
           1 1 1 1
           0 1 0 0
           0 1 2 3];
        x = [E(1), E(2), 0, dEda(2)*da];
        c = A\x';
    end
end
end
