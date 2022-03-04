function [Coeff] = monocubic_offline(params, E, deriv)
%% Vorberechnung der Koeffizienten des Interpolationspolynom für die verschiedenen Intervalle
a = params{1};
deriv_a = deriv{1};
m = length(a);

dEda = monocubic_partialderiv(a, E, deriv_a);

Coeff = zeros(m-1,4);
for j=1:m-1
    al = a(j);
    au = a(j+1);
    da = au-al;
    % Required Data in the corners of the chosen element
    Eint(1) = E(j);
    Eint(2) = E(j+1);
    Eda(1) = dEda(j);
    Eda(2) = dEda(j+1);
    if j == 1
        Coeff(j,:) = monocubic_coeff(Eint, Eda, da, -1);
    elseif j == m-1
        Coeff(j,:) = monocubic_coeff(Eint, Eda, da, 1);
    else
        Coeff(j,:) = monocubic_coeff(Eint, Eda, da);
    end
end
end
