function [Coeff] = bicubic_offline(a,b,E,E_number)
%% Vorberechnung der Koeffizienten des Interpolationspolynom für die verschiedenen Intervalle
m = length(a);
n = length(b);
[dEda, dEdb, dEdadb] = bicubic_partialderiv(a,b,E,E_number);
Coeff = zeros((m-1)*(n-1),16);
for j=1:m-1
    for k=1:n-1
        al = a(j);
        au = a(j+1);
        bl = b(k);
        bu = b(k+1);
        da=au-al;
        db=bu-bl;
        % Required Data in the corners of the chosen element
        Eint(1) = E(j,k);
        Eint(2) = E(j+1,k);
        Eint(3) = E(j+1,k+1);
        Eint(4) = E(j,k+1);
        Eda(1) = dEda(j,k);
        Eda(2) = dEda(j+1,k);
        Eda(3) = dEda(j+1,k+1);
        Eda(4) = dEda(j,k+1);
        Edb(1) = dEdb(j,k);
        Edb(2) = dEdb(j+1,k);
        Edb(3) = dEdb(j+1,k+1);
        Edb(4) = dEdb(j,k+1);
        Edadb(1) = dEdadb(j,k);
        Edadb(2) = dEdadb(j+1,k);
        Edadb(3) = dEdadb(j+1,k+1);
        Edadb(4) = dEdadb(j,k+1);
        Coeff((n-1)*(j-1)+k,:) = bicubic_coeff(Eint, Eda, Edb, Edadb, da, db);
    end
end
end
