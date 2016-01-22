function [Coeff] = tricubic_offline(a,b,c,E)
m = length(a);
n = length(b);
o = length(c);
da = (a(2)-a(1));
db = (b(2)-b(1));
dc = (c(2)-c(1));
%[dEda, dEdb,dEdc, dEdadb,dEdadc,dEdbdc,dEdadbdc] = tricubic_partialderiv(a,b,c,E,da,db,dc);
Coeff = zeros((m-1)*(n-1)*(o-1),64);
for j=1:m-1
    for k=1:n-1
        for l=1:o-1
        al = a(j);
        au = a(j+1);
        bl = b(k);
        bu = b(k+1);
        cl = c(j);
        cu = c(j+1);
        da=au-al;
        db=bu-bl;
        dc = cu-cl;
        % Required Data in the corners of the chosen element
%         Eint(1) = E(j,k,l);
%         Eint(2) = E(j+1,k,l);
%         Eint(3) = E(j,k+1,l);
%         Eint(4) = E(j+1,k+1,l);
%         Eint(5) = E(j,k,l+1);
%         Eint(6) = E(j+1,k,l+1);
%         Eint(7) = E(j,k+1,l+1);
%         Eint(8) = E(j+1,k+1,l+1);
%         
%         Eda(1) = dEda(j,k,l);
%         Eda(2) = dEda(j+1,k,l);
%         Eda(3) = dEda(j,k+1,l);
%         Eda(4) = dEda(j+1,k+1,l);
%         Eda(5) = dEda(j,k,l+1);
%         Eda(6) = dEda(j+1,k,l+1);
%         Eda(7) = dEda(j,k+1,l+1);
%         Eda(8) = dEda(j+1,k+1,l+1);
%         
%         Edb(1) = dEdb(j,k,l);
%         Edb(2) = dEdb(j+1,k,l);
%         Edb(3) = dEdb(j,k+1,l);
%         Edb(4) = dEdb(j+1,k+1,l);
%         Edb(5) = dEdb(j,k,l+1);
%         Edb(6) = dEdb(j+1,k,l+1);
%         Edb(7) = dEdb(j,k+1,l+1);
%         Edb(8) = dEdb(j+1,k+1,l+1);
%         
%         Edc(1) = dEdc(j,k,l);
%         Edc(2) = dEdc(j+1,k,l);
%         Edc(3) = dEdc(j,k+1,l);
%         Edc(4) = dEdc(j+1,k+1,l);
%         Edc(5) = dEdc(j,k,l+1);
%         Edc(6) = dEdc(j+1,k,l+1);
%         Edc(7) = dEdc(j,k+1,l+1);
%         Edc(8) = dEdc(j+1,k+1,l+1);
%         
%         Edadb(1) = dEdadb(j,k,l);
%         Edadb(2) = dEdadb(j+1,k,l);
%         Edadb(3) = dEdadb(j,k+1,l);
%         Edadb(4) = dEdadb(j+1,k+1,l);
%         Edadb(5) = dEdadb(j,k,l+1);
%         Edadb(6) = dEdadb(j+1,k,l+1);
%         Edadb(7) = dEdadb(j,k+1,l+1);
%         Edadb(8) = dEdadb(j+1,k+1,l+1);
%         
%         Edadc(1) = dEdadc(j,k,l);
%         Edadc(2) = dEdadc(j+1,k,l);
%         Edadc(3) = dEdadc(j,k+1,l);
%         Edadc(4) = dEdadc(j+1,k+1,l);
%         Edadc(5) = dEdadc(j,k,l+1);
%         Edadc(6) = dEdadc(j+1,k,l+1);
%         Edadc(7) = dEdadc(j,k+1,l+1);
%         Edadc(8) = dEdadc(j+1,k+1,l+1);
%         
%         Edbdc(1) = dEdbdc(j,k,l);
%         Edbdc(2) = dEdbdc(j+1,k,l);
%         Edbdc(3) = dEdbdc(j,k+1,l);
%         Edbdc(4) = dEdbdc(j+1,k+1,l);
%         Edbdc(5) = dEdbdc(j,k,l+1);
%         Edbdc(6) = dEdbdc(j+1,k,l+1);
%         Edbdc(7) = dEdbdc(j,k+1,l+1);
%         Edbdc(8) = dEdbdc(j+1,k+1,l+1);
%         
%         Edadbdc(1) = dEdadbdc(j,k,l);
%         Edadbdc(2) = dEdadbdc(j+1,k,l);
%         Edadbdc(3) = dEdadbdc(j,k+1,l);
%         Edadbdc(4) = dEdadbdc(j+1,k+1,l);
%         Edadbdc(5) = dEdadbdc(j,k,l+1);
%         Edadbdc(6) = dEdadbdc(j+1,k,l+1);
%         Edadbdc(7) = dEdadbdc(j,k+1,l+1);
%         Edadbdc(8) = dEdadbdc(j+1,k+1,l+1);
        [Eint Eda Edb Edc Edadb Edadc Edbdc Edadbdc] = tricubic_parderiv(j,k,l,E,da,db,dc);  
        Coeff((n-1)*(o-1)*(j-1)+(o-1)*(k-1)+l,:) = tricubic_coeff(Eint, Eda, Edb,Edc, Edadb,Edadc,Edbdc,Edadbdc, da, db,dc);
        end
    end
end
end