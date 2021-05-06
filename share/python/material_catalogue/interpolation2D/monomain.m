function [Coeff11,Coeff22,Coeff33,Coeff23,Coeff13,Coeff12,params,E11,E22,E33,E23,E13,E12] = monomain(inputfile,outputfile)
%Einlesen_der Materialtensoren aus dem Materialkatalog
% Rotation des Materialkatalogs
%angle = 0;%pi/4;
%list = detailed_stats_10;
%list = rotate_list(detailed_stats_10,angle);
list = load(inputfile);

nparam = 1;

m = list(1,1)-1;
assert(m > 1);
params{1} = list(2:end,1);
da = list(3,1) - list(2,1);

n = max(0, list(1,2)-1);
if n > 1
    params{2} = list(2:end,2);
    db = list(3,2) - list(2,2);
else
    db = 1;
end
o = max(0, list(1,3)-1);
if o > 1
    params{3} = list(2:end,3);
    dc = list(3,3) - list(2,3);
else
    dc = 1;
end

E11 = zeros(m+1,n+1,o+1);
E22 = zeros(m+1,n+1,o+1);
E33 = zeros(m+1,n+1,o+1);
E23 = zeros(m+1,n+1,o+1);
E13 = zeros(m+1,n+1,o+1);
E12 = zeros(m+1,n+1,o+1);
for i=2:size(list,1)
%     idx = sub2ind(size(E11), list(i,1)+1, min(1,n)*list(i,2)+1, min(1,o)*list(i,3)+1);
    suba =            round((list(i,1)-list(2,1)) / da) + 1;
    subb = min(1,n) * round((list(i,2)-list(2,2)) / db) + 1;
    subc = min(1,o) * round((list(i,3)-list(2,3)) / dc) + 1;
    idx = sub2ind(size(E11), suba, subb, subc);
    E11(idx) = list(i,1+nparam);
    E22(idx) = list(i,2+nparam);
    E33(idx) = list(i,3+nparam);
    E23(idx) = list(i,4+nparam);
    E13(idx) = list(i,5+nparam);
    E12(idx) = list(i,6+nparam);
end


%optional: only needed for penalization
interpolation_func = @monocubic_offline;
deriv{1} = [];
deriv{2} = [];
deriv{3} = [];
if nparam > 1
    interpolation_func = @bicubic_offline;
    deriv{4} = [];
    deriv{5} = [];
    deriv{6} = [];
end
if nparam > 2
    interpolation_func = @tricubic_offline;
    deriv{7} = [];
    deriv{8} = [];
    deriv{9} = [];
end
% Coefficients for bicubic interpolation polynomial
[Coeff11] = interpolation_func(params, E11, deriv);
[Coeff22] = interpolation_func(params, E22, deriv);
[Coeff33] = interpolation_func(params, E33, deriv);
[Coeff23] = interpolation_func(params, E23, deriv);
[Coeff13] = interpolation_func(params, E13, deriv);
[Coeff12] = interpolation_func(params, E12, deriv);

%test(a,b,E11);
write_to_xml(outputfile, params, Coeff11, Coeff22, Coeff33, Coeff23, Coeff13, Coeff12, []);

