function [coeffs,params,E] = main(inputfile,outputfile,opt,outputfile2)
%Einlesen_der Materialtensoren aus dem Materialkatalog
% Rotation des Materialkatalogs
%angle = 0;%pi/4;
%list = detailed_stats_10;
%list = rotate_list(detailed_stats_10,angle);
data = load(inputfile);

[params, nparamintervals, dparams] = getSampling(data);

nparam = sum(nparamintervals > 0);

% normalize parameters
for d=1:nparam
    params{d} = params{d} / params{d}(end);
    if max(data(:,d)) ~= 1
        data(:,d) = data(:,d)/max(data(:,d));
    end
end

if nparam == 1
    assert(nparamintervals(1) > 1);
    assert(nparamintervals(2) == 0);
    assert(nparamintervals(3) == 0);
elseif nparam == 2
    assert(nparamintervals(1) > 1);
    assert(nparamintervals(2) > 1);
    assert(nparamintervals(3) == 0);
else
    assert(nparamintervals(1) > 1);
    assert(nparamintervals(2) > 1);
    assert(nparamintervals(3) > 1);
end

regular = true;
if any( isnan(dparams) )
    regular = false;
end

E11 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E22 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E33 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E23 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E12 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E13 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E14 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E15 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E16 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E24 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E25 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E26 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E34 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E35 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E36 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E44 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E45 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E46 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E55 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E56 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
E66 = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);
Vol = zeros(nparamintervals(1)+1, nparamintervals(2)+1, nparamintervals(3)+1);

for i=2:size(data,1)
%     idx = sub2ind(size(E11), list(i,1)+1, min(1,n)*list(i,2)+1, min(1,o)*list(i,3)+1);
    sub = ones(3,1);
    for d=1:nparam
        if regular
            sub(d) = min(1, nparamintervals(d)) * round((data(i,d)-data(2,d)) / dparams(d)) + 1;
        else
            sub(d) = find(params{d} == data(i,d));
        end
    end
    idx = sub2ind(size(E11), sub(1), sub(2), sub(3));
    if size(data,2) == 4+nparam || size(data,2) == 5+nparam
        E11(idx) = data(i, 1+nparam);
        E12(idx) = data(i, 2+nparam);
        E22(idx) = data(i, 3+nparam);
        E33(idx) = data(i, 4+nparam);
        if size(data,2) == 5+nparam
            Vol(idx) = data(i, 5+nparam);
        end
    elseif size(data,2) == 6+nparam || size(data,2) == 7+nparam
        E11(idx) = data(i, 1+nparam);
        E22(idx) = data(i, 2+nparam);
        E33(idx) = data(i, 3+nparam);
        E23(idx) = data(i, 4+nparam);
        E13(idx) = data(i, 5+nparam);
        E12(idx) = data(i, 6+nparam);
        if size(data,2) == 7+nparam
            Vol(idx) = data(i, 7+nparam);
        end
    else
        assert(size(data,2) == 21+nparam || size(data,2) == 22+nparam);
        E11(idx) = data(i, 1+nparam);
        E12(idx) = data(i, 2+nparam);
        E13(idx) = data(i, 3+nparam);
        E14(idx) = data(i, 4+nparam);
        E15(idx) = data(i, 5+nparam);
        E16(idx) = data(i, 6+nparam);
        E22(idx) = data(i, 7+nparam);
        E23(idx) = data(i, 8+nparam);
        E24(idx) = data(i, 9+nparam);
        E25(idx) = data(i,10+nparam);
        E26(idx) = data(i,11+nparam);
        E33(idx) = data(i,12+nparam);
        E34(idx) = data(i,13+nparam);
        E35(idx) = data(i,14+nparam);
        E36(idx) = data(i,15+nparam);
        E44(idx) = data(i,16+nparam);
        E45(idx) = data(i,17+nparam);
        E46(idx) = data(i,18+nparam);
        E55(idx) = data(i,19+nparam);
        E56(idx) = data(i,20+nparam);
        E66(idx) = data(i,21+nparam);
        if size(data,2) == 22+nparam
            Vol(idx) = data(i, 22+nparam);
        end
    end
end

E = struct(...
'E11', E11,...
'E12', E12,...
'E13', E13,...
'E14', E14,...
'E15', E15,...
'E16', E16,...
'E22', E22,...
'E23', E23,...
'E24', E24,...
'E25', E25,...
'E26', E26,...
'E33', E33,...
'E34', E34,...
'E35', E35,...
'E36', E36,...
'E44', E44,...
'E45', E45,...
'E46', E46,...
'E55', E55,...
'E56', E56,...
'E66', E66,...
'Vol', Vol...
);

interpolation_func = @monocubic_offline;
% optional: deriv is only needed for penalization
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

temp = zeros(prod(nparamintervals(nparamintervals > 0)), 4^nparam, 22);

% Coefficients for interpolation polynomial
names = fieldnames(E);
for k = 1:length(names)
    [temp(:,:,k)] = interpolation_func(params, E.(names{k}), deriv);
end

coeffs = struct(...
'coeff11', temp(:,:,1),...
'coeff12', temp(:,:,2),...
'coeff13', temp(:,:,3),...
'coeff14', temp(:,:,4),...
'coeff15', temp(:,:,5),...
'coeff16', temp(:,:,6),...
'coeff22', temp(:,:,7),...
'coeff23', temp(:,:,8),...
'coeff24', temp(:,:,9),...
'coeff25', temp(:,:,10),...
'coeff26', temp(:,:,11),...
'coeff33', temp(:,:,12),...
'coeff34', temp(:,:,13),...
'coeff35', temp(:,:,14),...
'coeff36', temp(:,:,15),...
'coeff44', temp(:,:,16),...
'coeff45', temp(:,:,17),...
'coeff46', temp(:,:,18),...
'coeff55', temp(:,:,19),...
'coeff56', temp(:,:,20),...
'coeff66', temp(:,:,21)...
);

if norm(temp(:,:,22)) ~= 0
    coeffs = setfield(coeffs, 'coeffvol',temp(:,:,22));
end

test(params{1},params{2},E11);
write_to_xml(outputfile, params, coeffs);

% % Calculate penalization material catalogue in 2D for interval
% % [0,da]x[0,db]
% % penalization: scale material tensor entry E(da,db) by function (x/a(2))^3 * (y/b(2))^3
% if opt
%     %list2 = load(inputfile2);
%     m_p = list(1,1);
%     n_p = list(1,2);
%     da_p = a(2)/m_p;
%     db_p = b(2)/n_p;
%     %a(2) and b(2) is e.g. 0.1 if material catalogue is [0:0.1:1]
%     a_p = [0:da_p:1];
%     b_p = [0:db_p:1];
%     E11_p = zeros(m*m_p+1,n*n_p+1);
%     for i=1:m_p+1
%        for j=1:n_p+1
%           E11_p(i,j) = E11(2,2)*(a_p(i)/a(2))^3*(b_p(j)/b(2))^3;
%           % check if penalized value is lower than void tensor
%           if E11_p(i,j) < E11(1,1)
%              E11_p(i,j) = E11(1,1); 
%           end
%        end
%     end
%     E12_p = zeros(m_p+1,n_p+1);
%     for i=1:m_p+1
%        for j=1:n_p+1
%           E12_p(i,j) = E12(2,2)*(a_p(i)/a(2))^3*(b_p(j)/b(2))^3;
%           if E12_p(i,j) < E12(1,1)
%              E12_p(i,j) = E12(1,1); 
%           end
%        end
%     end
%     E22_p = zeros(m_p+1,n_p+1);
%     for i=1:m_p+1
%        for j=1:n_p+1
%           E22_p(i,j) = E22(2,2)*(a_p(i)/a(2))^3*(b_p(j)/b(2))^3;
%           if E22_p(i,j) < E22(1,1)
%              E22_p(i,j) = E22(1,1); 
%           end
%        end
%     end
%     E33_p = zeros(m_p+1,n_p+1);
%     for i=1:m_p+1
%        for j=1:n_p+1
%           E33_p(i,j) = E33(2,2)*(a_p(i)/a(2))^3*(b_p(j)/b(2))^3;
%           if E33_p(i,j) < E33(1,1)
%              E33_p(i,j) = E33(1,1); 
%           end
%        end
%     end
% %     for i=2:size(list2,1)
% %         E11_p(list2(i,1)*m_p+1,list2(i,2)*n_p+1) = list2(i,3);
% %     end
% %     E12_p = zeros(m_p+1,n_p+1);
% %     for i=2:size(list,1)
% %         E12_p(list2(i,1)*m_p+1,list2(i,2)*n_p+1) = list2(i,4);
% %     end
% % 
% %     E22_p = zeros(m_p+1,n_p+1);
% %     for i=2:size(list2,1)
% %         E22_p(list2(i,1)*m_p+1,list2(i,2)*n_p+1) = list2(i,5);
% %     end
% %     E33_p = zeros(m_p+1,n_p+1);
% %     for i=2:size(list2,1)
% %         E33_p(list2(i,1)*m_p+1,list2(i,2)*n_p+1) = list2(i,6);
% %     end
%     % Calculate derivatives from material catalogue [0,1] at end points 
%     % and calculate penalization interpolation coefficients Coeff_p
%     deriv_a = zeros(n_p,1);
%     deriv_b = zeros(n_p,1);
%     for i = 1:n_p+1
%         [deriv_a(i),deriv_b(i)] = bicubic_deriv(Coeff11,a,b,a_p(m),b_p(i));
%     end
%     for i = 1:m_p+1
%         [deriv_a2(i),deriv_b2(i)] = bicubic_deriv(Coeff11,a,b,a_p(i),b_p(m)); 
%     end
%     [Coeff11_p] = bicubic_offline(a_p,b_p,E11_p, deriv_a,deriv_b,deriv_a2,deriv_b2);
%     for i = 1:n_p+1
%         [deriv_a(i),deriv_b(i)] = bicubic_deriv(Coeff12,a,b,a_p(m),b_p(i)); 
%     end
%     for i = 1:m_p+1
%         [deriv_a2(i),deriv_b2(i)] = bicubic_deriv(Coeff12,a,b,a_p(i),b_p(m)); 
%     end
%     [Coeff12_p] = bicubic_offline(a_p,b_p,E12_p, deriv_a,deriv_b,deriv_a2,deriv_b2);
%     for i = 1:n_p+1
%         [deriv_a(i),deriv_b(i)] = bicubic_deriv(Coeff22,a,b,a_p(m),b_p(i)); 
%     end
%     for i = 1:m_p+1
%         [deriv_a2(i),deriv_b2(i)] = bicubic_deriv(Coeff22,a,b,a_p(i),b_p(m)); 
%     end
%     [Coeff22_p] = bicubic_offline(a_p,b_p,E22_p, deriv_a,deriv_b,deriv_a2,deriv_b2);
%     for i = 1:n_p
%         [deriv_a(i),deriv_b(i)] = bicubic_deriv(Coeff33,a,b,a_p(m),b_p(i)); 
%     end
%     for i = 1:m_p+1
%         [deriv_a2(i),deriv_b2(i)] = bicubic_deriv(Coeff11,a,b,a_p(i),b_p(m)); 
%     end
%     [Coeff33_p] = bicubic_offline(a_p,b_p,E33_p,deriv_a,deriv_b,deriv_a2,deriv_b2);
%     write_to_xml(outputfile2,m_p,n_p,a_p,b_p,Coeff11_p,Coeff12_p,Coeff22_p,Coeff33_p);
% end
end


function [params, nparamintervals, dparams] = getSampling(list)

nparamintervals = zeros(3,1);
dparams = zeros(3,1);

for d=1:3
    nparamintervals(d) = max(0, list(1,d));
    if nparamintervals(d) > 1
        params{d} = unique(list(2:end,d));
        dparams(d) = list(3,d) - list(2,d);

        equi = params{d}(1):dparams(d):params{d}(end);
        if numel(equi) ~= numel(params{d}) || any( abs(equi - params{d}') > 1e-6*dparams(d))
            dparams(d) = NaN;
        end
    else
        dparams(d) = 1;
        params{d} = [];
    end
end
% in detailed_stats the numbers in the first line are the numbers of
% intervals for each parameter
% in other catalogue files the numbers correspond to the numbers of points
if size(list,1) - 1 < prod( nparamintervals(nparamintervals > 0)+1 )
    nparamintervals = max(0, nparamintervals - 1);
end
end
