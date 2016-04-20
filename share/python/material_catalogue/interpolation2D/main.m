function [Coeff11,Coeff12,Coeff22,Coeff33,a,b,E11,E12,E22,E33] = main(inputfile,outputfile,opt,outputfile2)
%Einlesen_der Materialtensoren aus dem Materialkatalog
% Rotation des Materialkatalogs
%angle = 0;%pi/4;
%list = detailed_stats_10;
%list = rotate_list(detailed_stats_10,angle);
list = load(inputfile);
m = list(1,1);
n = list(1,2);
da = 1/m;
db = 1/n;
a = [0:da:1];
b = [0:db:1];
E11 = zeros(m+1,n+1);
for i=2:size(list,1)
    E11(list(i,1)+1,list(i,2)+1) = list(i,3);
end
E12 = zeros(m+1,n+1);
for i=2:size(list,1)
    E12(list(i,1)+1,list(i,2)+1) = list(i,4);
end

E22 = zeros(m+1,n+1);
for i=2:size(list,1)
    E22(list(i,1)+1,list(i,2)+1) = list(i,5);
end
E33 = zeros(m+1,n+1);
for i=2:size(list,1)
    E33(list(i,1)+1,list(i,2)+1) = list(i,6);
end
%optional: only needed for penalization
deriv_a = [];
deriv_b = [];
deriv_a2 = [];
deriv_b2 = [];
% Coefficients for bicubic interpolation polynomial
[Coeff11] = bicubic_offline(a,b,E11, deriv_a,deriv_b,deriv_a2,deriv_b2);
[Coeff12] = bicubic_offline(a,b,E12, deriv_a,deriv_b,deriv_a2,deriv_b2);
[Coeff22] = bicubic_offline(a,b,E22, deriv_a,deriv_b,deriv_a2,deriv_b2);
[Coeff33] = bicubic_offline(a,b,E33, deriv_a,deriv_b,deriv_a2,deriv_b2);
%test(a,b,E11);
write_to_xml(outputfile,m,n,a,b,Coeff11,Coeff12,Coeff22,Coeff33);

% Calculate penalization material catalogue in 2D for interval
% [0,da]x[0,db]
% penalization: scale material tensor entry E(da,db) by function (x/a(2))^3 * (y/b(2))^3
if opt
    %list2 = load(inputfile2);
    m_p = list(1,1);
    n_p = list(1,2);
    da_p = a(2)/m_p;
    db_p = b(2)/n_p;
    %a(2) and b(2) is e.g. 0.1 if material catalogue is [0:0.1:1]
    a_p = [0:da_p:1];
    b_p = [0:db_p:1];
    E11_p = zeros(m*m_p+1,n*n_p+1);
    for i=1:m_p+1
       for j=1:n_p+1
          E11_p(i,j) = E11(2,2)*(a_p(i)/a(2))^3*(b_p(j)/b(2))^3;
          % check if penalized value is lower than void tensor
          if E11_p(i,j) < E11(1,1)
             E11_p(i,j) = E11(1,1); 
          end
       end
    end
    E12_p = zeros(m_p+1,n_p+1);
    for i=1:m_p+1
       for j=1:n_p+1
          E12_p(i,j) = E12(2,2)*(a_p(i)/a(2))^3*(b_p(j)/b(2))^3;
          if E12_p(i,j) < E12(1,1)
             E12_p(i,j) = E12(1,1); 
          end
       end
    end
    E22_p = zeros(m_p+1,n_p+1);
    for i=1:m_p+1
       for j=1:n_p+1
          E22_p(i,j) = E22(2,2)*(a_p(i)/a(2))^3*(b_p(j)/b(2))^3;
          if E22_p(i,j) < E22(1,1)
             E22_p(i,j) = E22(1,1); 
          end
       end
    end
    E33_p = zeros(m_p+1,n_p+1);
    for i=1:m_p+1
       for j=1:n_p+1
          E33_p(i,j) = E33(2,2)*(a_p(i)/a(2))^3*(b_p(j)/b(2))^3;
          if E33_p(i,j) < E33(1,1)
             E33_p(i,j) = E33(1,1); 
          end
       end
    end
%     for i=2:size(list2,1)
%         E11_p(list2(i,1)*m_p+1,list2(i,2)*n_p+1) = list2(i,3);
%     end
%     E12_p = zeros(m_p+1,n_p+1);
%     for i=2:size(list,1)
%         E12_p(list2(i,1)*m_p+1,list2(i,2)*n_p+1) = list2(i,4);
%     end
% 
%     E22_p = zeros(m_p+1,n_p+1);
%     for i=2:size(list2,1)
%         E22_p(list2(i,1)*m_p+1,list2(i,2)*n_p+1) = list2(i,5);
%     end
%     E33_p = zeros(m_p+1,n_p+1);
%     for i=2:size(list2,1)
%         E33_p(list2(i,1)*m_p+1,list2(i,2)*n_p+1) = list2(i,6);
%     end
    % Calculate derivatives from material catalogue [0,1] at end points 
    % and calculate penalization interpolation coefficients Coeff_p
    deriv_a = zeros(n_p,1);
    deriv_b = zeros(n_p,1);
    for i = 1:n_p+1
        [deriv_a(i),deriv_b(i)] = bicubic_deriv(Coeff11,a,b,a_p(m),b_p(i));
    end
    for i = 1:m_p+1
        [deriv_a2(i),deriv_b2(i)] = bicubic_deriv(Coeff11,a,b,a_p(i),b_p(m)); 
    end
    [Coeff11_p] = bicubic_offline(a_p,b_p,E11_p, deriv_a,deriv_b,deriv_a2,deriv_b2);
    for i = 1:n_p+1
        [deriv_a(i),deriv_b(i)] = bicubic_deriv(Coeff12,a,b,a_p(m),b_p(i)); 
    end
    for i = 1:m_p+1
        [deriv_a2(i),deriv_b2(i)] = bicubic_deriv(Coeff12,a,b,a_p(i),b_p(m)); 
    end
    [Coeff12_p] = bicubic_offline(a_p,b_p,E12_p, deriv_a,deriv_b,deriv_a2,deriv_b2);
    for i = 1:n_p+1
        [deriv_a(i),deriv_b(i)] = bicubic_deriv(Coeff22,a,b,a_p(m),b_p(i)); 
    end
    for i = 1:m_p+1
        [deriv_a2(i),deriv_b2(i)] = bicubic_deriv(Coeff22,a,b,a_p(i),b_p(m)); 
    end
    [Coeff22_p] = bicubic_offline(a_p,b_p,E22_p, deriv_a,deriv_b,deriv_a2,deriv_b2);
    for i = 1:n_p
        [deriv_a(i),deriv_b(i)] = bicubic_deriv(Coeff33,a,b,a_p(m),b_p(i)); 
    end
    for i = 1:m_p+1
        [deriv_a2(i),deriv_b2(i)] = bicubic_deriv(Coeff11,a,b,a_p(i),b_p(m)); 
    end
    [Coeff33_p] = bicubic_offline(a_p,b_p,E33_p,deriv_a,deriv_b,deriv_a2,deriv_b2);
    write_to_xml(outputfile2,m_p,n_p,a_p,b_p,Coeff11_p,Coeff12_p,Coeff22_p,Coeff33_p);
end
end
