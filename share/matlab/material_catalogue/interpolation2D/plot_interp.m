% [Coeff11,Coeff12,Coeff22,Coeff33,a,b,E11,E12,E22,E33] = main('coeff.xml');
% load detailed_stats_10
% list = detailed_stats_10;
% m = list(1,1);
% n = list(1,2);
% da = 1/m;
% db = 1/n;
% a2 = a(1):da:a(end);
% b2 = b(1):db:b(end);
% E11_grid = zeros(m+1,n+1);
% for i=2:size(list,1)
%     E11_grid(list(i,1)+1,list(i,2)+1) = list(i,3);
% end
% E12_grid = zeros(m+1,n+1);
% for i=2:size(list,1)
%     E12_grid(list(i,1)+1,list(i,2)+1) = list(i,4);
% end
% E22_grid = zeros(m+1,n+1);
% for i=2:size(list,1)
%     E22_grid(list(i,1)+1,list(i,2)+1) = list(i,5);
% end
% E33_grid = zeros(m+1,n+1);
% for i=2:size(list,1)
%     E33_grid(list(i,1)+1,list(i,2)+1) = list(i,6);
% end
[XX,YY]=meshgrid(a(1):1:a(end),b(1):1:b(end));
ZZ = zeros(size(XX));
for ii=1:size(XX,1)
    for jj=1:size(XX,2)
        ZZ(ii,jj) = bicubic_int(Coeff11,a,b,XX(ii,jj),YY(ii,jj));
    end
end
surf(XX,YY,ZZ);
hold on;
% [A,B]=meshgrid(a(1):da:a(end),b(1):db:b(end));
% plot3(A',B',E11_grid,'*');
xlabel('stiff1')
ylabel('stiff2')
hold off;

% figure;
% for ii=1:size(XX,1)
%     for jj=1:size(XX,2)
%         ZZ(ii,jj) = bicubic_int(Coeff12,a,b,XX(ii,jj),YY(ii,jj));
%     end
% end
% surf(XX,YY,ZZ);
% hold on;
% [A,B]=meshgrid(a(1):da:a(end),b(1):db:b(end));
% plot3(A',B',E12_grid,'*');
% xlabel('stiff1')
% ylabel('stiff2')
% hold off;
% 
% figure;
% for ii=1:size(XX,1)
%     for jj=1:size(XX,2)
%         ZZ(ii,jj) = bicubic_int(Coeff22,a,b,XX(ii,jj),YY(ii,jj));
%     end
% end
% surf(XX,YY,ZZ);
% hold on;
% [A,B]=meshgrid(a(1):da:a(end),b(1):db:b(end));
% plot3(A',B',E22_grid,'*');
% xlabel('stiff1')
% ylabel('stiff2')
% hold off;
% 
% figure;
% for ii=1:size(XX,1)
%     for jj=1:size(XX,2)
%         ZZ(ii,jj) = bicubic_int(Coeff33,a,b,XX(ii,jj),YY(ii,jj));
%     end
% end
% surf(XX,YY,ZZ);
% hold on;
% [A,B]=meshgrid(a(1):da:a(end),b(1):db:b(end));
% plot3(A',B',E33_grid,'*');
% xlabel('stiff1')
% ylabel('stiff2')
% hold off;
% 
