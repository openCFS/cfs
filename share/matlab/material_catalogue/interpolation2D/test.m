function test(a,b,E)
%Test bicubic interpolation if it interpolates data points
%In addition the derivative ist checkt by finite differences
params{1} = a;
params{2} = a;
[Coeff] = bicubic_offline(params,E);
for i=1:length(a)
    for j=1:length(b)
        [res1] = bicubic_int(Coeff,a,b,a(i),b(j));
        error = abs(E(i,j)-res1);
        if error > 1.e-6
            display(sprintf('error: %f for index: %d , %d',error,i,j));
        end
    end
end
% derivative test
% epsi = 1.e-10;
% maxer1 = -1;
% maxer2 = -1;
% for i=1:length(a)-1
%     for j=1:length(b)-1
%         [fd1e] = bicubic_int(Coeff,a,b,a(i)+epsi,b(j));
%         [fd2e] = bicubic_int(Coeff,a,b,a(i),b(j)+epsi);
%         [fd] = bicubic_int(Coeff,a,b,a(i),b(j));
%         fd1 = (fd1e-fd)/epsi;
%         fd2 = (fd2e-fd)/epsi;
%         [deriv1,deriv2] = bicubic_deriv(Coeff,a,b,0,0);
%         error1 = (abs(deriv1-fd1));
%         maxer1 = max(error1,maxer1);
%         if error1 > 1.e-6
%             display(sprintf('error: %f, deriv 1 for index: %d , %d',error1,i,j));
%         end
%         error2 = abs(deriv2-fd2);
%         maxer2 = max(error2,maxer2);
%         if  error2 > 1.e-6
%             display(sprintf('error: %f deriv 2 for index: %d , %d',error2, i,j));
%         end
%     end
% end
% display(sprintf('maxerror: 1: %f 2:%f ',maxer1,maxer2));
