function [dEda, dEdb, dEdc, dEdadb, dEdbdc, dEdadc,dEdadbdc] = tricubic_partialderiv(a,b,c,E,da,db,dc)
dEda = zeros(length(a),length(b),length(c));
dEdb = zeros(length(a),length(b),length(c));
dEdc = zeros(length(a),length(b),length(c));
dEdadb = zeros(length(a),length(b),length(c));
dEdadc = zeros(length(a),length(b),length(c));
dEdbdc = zeros(length(a),length(b),length(c));
dEdadbdc = zeros(length(a),length(b),length(c));

tmp = zeros(5,5,5);
tmp(2:end-1,2:end-1,2:end-1) = E;
E = tmp;
%Central Difference quotient for inner points
for j=2:length(a)+1
    for k=2:length(b)+1
        for l = 2:length(c)+1
            %if j~=1 && j~=length(a)
            dEda(j-1,k-1,l-1)=(E(j+1,k,l)-E(j-1,k,l))/(2*da);
            %end
            %if k~=1 && k~=length(b)
            dEdb(j-1,k-1,l-1)=(E(j,k+1,l)-E(j,k-1,l))/(2*db);
            %end
            %if l~=1 && l~=length(c)
            dEdc(j-1,k-1,l-1)=(E(j,k,l+1)-E(j,k,l-1))/(2*dc);
            %end
            %if j~=1 && j~=length(a) && k~=1 && k~=length(b)
            dEdadb(j-1,k-1,l-1)=(E(j+1,k+1)-E(j+1,k-1)-E(j-1,k+1)+E(j-1,k-1))...
            /(4*da*db);
            %end
            %if j~=1 && j~=length(a) && l~=1 && l~=length(c)
            dEdadc(j-1,k-1,l-1)=(E(j+1,k,l+1)-E(j+1,k,l-1)-E(j-1,k,l+1)+E(j-1,k,l-1))...
            /(4*da*dc);
            %end
            %if k~=1 && k~=length(b) && l~=1 && l~=length(c)
            dEdbdc(j-1,k-1,l-1)=(E(j,k+1,l+1)-E(j,k-1,l+1)-E(j,k+1,l-1)+E(j,k-1,l-1))...
            /(4*db*dc);
            %end
             %if j~=1 && j~=length(a) && k~=1 && k~=length(b) && l~=1 && l~=length(c)
            dEdadbdc(j-1,k-1,l-1) = (E(j+1,k+1,l+1)-E(j+1,k-1,l+1)-E(j-1,k+1,l+1)+E(j-1,k-1,l+1)...
                -E(j+1,k+1,l-1)+E(j+1,k-1,l-1)+E(j-1,k+1,l-1)-E(j-1,k-1,l-1))...
            /(8*da*db*dc);
             %end
        end
    end
end



% %Border points difference quotient 
% dEda(1,:,:) = (E(2,:,:)-E(1,:,:))/(a(2)-a(1));
% dEda(length(a),:,:) = (E(length(a),:,:)-E(length(a)-1,:,:))/(a(length(a))-a(length(a)-1));
% dEdb(:,1,:) = (E(:,2,:)-E(:,1,:))/(b(2)-b(1));
% dEdb(:,length(b),:)= (E(:,length(b),:)-E(:,length(b)-1,:))/(b(length(b))-b(length(b)-1));
% dEdc(:,:,1) = (E(:,:,2)-E(:,:,1))/(c(2)-c(1));
% dEdc(:,:,length(c))= (E(:,:,length(c))-E(:,:,length(c)-1))/(c(length(c))-c(length(c)-1));
% 
% 
% for i=2:length(b)-1
% dEdadb(1,i,:) = (E(2,i+1,:)-E(1,i+1,:)-E(2,i,:)+E(1,i,:))/((a(2)-a(1))*(b(i+1)-b(i)));
% dEdadb(length(a),i,:) = (E(length(a),i+1,:)-E(length(a)-1,i+1,:)-E(length(a),i,:)+E(length(a)-1,i,:))/(a(length(a))-a(length(a)-1)*(b(i+1)-b(i)));
% end
% for i=1:length(a)-1
% dEdadb(i,1,:) = (E(i+1,2,:)-E(i+1,1,:)-E(i,2,:)+E(i,1,:))/((a(i+1)-a(i))*(b(2)-b(1)));
% dEdadb(i,length(b),:) = (E(i+1,length(b),:)-E(i+1,length(b)-1,:)-E(i,length(b),:)+E(i,length(b)-1,:))/(b(length(b))-b(length(b)-1)*(a(i+1)-a(i)));
% end
% dEdadb(length(a),1,:) = (E(length(a),2,:)-E(length(a),1,:)-E(length(a)-1,2,:)+E(length(a)-1,1,:))/((a(length(a))-a(length(a)-1))*(b(2)-b(1)));
% dEdadb(length(a),length(b),:) = (E(length(a),length(b),:)-E(length(a),length(b)-1,:)-E(length(a)-1,length(b),:)+E(length(a)-1,length(b)-1,:))/((a(length(a))-a(length(a)-1))*(b(length(b))-b(length(b)-1)));
% 
% for i=2:length(c)-1
% dEdadc(1,:,i) = (E(2,:,i+1)-E(1,:,i+1)-E(2,:,i)+E(1,:,i))/((a(2)-a(1))*(c(i+1)-c(i)));
% dEdadc(length(a),:,i) = (E(length(a),:,i+1)-E(length(a)-1,:,i+1)-E(length(a),:,i)+E(length(a)-1,:,i))/(a(length(a))-a(length(a)-1)*(c(i+1)-c(i)));
% end
% for i=1:length(a)-1
% dEdadc(i,:,1) = (E(i+1,:,2)-E(i+1,:,1)-E(i,:,2)+E(i,:,1))/((a(i+1)-a(i))*(c(2)-c(1)));
% dEdadc(i,:,length(c)) = (E(i+1,:,length(c))-E(i+1,:,length(c)-1)-E(i,:,length(c))+E(i,:,length(c)-1))/(c(length(c))-c(length(c)-1)*(a(i+1)-a(i)));
% end
% dEdadc(length(a),:,1) = (E(length(a),:,2)-E(length(a),:,1)-E(length(a)-1,:,2)+E(length(a)-1,:,1))/((a(length(a))-a(length(a)-1))*(c(2)-c(1)));
% dEdadc(length(a),:,length(c)) = (E(length(a),:,length(c))-E(length(a),:,length(c)-1)-E(length(a)-1,:,length(c))+E(length(a)-1,:,length(c)-1))/((a(length(a))-a(length(a)-1))*(c(length(c))-c(length(c)-1)));
% 
% for i=2:length(c)-1
% dEdbdc(:,1,i) = (E(:,2,i+1)-E(:,1,i+1)-E(:,2,i)+E(:,1,i))/((b(2)-b(1))*(c(i+1)-c(i)));
% dEdbdc(:,length(b),i) = (E(:,length(b),i+1)-E(:,length(b)-1,i+1)-E(:,length(b),i)+E(:,length(b)-1,i))/(b(length(b))-b(length(b)-1)*(c(i+1)-c(i)));
% end
% for i=1:length(b)-1
% dEdbdc(i,:,1) = (E(:,i+1,2)-E(:,i+1,1)-E(:,i,2)+E(:,i,1))/((b(i+1)-b(i))*(c(2)-c(1)));
% dEdbdc(:,i,length(c)) = (E(:,i+1,length(c))-E(:,i+1,length(c)-1)-E(:,i,length(c))+E(:,i,length(c)-1))/(c(length(c))-c(length(c)-1)*(b(i+1)-b(i)));
% end
% dEdbdc(:,length(b),1) = (E(:,length(b),2)-E(:,length(b),1)-E(:,length(b)-1,2)+E(:,length(b)-1,1))/((b(length(b))-b(length(b)-1))*(c(2)-c(1)));
% dEdbdc(:,length(b),length(c)) = (E(:,length(b),length(c))-E(:,length(b),length(c)-1)-E(:,length(b)-1,length(c))+E(:,length(b)-1,length(c)-1))/((b(length(b))-b(length(b)-1))*(c(length(c))-c(length(c)-1)));
% 
% 
% for j=2:length(a)-1
%     for k=2:length(b)-1
%         dEdadbdc(j,k,1) = ((E(j+1,k+1,2)-E(j+1,k-1,2)-E(j-1,k+1,2)+E(j-1,k-1,2))-(E(j+1,k+1,1)-E(j+1,k-1,1)-E(j-1,k+1,1)+E(j-1,k-1,1)))...
%             /((a(j+1)-a(j-1))*(b(k+1)-b(k-1))*(c(2)-c(1)));
%         dEdadbdc(j,k,length(c)) = ((E(j+1,k+1,length(c))-E(j+1,k-1,length(c))-E(j-1,k+1,length(c))+E(j-1,k-1,length(c)))-(E(j+1,k+1,length(c)-1)-E(j+1,k-1,length(c)-1)-E(j-1,k+1,length(c)-1)+E(j-1,k-1,length(c)-1)))...
%             /((a(j+1)-a(j-1))*(b(k+1)-b(k-1))*(c(length(c))-c(length(c)-1)));
%     end
% end
% for j=2:length(a)-1
%     for l=2:length(c)-1
%         dEdadbdc(j,1,l) = ((E(j+1,2,l+1)-E(j+1,1,l+1)-E(j-1,2,l+1)+E(j-1,1,l+1))-(E(j+1,2,l-1)-E(j+1,1,l-1)-E(j-1,2,l-1)+E(j-1,1,l-1)))...
%             /((a(j+1)-a(j-1))*(b(2)-b(1))*(c(l+1)-c(l-1)));
%         dEdadbdc(j,length(b),l) = ((E(j+1,length(b),l+1)-E(j+1,length(b)-1,l+1)-E(j-1,length(b),l+1)+E(j-1,length(b)-1,l+1))-(E(j+1,length(b),l-1)-E(j+1,length(b)-1,l-1)-E(j-1,length(b),l-1)+E(j-1,length(b)-1,l-1)))...
%             /((a(j+1)-a(j-1))*(b(length(b))-b(length(b)-1))*(c(l+1)-c(l-1)));
%     end
% end
% 
% for k=2:length(b)-1
%     for l=2:length(c)-1
%         dEdadbdc(1,k,l) = ((E(2,k+1,l+1)-E(2,k-1,l+1)-E(1,k+1,l+1)+E(1,k-1,l+1))-(E(2,k+1,l-1)-E(2,k-1,l-1)-E(1,k+1,l-1)+E(1,k-1,l-1)))...
%             /((a(2)-a(1))*(b(k+1)-b(k-1))*(c(l+1)-c(l-1)));
%         dEdadbdc(length(a),k,l) = ((E(length(a),k+1,l+1)-E(length(a),k-1,l+1)-E(length(a)-1,k+1,l+1)+E(length(a)-1,k-1,l+1))-(E(length(a),k+1,l-1)-E(length(a),k-1,l-1)-E(length(a)-1,k+1,l-1)+E(length(a)-1,k-1,l-1)))...
%             /((a(length(a))-a(length(a)-1))*(b(k+1)-b(k-1))*(c(l+1)-c(l-1)));
%     end
% end
% 
% for j=2:length(a)-1
%    dEdadbdc(j,1,1) = (E(j+1,2,2)-E(j+1,1,2)-E(j-1,2,2)+E(j-1,1,2)...
%                 -E(j+1,2,1)+E(j+1,1,1)+E(j-1,2,1)-E(j-1,1,1))...
%             /((a(j+1)-a(j-1))*(b(2)-b(1))*(c(2)-c(1)));
%    dEdadbdc(j,length(b),1) = (E(j+1,length(b),2)-E(j+1,length(b)-1,2)-E(j-1,length(b),2)+E(j-1,length(b)-1,2)...
%                 -E(j+1,length(b),1)+E(j+1,length(b)-1,1)+E(j-1,length(b),1)-E(j-1,length(b)-1,1))...
%    /((a(j+1)-a(j-1))*(b(length(b))-b(length(b)-1))*(c(2)-c(1)));
%    dEdadbdc(j,1,length(c)) = (E(j+1,2,length(c))-E(j+1,1,length(c))-E(j-1,2,length(c))+E(j-1,1,length(c))...
%                 -E(j+1,2,length(c)-1)+E(j+1,1,length(c)-1)+E(j-1,2,length(c)-1)-E(j-1,1,length(c)-1))...
%             /((a(j+1)-a(j-1))*(b(2)-b(1))*(c(length(c))-c(length(c)-1)));
%    dEdadbdc(j,length(b),length(c)) = (E(j+1,length(b),length(c))-E(j+1,length(b)-1,length(c))-E(j-1,length(b),length(c))+E(j-1,length(b)-1,length(c))...
%            -E(j+1,length(b),length(c)-1)+E(j+1,length(b)-1,length(c)-1)+E(j-1,length(b),length(c)-1)-E(j-1,length(b)-1,length(c)-1))...
%            /((a(j+1)-a(j-1))*(b(length(b))-b(length(b)-1))*(c(length(c))-c(length(c)-1)));
% end
% 
% for k=2:length(b)-1
%    dEdadbdc(1,k,1) = (E(2,k+1,2)-E(2,k-1,2)-E(1,k+1,2)+E(1,k-1,2)...
%                 -E(2,k+1,1)+E(2,k-1,1)+E(1,k+1,1)-E(1,k-1,1))...
%             /((a(2)-a(1))*(b(k+1)-b(k-1))*(c(2)-c(1)));
%    dEdadbdc(length(a),k,1) = (E(length(a),k+1,2)-E(length(a),k-1,2)-E(length(a)-1,k+1,2)+E(length(a)-1,k-1,2)...
%                 -E(length(a),k+1,1)+E(length(a),k-1,1)+E(length(a)-1,k+1,1)-E(length(a)-1,k-1,1))...
%             /((a(length(a))-a(length(a)-1))*(b(k+1)-b(k-1))*(c(2)-c(1)));
%    dEdadbdc(1,k,length(c)) = (E(2,k+1,length(c))-E(2,k-1,length(c))-E(1,k+1,length(c))+E(1,k-1,length(c))...
%                 -E(2,k+1,length(c)-1)+E(2,k-1,length(c)-1)+E(1,k+1,length(c)-1)-E(1,k-1,length(c)-1))...
%             /((a(2)-a(1))*(b(k+1)-b(k-1))*(c(length(c))-c(length(c)-1)));
%    dEdadbdc(length(a),k,length(c)) =  (E(length(a),k+1,length(c))-E(length(a),k-1,length(c))-E(length(a)-1,k+1,length(c))+E(length(a)-1,k-1,length(c))...
%                 -E(length(a),k+1,length(c)-1)+E(length(a),k-1,length(c)-1)+E(length(a)-1,k+1,length(c)-1)-E(length(a)-1,k-1,length(c)-1))...
%             /((a(length(a))-a(length(a)-1))*(b(k+1)-b(k-1))*(c(length(c))-c(length(c)-1)));
% end
% 
% for l=2:length(c)-1
%    dEdadbdc(1,1,l) = (E(2,2,l+1)-E(2,1,l+1)-E(1,2,l+1)+E(1,1,l+1)...
%                 -E(2,2,l-1)+E(2,1,l-1)+E(1,1,l-1)-E(1,1,l-1))...
%             /((a(2)-a(1))*(b(2)-b(1))*(c(l+1)-c(l-1)));
%    dEdadbdc(1,length(b),l) = (E(2,length(b),l+1)-E(2,length(b)-1,l+1)-E(1,length(b),l+1)+E(1,length(b)-1,l+1)...
%                 -E(2,length(b),l-1)+E(2,length(b)-1,l-1)+E(1,length(b)-1,l-1)-E(1,length(b)-1,l-1))...
%             /((a(2)-a(1))*(b(length(b))-b(length(b)-1))*(c(l+1)-c(l-1)));
%    dEdadbdc(length(a),1,l) = (E(length(a),2,l+1)-E(length(a),1,l+1)-E(length(a)-1,2,l+1)+E(length(a)-1,1,l+1)...
%                 -E(length(a),2,l-1)+E(length(a),1,l-1)+E(length(a)-1,1,l-1)-E(length(a)-1,1,l-1))...
%             /((a(length(a))-a(length(a)-1))*(b(2)-b(1))*(c(l+1)-c(l-1)));
%    dEdadbdc(length(a),length(b),l) = (E(length(a),length(b),l+1)-E(length(a),length(b)-1,l+1)-E(length(a)-1,length(b),l+1)+E(length(a)-1,length(b)-1,l+1)...
%                 -E(length(a),length(b),l-1)+E(length(a),length(b)-1,l-1)+E(length(a)-1,length(b)-1,l-1)-E(length(a)-1,length(b)-1,l-1))...
%             /((a(length(a))-a(length(a)-1))*(b(length(b))-b(length(b)-1))*(c(l+1)-c(l-1)));
% end
%         
% dEdadbdc(1,1,1) = ((E(2,2,2)-E(2,1,2)-E(1,2,2)+E(1,1,2))-(E(2,2,1)-E(2,1,1)-E(1,2,1)+E(1,1,1)))...
%             /((a(2)-a(1))*(b(2)-b(1))*(c(2)-c(1)));
% dEdadbdc(length(a),1,1) = ((E(length(a),2,2)-E(length(a),1,2)-E(length(a)-1,2,2)+E(length(a)-1,1,2))-(E(length(a),2,1)-E(length(a),1,1)-E(length(a)-1,2,1)+E(length(a)-1,1,1)))...
%             /((a(length(a))-a(length(a)-1))*(b(2)-b(1))*(c(2)-c(1)));
% dEdadbdc(1,length(b),1) = ((E(2,length(b),2)-E(2,length(b)-1,2)-E(1,length(b),2)+E(1,length(b)-1,2))-(E(2,length(b),1)-E(2,length(b)-1,1)-E(1,length(b),1)+E(1,length(b)-1,1)))...
%             /((a(2)-a(1))*(b(length(b))-b(length(b)-1))*(c(2)-c(1)));
% dEdadbdc(1,1,length(c)) = ((E(2,2,length(c))-E(2,1,length(c))-E(1,2,length(c))+E(1,1,length(c)))-(E(2,2,length(c)-1)-E(2,1,length(c)-1)-E(1,2,length(c)-1)+E(1,1,length(c)-1)))...
%             /((a(2)-a(1))*(b(2)-b(1))*(c(length(c))-c(length(c)-1)));
% dEdadbdc(1,length(b),length(c)) = ((E(2,length(b),length(c))-E(2,length(b)-1,length(c))-E(1,length(b),length(c))+E(1,length(b)-1,length(c)))-(E(2,length(b),length(c)-1)-E(2,length(b)-1,length(c)-1)-E(1,length(b),length(c)-1)+E(1,length(b)-1,length(c)-1)))...
%             /((a(2)-a(1))*(b(length(b))-b(length(b)-1))*(c(length(c))-c(length(c)-1)));
% dEdadbdc(length(a),length(b),1) = ((E(length(a),length(b),2)-E(length(a),length(b)-1,2)-E(length(a)-1,length(b),2)+E(length(a)-1,length(b)-1,2))-(E(length(a),length(b),1)-E(length(a),length(b)-1,1)-E(length(a)-1,length(b),1)+E(length(a)-1,length(b)-1,1)))...
%             /((a(length(a))-a(length(a)-1))*(b(length(b))-b(length(b)-1))*(c(2)-c(1)));
% dEdadbdc(length(a),1,length(c)) = ((E(length(a),2,length(c))-E(length(a),1,length(c))-E(length(a)-1,2,length(c))+E(length(a)-1,1,length(c)))-(E(length(a),2,length(c)-1)-E(length(a),1,length(c)-1)-E(length(a)-1,2,length(c)-1)+E(length(a)-1,1,length(c)-1)))...
%             /((a(length(a))-a(length(a)-1))*(b(2)-b(1))*(c(length(c))-c(length(c)-1)));
% dEdadbdc(length(a),length(b),length(c)) = ((E(length(a),length(b),length(c))...
% -E(length(a),length(b)-1,length(c))-E(length(a)-1,length(b),length(c))...
% +E(length(a)-1,length(b)-1,length(c)))-(E(length(a),length(b),length(c)-1)...
% -E(length(a),length(b)-1,length(c)-1)-E(length(a)-1,length(b),length(c)-1)...
% +E(length(a)-1,length(b)-1,length(c)-1)))...
%             /((a(length(a))-a(length(a)-1))*(b(length(b))-b(length(b)-1))*(c(length(c))-c(length(c)-1)));
 end