% function [dEda, dEdb, dEdadb] = bicubic_partialderiv(a,b,E,E_number)
% % Approximation der Ableitungen in den Interpolationspunkten mit Finiten Differenzen
% dEda = zeros(length(a),length(b));
% dEdb = zeros(length(a),length(b));
% dEdadb = zeros(length(a),length(b));
% %Central Difference quotient for inner points
% for j=2:length(a)-1
%     for k=2:length(b)-1
%         dEda(j,k)=(E(j+1,k)-E(j-1,k))/(a(j+1)-a(j-1));
%         dEdb(j,k)=(E(j,k+1)-E(j,k-1))/(b(k+1)-b(k-1));
%         dEdadb(j,k)=(E(j+1,k+1)-E(j+1,k-1)-E(j-1,k+1)+E(j-1,k-1))...
%         /((a(j+1)-a(j-1))*(b(k+1)-b(k-1)));
%     end
% end
% if E_number == 1
%     dEda(:,1) = 0.;
%     dEda(:,length(b)) = 0.;
%     dEda(1,:) =  (E(2,:)-E(1,:))/(a(2)-a(1));
%     dEda(length(a),:) = (E(length(a),:)-E(length(a)-1,:))/(a(length(a))-a(length(a)-1));
%     
%     dEdb(1,:) = (E(1,2)-E(1,1))/(b(2)-b(1));
%     dEdb(length(a),:) = 0.;
%     dEdb(:,length(b)) = (E(:,length(b))-E(:,length(b)-1))/(b(length(b))-b(length(b)-1));
%     dEdb(:,1) = (E(:,2)-E(:,1))/(b(2)-b(1));
% elseif E_number ==2
%     dEda(:,1) = (E(2,1)-E(1,1))/(a(2)-a(1));
%     dEda(:,length(b)) = 0.;
%     dEda(1,:) =  (E(2,:)-E(1,:))/(a(2)-a(1));
%     dEda(length(a),:) = 2*(E(length(a),:)-E(length(a)-1,:))/(a(length(a))-a(length(a)-1));
%     
%     dEdb(1,:) = 0;
%     dEdb(length(a),:) = 0.;
%     dEdb(:,length(b)) = (E(:,length(b))-E(:,length(b)-1))/(b(length(b))-b(length(b)-1));
%     dEdb(:,1) = (E(:,2)-E(:,1))/(b(2)-b(1));
% end
% end

function [dEda, dEdb, dEdadb] = bicubic_partialderiv(a,b,E,deriv_a,deriv_b,deriv_a2,deriv_b2)
% Approximation der Ableitungen in den Interpolationspunkten mit Finiten Differenzen
% vectors deriv_a and deriv_b are fixed derivatives used by the penalization
% approach
dEda = zeros(length(a),length(b));
dEdb = zeros(length(a),length(b));
dEdadb = zeros(length(a),length(b));
% Central difference quotient for inner points
for j=2:length(a)-1
    for k=2:length(b)-1
        dEdadb(j,k)=(E(j+1,k+1)-E(j+1,k-1)-E(j-1,k+1)+E(j-1,k-1))...
        /((a(j+1)-a(j-1))*(b(k+1)-b(k-1)));
    end
end
% Central difference quotient for inner points
for j=2:length(a)-2
    dEda(j,:)=(E(j+1,:)-E(j-1,:))/(a(j+1)-a(j-1));
end
for j=2:length(b)-2
    dEdb(:,j)=(E(:,j+1)-E(:,j-1))/(b(j+1)-b(j-1));
end

%Backwards difference quotient for outer points and next to outer points.
%As a=0 & b=1 and a=1 & b=0 are treated as full material, these points can
%not be used for central difference quotients.
dEda(1,:) = (E(2,:)-E(1,:))/(a(2)-a(1));
dEda(length(a)-1,:) = (E(length(a)-1,:)-E(length(a)-2,:))/(a(length(a)-1)-a(length(a)-2));
dEdb(:,1) = (E(:,2)-E(:,1))/(b(2)-b(1));
dEdb(:,length(b)-1) = (E(:,length(b)-1)-E(:,length(b)-2))/(b(length(b)-1)-b(length(b)-2));
if isempty(deriv_a)
    dEda(length(a),:) = (E(length(a),:)-E(length(a)-1,:))/(a(length(a))-a(length(a)-1));
else
    dEda(length(a),:) = deriv_a(:);
    dEda(:,length(b)) = deriv_a2(:);
end
if isempty(deriv_b)
    dEdb(:,length(b)) = (E(:,length(b))-E(:,length(b)-1))/(b(length(b))-b(length(b)-1));
else
    dEdb(length(a),:) = deriv_b(:);
    dEdb(:,length(b)) = deriv_b2(:);
end
