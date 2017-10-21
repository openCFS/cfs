function [dEda, dEdb, dEdc, dEdadb, dEdbdc, dEdadc,dEdadbdc] = tricubic_partialderiv(a,b,c,E)
dEda = zeros(length(a),length(b),length(c));
dEdb = zeros(length(a),length(b),length(c));
dEdc = zeros(length(a),length(b),length(c));
dEdadb = zeros(length(a),length(b),length(c));
dEdadc = zeros(length(a),length(b),length(c));
dEdbdc = zeros(length(a),length(b),length(c));
dEdadbdc = zeros(length(a),length(b),length(c));

%Central Difference quotient for inner points
for j=2:length(a)-1
    for k=2:length(b)-1
        for l = 2:length(c)-1
            dEda(j,k,l)=(E(j+1,k,l)-E(j-1,k,l))/2.;
            dEdb(j,k,l)=(E(j,k+1,l)-E(j,k-1,l))/2.;
            dEdc(j,k,l)=(E(j,k,l+1)-E(j,k,l-1))/2.;
            dEdadb(j,k,l)=(E(j+1,k+1)-E(j+1,k-1)-E(j-1,k+1)+E(j-1,k-1))...
            /4.;
            dEdadc(j,k,l)=(E(j+1,k,l+1)-E(j+1,k,l-1)-E(j-1,k,l+1)+E(j-1,k,l-1))...
            /4.;
            dEdbdc(j,k,l)=(E(j,k+1,l+1)-E(j,k-1,l+1)-E(j,k+1,l-1)+E(j,k-1,l-1))...
            /4.;
            dEdadbdc(j,k,l) = (E(j+1,k+1,l+1)-E(j+1,k-1,l+1)-E(j-1,k+1,l+1)+E(j-1,k-1,l+1)...
                -E(j+1,k+1,l-1)+E(j+1,k-1,l-1)+E(j-1,k+1,l-1)-E(j-1,k-1,l-1))...
            /8.;
        end
    end
end





%Border points difference quotient 
dEda(1,:,:) = (E(2,:,:)-E(1,:,:))/(a(2)-a(1));
dEda(length(a),:,:) = (E(length(a),:,:)-E(length(a)-1,:,:))/2.;
dEdb(:,1,:) = (E(:,2,:)-E(:,1,:))/(b(2)-b(1));
dEdb(:,length(b),:)= (E(:,length(b),:)-E(:,length(b)-1,:))/2.;
dEdc(:,:,1) = (E(:,:,2)-E(:,:,1))/(c(2)-c(1));
dEdc(:,:,length(c))= (E(:,:,length(c))-E(:,:,length(c)-1))/2.;


for i=2:length(b)-1
dEdadb(1,i,:) = (E(2,i+1,:)-E(1,i+1,:)-E(2,i,:)+E(1,i,:))/4.;
dEdadb(length(a),i,:) = (E(length(a),i+1,:)-E(length(a)-1,i+1,:)-E(length(a),i,:)+E(length(a)-1,i,:))/4.;
end
for i=1:length(a)-1
dEdadb(i,1,:) = (E(i+1,2,:)-E(i+1,1,:)-E(i,2,:)+E(i,1,:))/4.;
dEdadb(i,length(b),:) = (E(i+1,length(b),:)-E(i+1,length(b)-1,:)-E(i,length(b),:)+E(i,length(b)-1,:))/4.;
end
dEdadb(length(a),1,:) = (E(length(a),2,:)-E(length(a),1,:)-E(length(a)-1,2,:)+E(length(a)-1,1,:))/4.;
dEdadb(length(a),length(b),:) = (E(length(a),length(b),:)-E(length(a),length(b)-1,:)-E(length(a)-1,length(b),:)+E(length(a)-1,length(b)-1,:))/4.;

for i=2:length(c)-1
dEdadc(1,:,i) = (E(2,:,i+1)-E(1,:,i+1)-E(2,:,i)+E(1,:,i))/4.;
dEdadc(length(a),:,i) = (E(length(a),:,i+1)-E(length(a)-1,:,i+1)-E(length(a),:,i)+E(length(a)-1,:,i))/4.;
end
for i=1:length(a)-1
dEdadc(i,:,1) = (E(i+1,:,2)-E(i+1,:,1)-E(i,:,2)+E(i,:,1))/4.;
dEdadc(i,:,length(c)) = (E(i+1,:,length(c))-E(i+1,:,length(c)-1)-E(i,:,length(c))+E(i,:,length(c)-1))/4.;
end
dEdadc(length(a),:,1) = (E(length(a),:,2)-E(length(a),:,1)-E(length(a)-1,:,2)+E(length(a)-1,:,1))/4.;
dEdadc(length(a),:,length(c)) = (E(length(a),:,length(c))-E(length(a),:,length(c)-1)-E(length(a)-1,:,length(c))+E(length(a)-1,:,length(c)-1))/4.;

for i=2:length(c)-1
dEdbdc(:,1,i) = (E(:,2,i+1)-E(:,1,i+1)-E(:,2,i)+E(:,1,i))/((b(2)-b(1))*(c(i+1)-c(i)));
dEdbdc(:,length(b),i) = (E(:,length(b),i+1)-E(:,length(b)-1,i+1)-E(:,length(b),i)+E(:,length(b)-1,i))/4.;
end
for i=1:length(b)-1
dEdbdc(i,:,1) = (E(:,i+1,2)-E(:,i+1,1)-E(:,i,2)+E(:,i,1))/((b(i+1)-b(i))*(c(2)-c(1)));
dEdbdc(:,i,length(c)) = (E(:,i+1,length(c))-E(:,i+1,length(c)-1)-E(:,i,length(c))+E(:,i,length(c)-1))/4.;
end
dEdbdc(:,length(b),1) = (E(:,length(b),2)-E(:,length(b),1)-E(:,length(b)-1,2)+E(:,length(b)-1,1))/4.;
dEdbdc(:,length(b),length(c)) = (E(:,length(b),length(c))-E(:,length(b),length(c)-1)-E(:,length(b)-1,length(c))+E(:,length(b)-1,length(c)-1))/4.;


for j=2:length(a)-1
    for k=2:length(b)-1
        dEdadbdc(j,k,1) = ((E(j+1,k+1,2)-E(j+1,k-1,2)-E(j-1,k+1,2)+E(j-1,k-1,2))-(E(j+1,k+1,1)-E(j+1,k-1,1)-E(j-1,k+1,1)+E(j-1,k-1,1)))...
            /8.;
        dEdadbdc(j,k,length(c)) = ((E(j+1,k+1,length(c))-E(j+1,k-1,length(c))-E(j-1,k+1,length(c))+E(j-1,k-1,length(c)))-(E(j+1,k+1,length(c)-1)-E(j+1,k-1,length(c)-1)-E(j-1,k+1,length(c)-1)+E(j-1,k-1,length(c)-1)))...
            /8.;
    end
end
for j=2:length(a)-1
    for l=2:length(c)-1
        dEdadbdc(j,1,l) = ((E(j+1,2,l+1)-E(j+1,1,l+1)-E(j-1,2,l+1)+E(j-1,1,l+1))-(E(j+1,2,l-1)-E(j+1,1,l-1)-E(j-1,2,l-1)+E(j-1,1,l-1)))...
            /8.;
        dEdadbdc(j,length(b),l) = ((E(j+1,length(b),l+1)-E(j+1,length(b)-1,l+1)-E(j-1,length(b),l+1)+E(j-1,length(b)-1,l+1))-(E(j+1,length(b),l-1)-E(j+1,length(b)-1,l-1)-E(j-1,length(b),l-1)+E(j-1,length(b)-1,l-1)))...
            /8.;
    end
end

for k=2:length(b)-1
    for l=2:length(c)-1
        dEdadbdc(1,k,l) = ((E(2,k+1,l+1)-E(2,k-1,l+1)-E(1,k+1,l+1)+E(1,k-1,l+1))-(E(2,k+1,l-1)-E(2,k-1,l-1)-E(1,k+1,l-1)+E(1,k-1,l-1)))...
            /8.;
        dEdadbdc(length(a),k,l) = ((E(length(a),k+1,l+1)-E(length(a),k-1,l+1)-E(length(a)-1,k+1,l+1)+E(length(a)-1,k-1,l+1))-(E(length(a),k+1,l-1)-E(length(a),k-1,l-1)-E(length(a)-1,k+1,l-1)+E(length(a)-1,k-1,l-1)))...
            /8.;
    end
end

dEdadbdc(1,1,1) = ((E(2,2,2)-E(2,1,2)-E(1,2,2)+E(1,1,2))-(E(2,2,1)-E(2,1,1)-E(1,2,1)+E(1,1,1)))...
            /8.;
dEdadbdc(length(a),1,1) = ((E(length(a),2,2)-E(length(a),1,2)-E(length(a)-1,2,2)+E(length(a)-1,1,2))-(E(length(a),2,1)-E(length(a),1,1)-E(length(a)-1,2,1)+E(length(a)-1,1,1)))...
            /8.;
dEdadbdc(1,length(b),1) = ((E(2,length(b),2)-E(2,length(b)-1,2)-E(1,length(b),2)+E(1,length(b)-1,2))-(E(2,length(b),1)-E(2,length(b)-1,1)-E(1,length(b),1)+E(1,length(b)-1,1)))...
            /8.;
dEdadbdc(1,1,length(c)) = ((E(2,2,length(c))-E(2,1,length(c))-E(1,2,length(c))+E(1,1,length(c)))-(E(2,2,length(c)-1)-E(2,1,length(c)-1)-E(1,2,length(c)-1)+E(1,1,length(c)-1)))...
            /8.;
dEdadbdc(1,length(b),length(c)) = ((E(2,length(b),length(c))-E(2,length(b)-1,length(c))-E(1,length(b),length(c))+E(1,length(b)-1,length(c)))-(E(2,length(b),length(c)-1)-E(2,length(b)-1,length(c)-1)-E(1,length(b),length(c)-1)+E(1,length(b)-1,length(c)-1)))...
            /8.;
dEdadbdc(length(a),length(b),1) = ((E(length(a),length(b),2)-E(length(a),length(b)-1,2)-E(length(a)-1,length(b),2)+E(length(a)-1,length(b)-1,2))-(E(length(a),length(b),1)-E(length(a),length(b)-1,1)-E(length(a)-1,length(b),1)+E(length(a)-1,length(b)-1,1)))...
            /8.;
dEdadbdc(length(a),1,length(c)) = ((E(length(a),2,length(c))-E(length(a),1,length(c))-E(length(a)-1,2,length(c))+E(length(a)-1,1,length(c)))-(E(length(a),2,length(c)-1)-E(length(a),1,length(c)-1)-E(length(a)-1,2,length(c)-1)+E(length(a)-1,1,length(c)-1)))...
            /8.;
dEdadbdc(length(a),length(b),length(c)) = ((E(length(a),length(b),length(c))...
-E(length(a),length(b)-1,length(c))-E(length(a)-1,length(b),length(c))...
+E(length(a)-1,length(b)-1,length(c)))-(E(length(a),length(b),length(c)-1)...
-E(length(a),length(b)-1,length(c)-1)-E(length(a)-1,length(b),length(c)-1)...
+E(length(a)-1,length(b)-1,length(c)-1)))...
            /8.;
 end