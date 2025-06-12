%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%
function [dcn]=filterE(nelx,nely,rmin,x,dc)
dcn=zeros(3,3,nelx,nely);
for i = 1:nelx
    for j = 1:nely
        sum=0.0;
        for k = max(i-round(rmin),1):min(i+round(rmin),nelx)
            for l = max(j-round(rmin),1):min(j+round(rmin), nely)
                fac = rmin-sqrt((i-k)^2+(j-l)^2);
                sum = sum+max(0,fac);
                dcn(:,:,i,j) = dcn(:,:,i,j) + max(0,fac)*x(k,l)*dc(:,:,k,l);
            end
        end
        dcn(:,:,i,j) = dcn(:,:,i,j)/(x(i,j)*sum);
    end
end