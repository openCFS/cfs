function [pp] = computePP(fhandle,nelx,nely,lx,ly)

str = strcat('bmat',num2str(nelx));
load(str);

nnodx = nelx+1;
nnody = nely+1;

[E0] = createISO(10,.3);
[u0] = createu0(nnodx, nnody, lx, ly);
Ttrans=td(nnodx,nnody);
createA2;

% Compute some homogenized elasticity tensors
presets = 1e-6:.1:1;
npresets = length(presets);
loop = 1;
for i = presets
    a = i;
    b = i;
    x = fhandle(nelx,nely,lx,ly,a,b);
%     imagesc(reshape(x,nelx,nely));axis tight;axis equal;drawnow;
    Eh{loop} = computeEh(reshape(x,nelx*nely,1),nnodx*nnody*2,BB,pVol,Ttrans,u0);
    loop = loop+1;
end

% Interpolation
for i=1:3
    for j=1:3
        val = zeros(1,npresets);
        for k=1:npresets
            val(k) = Eh{k}(i,j);
        end
        pp(i,j) = pchip(presets,val);
    end
end

% E0
% 
% disp ('Poisson');
% E0(2,1)/E0(1,1)
% 
% disp ('Y.M.');
% E0(1,1)*(1-E0(2,1)^2/E0(1,1)^2)
