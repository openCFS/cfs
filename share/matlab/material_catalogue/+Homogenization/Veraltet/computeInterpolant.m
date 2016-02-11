function [Interpolant,Eh11,Eh12,Eh22,Eh33] = computeInterpolant(fhandle,nelx,nely,lx,ly)

str = strcat('bmat',num2str(nelx));
load(str);

nnodx = nelx+1;
nnody = nely+1;

[E0] = createISO(1,.3);
[u0] = createu0(nnodx, nnody, lx, ly);
Ttrans=td(nnodx,nnody);
createA2;

% Compute some homogenized elasticity tensors
presetsa = 0:.1:1;
presetsb = presetsa; % -> Eh symmetrisch
npresetsa = numel(presetsa);
npresetsb = numel(presetsb);
Eh11 = zeros(npresetsa,npresetsb);
Eh12 = zeros(npresetsa,npresetsb);
Eh22 = zeros(npresetsa,npresetsb);
Eh33 = zeros(npresetsa,npresetsb);
loopa = 1;
for a = presetsa
    loopb = 1;
    for b = presetsb
        x = fhandle(nelx,nely,a,b);
%         imagesc(reshape(x,nelx,nely));axis tight;axis
%         equal;drawnow;pause(.5);
        Eh = computeEh(reshape(x,nelx*nely,1),nnodx*nnody*2,BB,pVol,Ttrans,u0);
        Eh11(loopa,loopb) = Eh(1,1);
        Eh12(loopa,loopb) = Eh(1,2);
        Eh22(loopa,loopb) = Eh(2,2);
        Eh33(loopa,loopb) = Eh(3,3)/2;
        loopb = loopb+1;
    end
    loopa = loopa+1;
end

% Save Eh
filename = strcat('detailed_stats_',num2str(npresetsa-1));
fid = fopen(filename,'w');
fprintf(fid,'%d\t%d\t%e\t%e\t%e\t%e\r\n',npresetsa-1,npresetsa-1,0,0,0,0);
for i=1:npresetsa
    for j=1:i
        fprintf(fid,'%d\t%d\t%e\t%e\t%e\t%e\r\n',i-1,j-1,Eh11(i,j),Eh12(i,j),Eh22(i,j),Eh33(i,j));
        if i~=j
            fprintf(fid,'%d\t%d\t%e\t%e\t%e\t%e\r\n',j-1,i-1,Eh11(i,j),Eh12(i,j),Eh22(i,j),Eh33(i,j));
        end
    end
end
fclose(fid);

% Get Interpolants
Interpolant = cell(3);
for i=1:3
    for j=1:3
        if i==1 && j==1
            Val = Eh11;
        elseif (i==1 && j==2) || (i==2 && j==1)
            Val = Eh12;
        elseif i==2 && j==2
            Val = Eh22;
        elseif i==3 && j==3
            Val = Eh33;
        else
            Val = zeros(npresetsa,npresetsb);
        end
        [X,Y] = ndgrid(presetsa,presetsb);
        E = griddedInterpolant(X,Y,Val,'cubic');
        [Eb,Ea] = gradient(Val,presetsb,presetsa);
        dEa = griddedInterpolant(X,Y,Ea,'cubic');
        dEb = griddedInterpolant(X,Y,Eb,'cubic');
        Interpolant{i,j} = {E;dEa;dEb};
    end
end
