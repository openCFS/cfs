function [x] = createCrossShear(nelx,nely,a0,b0,phi)
% Creates a periodic micro structure consisting of sheared crosses.
% @param:
% nelx     number of elements in x direction
% nely     number of elements in y direction
% a0       width  of horizontal bars, 0<=a0<=1
% b0       width  of vertical bars, 0<=b0<=1
% phi      shearing angle, -90<=phi<=90
% 

%% Convert phi to radian
phi = phi/180*pi;

% if nx>5e3
%     warning('createCrossShear:MeshTooBig','Meshsize is over limit: nx = %d. Exiting',nx);
%     x = -1;
%     return;
% end
x = 1e-6*ones(nelx,nely);

% get barwidth
nela = round(a0/2*nelx);
nelb = round(b0/2*nely);

%% create horizontal bar
if mod(nelx,2)
    ela = ceil(nelx/2)-nela+1:ceil(nelx/2)+nela;
else
    ela = nelx/2-nela+1:nelx/2+nela;
end
for elx = ela
    for ely = 1:nely
        x(elx,ely) = 1;
    end
end

%% create rotated vertical bar
if abs(phi) <= pi/2-1e-6
    width = nelb/cos(phi); % half of horizontal barwidth
    for elx = 1:nelx
        mid = ceil(nely/2)+(elx-nelx/2-1)*tan(phi);
        % get min and max elements of main bar
        elbMin = mid-width+1;
        elbMax = mid+width;
        el = round(elbMin):round(elbMax);
        % select only feasible elements
        el = el(el>0);
        el = el(el<=nely);
        for ely = el
            x(elx,ely) = 1;
        end
    end
end
end