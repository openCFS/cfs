function [x] = createCross(nelx,nely,a,b)

if a > 1
    warning('MATLAB:createCross','parameter a too big');
    a = 1;
end
if b > 1
    warning('MATLAB:createCross','parameter b too big');
    b = 1;
end

x = 1e-6*ones(nely,nelx);

nela = round(a*nelx);
nelb = round(b*nely);

if nelx-nela == 1
    nela = nelx;
end
if nely-nelb == 1
    nelb = nely;
end

% Vertical bar
ela = ceil(nelx/2-(nela-1)/2):ceil(nelx/2+(nela-1)/2);
x(:,ela) = 1;

% Horizontal bar
elb = ceil(nely/2-(nelb-1)/2):ceil(nely/2+(nelb-1)/2);
x(elb,:) = 1;
