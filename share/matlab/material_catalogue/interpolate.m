function [val, newdata] = interpolate(point, data)
% INTERPOLATE interpolates the value at point given the matrix data
% in the [0,1]-hypercube, where
% data = [p1x p1y ... p1val1 p1val2 ...
%         p2x p2y ... p2val1 p2val2 ...
%         p3x p3y ... p3val1 p3val2 ...]

n = length(point);

% If the point is given in data take the interpolation values from there
datapoints = data(:,1:n);
[~,idx] = ismemberf(point,datapoints,'rows');
if idx > 0
    val = data(idx, n+1:end);
    newdata = data;
    return
end

values = zeros(n,size(data,2));
dist = zeros(n,1);

for dim = 1:n
    [~,idcs] = ismemberf(data(:,setdiff(1:n,dim)),point(setdiff(1:n,dim)),'rows');
    if sum(idcs) < 1.5
        dist(dim) = Inf;
        values(dim,:) = NaN * zeros(1,size(data,2));
        continue
    end
    B = data(idcs|zeros(size(data,1),1),:);
    [vec,I] = sort(B(:,dim));
    B = B(I,:);

    idx1 = find(vec < point(dim),1,'last');
    idx2 = find(vec > point(dim),1);

    if isempty(idx1)
        tmp = find(vec > point(dim),2);
        idx1 = tmp(2);
    end

    if isempty(idx2)
        tmp = find(vec < point(dim),2);
        idx2 = tmp(2);
    end
    
    dist(dim) = vec(idx1) - vec(idx2);

    m = (B(idx1,:) - B(idx2,:)) / dist(dim);

    values(dim,:) = B(idx1,:) + m * (point(dim) - vec(idx1));
    
end
mini = min(abs(dist));

intValue = mean(values((abs(dist) - mini) < eps,:),1);


newdata = [data; intValue];
val = intValue(n+1:end);
