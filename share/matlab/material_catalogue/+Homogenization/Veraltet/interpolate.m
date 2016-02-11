function [val, newdata] = interpolate(point, data)
% interpolate interpolates the value at point given the matrix data
% in the [0,1]-hypercube, where
% data = [p1x p1y ... p1val1 p1val2 ...
%         p2x p2y ... p2val1 p2val2 ...
%         p3x p3y ... p3val1 p3val2 ...]

n = length(point);

% If the point is given in data take the interpolation values from there
datapoints = data(:,1:n);
[~,idx] = ismember(point,datapoints,'rows');
if idx > 0
    val = data(idx, n+1:end);
    newdata = data;
    return
end

nonboundaryentries = n - sum( point == 0 | point == 1);

intval = zeros(nonboundaryentries * 2, size(data,2) - n + 1);

% Project point to the boundary of the current hyperplane
% -> two new points in each dimension on parallel sub-hyperplanes
c = 1;
for i = 1:n
    if point(i) > 0 && point(i) < 1
        newpoint = point;
        newpoint(i) = 0;
        [tmp, data] = interpolate(newpoint, data);
        intval(c,:) = [tmp, norm(newpoint - point)];
        c = c + 1;
        newpoint(i) = 1;
        [tmp, data] = interpolate(newpoint, data);
        intval(c,:) = [tmp, norm(newpoint - point)];
        c = c + 1;
    end
end

% Perform weighted linear interpolation
for i = 1:size(intval,1)
    intval(i,1:end-1) = intval(i,1:end-1) .* intval(i,end);
end
sumval = sum(intval);
val = sumval(1:end-1) ./ sumval(end);

% Add the new point to the data set
newdata = [data;point,val];
