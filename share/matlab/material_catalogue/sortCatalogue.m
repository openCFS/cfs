function sortCatalogue(arg1, file2)
% @args:
%   arg1    data as matrix or cataloguefile
%   file2   presets

if ischar(arg1)
    try
        cataloguefile = arg1;
        firstentry = dlmread(arg1,'\t',[0 0 0 0]);
        if firstentry > 1
            data = dlmread(arg1,'\t',1,0);
            dim = firstentry;
        else
            data = load(arg1);
            dim = size(data,2) - 7;
        end
    catch
        disp('Could not read file.');
        exit;
    end
else
    tmp = strsplit(file2,{'/','\'});
    cataloguefile = ['catalogues/detailed_stats_',tmp{end}];
    firstentry = arg1(1);
    data = arg1;
    dim = size(data,2) - 7;
end

tmp = strsplit(file2,'_');
tmp = tmp{end};
level = str2double(tmp(2:end));

points = load(file2);
points(1,:) = [];
if min(points(:,1)) < 0
    points = (points + 1)/2;
end
dataPoints = data(:,1:dim);
points = points(:,1:dim);

[~,indeces] = ismember(points,dataPoints,'rows');
idcs = find(indeces == 0);

for i = idcs
    [~,idx] = ismemberf(points(i,:),dataPoints,'rows');
    indeces(i) = idx;
end

% Sort data to fit mapping in Sparse Grid
data1 = data(indeces,:);

% Save copy of catalogue
copyfile(cataloguefile,[cataloguefile, '_copy_before_sort'],'f');

fid = fopen(cataloguefile,'wt');
if firstentry == .5
    format = [ '%dD\tL%d\t%d\t%s', repmat('\t%e',1,7-4+dim) , '\n'];
else
    format = [ '%dD\t%d\t%d\t%s', repmat('\t%e',1,7-4+dim) , '\n'];
end
fprintf(fid,format,size(points,2),level,size(points,1),'voigt',zeros(1,7-4+dim));
format = [ repmat('%.10f\t',1,size(points,2)) , repmat('%e\t',1,6), '%e\n'];
fprintf(fid,format,data1');
fclose(fid);

