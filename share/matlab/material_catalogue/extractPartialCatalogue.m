cataloguefile = 'catalogues/detailed_stats_presets3D_63';
presetsfile = 'presets/byCoords/presets3D_L5';

fid = fopen(cataloguefile);
tmp = fscanf(fid,'%dD\t%d');
dim = tmp(1);
fclose(fid);
if ~exist('fulldata','var')
    fulldata = dlmread(cataloguefile,'\t',1,0);
end

points = load(presetsfile);
points(1,:) = [];
if min(points(:,1)) < 0
    points = (points+1)/2;
end
points = points(:,1:dim);

idx = ismember(fulldata(:,1:dim),points,'rows');
partialdata = fulldata(idx,:);

assert( size(partialdata,1)==size(points,1) )

sortCatalogue(partialdata,presetsfile);

clear fid tmp nPoints idx ans
