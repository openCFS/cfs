% Move values without tensor from end of cataloguefile to new file then call this script

cataloguefile = 'detailed_stats_presets3D_10';
file1 = 'detailed_stats_1';
presetsfile = 'presets3D_10';
level = 7;

data = dlmread(cataloguefile,'\t',1,0);

% Copy original file
copyfile(cataloguefile,[cataloguefile, '_copy'],'f');

points = load(file1);

% Write missing tensors at the end of the cataloguefile copy
% The copy has the same format as returned from
% Homogenization.getInterpolationValuesFromFile
format = [repmat('%.10f,',1,size(points,2)-1),'%.10f\n'];
fid = fopen([cataloguefile, '_copy'],'at');
fprintf(fid,format,points');
fclose(fid);


% Scale points
points = (points + 1) / 2;

fid = fopen(cataloguefile,'at');

% Due to the symmetry we first check if the cell belonging to the current 
% point but rotated by pi/2 exists. In this case we rotate the material
% tensor.
n = size(points,2);
Tensor = zeros(3,3);
intPoints = [];
for i = 1:size(points,1)
    % check if rotated cell exists
%     [ism, idx] = ismemberf(points(i,[2 1 4 3]), data(:,1:n), 'rows');
    ism = 0;
    if ism
        % get tensor
        vec = data(idx,n+1:end);
        Tensor(1,2:3) = vec(2:3);
        Tensor(2,3) = vec(5);
        Tensor = Tensor + Tensor';
        Tensor(1,1) = vec(1);
        Tensor(2,2) = vec(4);
        Tensor(3,3) = vec(6);
        % rotate tensor
        Rot = hillmandeltovoigt( rotatehillmandel( voigttohillmandel(Tensor), pi/2) );
        % append volume
        val = [Rot(find(tril(Rot))); vec(end)];
        % append values to data
        data = [data; points(i,:), val'];
        % write values to catalogue file
        fprintf(fid, '%.10f\t', points(i,:));
        fprintf(fid, '%e\t', val(1:end-1));
        fprintf(fid, '%e\n', val(end));
    else
        intPoints = [intPoints; points(i,:)];
    end
end

% If no rotated cell exists we linearly interpolate the missing tensor.
for i = 1:size(intPoints,1)
    val = interpolate(intPoints(i,:), data);
    if isnan(val)
        disp( strcat('No interpolation for ', sprintf('%f\t',intPoints(i,:)), sprintf('\n')) );
    end
    if size(val,1) ~= 1
        val = val';
    end
    data = [data; points(i,:), val];
    fprintf(fid, '%.10f\t', intPoints(i,:));
    fprintf(fid, '%e\t', val(1:end-1));
    fprintf(fid, '%e\n', val(end));
end

fclose(fid);

sortCatalogue(data,presetsfile);