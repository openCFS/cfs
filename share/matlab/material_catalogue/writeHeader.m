function writeHeader(file)

% Load data from file
data = load(file);

if any(data(:,1) > 1) % -> point given by level and index
    givenByLevelAndIndex = 1;
    dim = ( find( abs( data(2,:) - round(data(2,:)) ), 1)-1 ) / 2;
    level = max(data(:,1));
    m = size(data,2)-2*dim;
else
    givenByLevelAndIndex = 0;
    dim = size(data,2)-7;
    level = round(size(data,1)^(1/dim),12);
    if ~(mod(level,1) == 0)
        level = log( numel(unique(data(:,1)))+1 ) / log(2);
    end
    m = size(data,2)-dim;
end

% Write header and data
fid = fopen(file,'wt');
if givenByLevelAndIndex
    % sparsegrid N d m voigt 0 0 0...
    format = 'sparsegrid_ver2\t%d\t%d\t%d\t%s\t%s\n';
    fprintf(fid,format,size(data,1),dim,m,'not_hierarchized','voigt');
    % L1 L2 L3 I1 I2 I3 E11 E12 ...
    format = [ repmat('%d\t',1,2*dim), repmat('%e\t',1,m-1), '%e\n'];
    fprintf(fid,format,data');
else
    if data(1,1) == .5
        % 3D L5 N voigt 0 0 0...
        format = [ '%dD\tL%d\t%d\t%s', repmat('\t%e',1,m-4+dim) , '\n'];
    else
        % 3D 1 N voigt 0 0 0...
        format = [ '%dD\t%d\t%d\t%s', repmat('\t%e',1,m-4+dim) , '\n'];
    end
    fprintf(fid,format,dim,level,size(data,1),'voigt',zeros(1,m-4+dim));
    format = [ repmat('%.10f\t',1,dim) , repmat('%e\t',1,m-1), '%e\n'];
    fprintf(fid,format,data');
end
fclose(fid);

