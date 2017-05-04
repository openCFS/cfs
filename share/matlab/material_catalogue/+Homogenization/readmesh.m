function [nodes,elems,BCs] = readmesh(meshfile)

[~,~,ext] = fileparts(meshfile);

if strcmp(ext,'.mesh')

    fid = fopen(meshfile,'r');
    
    for i=1:4
        line = fgetl(fid);
    end

    numNodes = str2double(line(10:end));

    fgetl(fid);

    line = fgetl(fid);
    num2DElems = str2double(line(15:end));

    fgetl(fid);

    line = fgetl(fid);
    numNodeBC = str2double(line(11:end));

    for i=1:7
        line = fgetl(fid);
    end
    numTriangle = str2double(line(20:end));
    line = fgetl(fid);
    numTriangleQuad = str2double(line(20:end));
    line = fgetl(fid);
    numQuadr = str2double(line(20:end));

    % Jump to node section
    while isempty(strfind(line,'[Nodes]')) && ~feof(fid)
        line = fgetl(fid);
    end
    fgetl(fid);

    % Read nodes
    % nodes = zeros(numNodes,2);
    % for i=1:numNodes
    %     line = fgetl(fid);
    %     node = str2num(line);
    %     nodes(i,:) = node(2:3);
    % end
    nodes = fscanf(fid,'%8d %e %e %e',[4,numNodes]);
    nodes = nodes(2:3,:)';

    % Jump to 2D elements section
    while isempty(strfind(line,'[2D Elements]')) && ~feof(fid)
        line = fgetl(fid);
    end
    fgetl(fid);
    fgetl(fid);

    % Read 2D elements
    if numTriangle > 0
        elems = fscanf(fid,'%d 4 3 mech\n%d %d %d',[4,numTriangle]);
        elems = elems(2:4,:)';
    end
    if numTriangleQuad > 0
        elems = fscanf(fid,'%d 5 3 mech\n%d %d %d',[4,numTriangleQuad]);
        elems = elems(2:4,:)';
    end
    if numQuadr > 0
        elems = fscanf(fid,'%d 6 4 mech\n%d %d %d %d',[5,numQuadr]);
        elems = elems(2:5,:)';
    %     elems = fscanf(fid,'%d 6 4 region%d\n %d %d %d %d',[6,num2DElems]);
    %     elems = elems(3:6,:)';
    end

    % Jump to BC section
    while isempty(strfind(line,'[Node BC]')) && ~feof(fid)
        line = fgetl(fid);
    end
    fgetl(fid);

    % Read BCs
    BCs = fscanf(fid,'%8d nodes%d\n',[2,numNodeBC]);
    BCs = BCs';
    
    fclose(fid);

else % density file
    
    fid = fopen(meshfile,'r');
    
    line = '';

    while isempty(strfind(line,'<mesh')) && ~feof(fid)
        line = fgetl(fid);
    end
    C = strsplit(line,'"');
    n = str2double(C{2});
    m = str2double(C{4});

    while isempty(strfind(line,'<set')) && ~feof(fid)
        line = fgetl(fid);
    end
    
    A = fscanf(fid,'    <element nr="%d" type="density" design="%e"/>\n',[2,n*m]);
    
    fclose(fid);
    
    nodes = reshape(A(2,:),n,m);
    nodes = flipud(nodes');
    elems = [];
    BCs = [];
    
end
