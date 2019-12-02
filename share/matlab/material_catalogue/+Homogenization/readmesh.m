function [nodes,elems,BCs,regions] = readmesh(meshfile)

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
    
    elems = [];
    regions = [];

    % Read 2D elements
    if numTriangle > 0
        elems = fscanf(fid,'%d 4 3 %d\n%d %d %d',[5,numTriangle]);
        regions = elems(2,:)';
        elems = elems(3:end,:)';
    end
    if numTriangleQuad > 0
        elems = fscanf(fid,'%d 5 3 %d\n%d %d %d',[5,numTriangleQuad]);
        regions = elems(2,:)';
        elems = elems(3:end,:)';
    end
    if numQuadr > 0
        elems = fscanf(fid,'%d 6 4 %d\n%d %d %d %d',[6,numQuadr]);
        regions = elems(2,:)';
        elems = elems(3:end,:)';
    end

    % Jump to BC section
    while isempty( strfind(line,'[Node BC]')) && ~feof(fid)
        line = fgetl(fid);
    end
    
    % Read BCs
    % bcs may be decoded in numbers or named like in bc_strings
    bc_strings = {'border','controlpoints','bottom','up','left','right','front','back'};
    
    pos_bc = ftell(fid);
    % node names: 1 border, 2 controlpoints, 3 bottom, 4 top, 5 left, 6 right
    fgetl(fid);
    %TODO: FIX
    BCs = fscanf(fid,'%8d nodes%d\n',[2,numNodeBC]);
    
    % : bottom, up, left, right, front, back
    if size(BCs,2) ~= numNodeBC
        BCs = [];
        for ii=1:length(bc_strings)
            fseek(fid,pos_bc,'bof');
            pos = ftell(fid);
            line = fgetl(fid);
            while isempty( strfind(line, bc_strings{ii})) && ~feof(fid)
                pos = ftell(fid);
                line = fgetl(fid);
            end
            fseek(fid,pos,'bof');
            [bc_nodes, count] = fscanf(fid,['%8d ',bc_strings{ii},'\n']);
            % Check if we read too far
            line = fgetl(fid);
            if any( strcmp(bc_strings, line))
                bc_nodes = bc_nodes(1:end-1);
                count = count - 1;
            end
            bc_nodes = [bc_nodes'; ii*ones(1,count)];
            BCs = [BCs, bc_nodes];
        end
    end
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
    regions = [];
    
end
