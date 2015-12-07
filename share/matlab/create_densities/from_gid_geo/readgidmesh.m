function [ nodes, elements ] = readgidmesh( filename )

mesh_file = fopen(filename, 'r');

line = fgetl(mesh_file);
if strcmp(line, '[Info]') != 1
  disp('wrong file type, please provide GiD-mesh-file');
  fclose(mesh_file);
  exit;
end

nnodes = 0;
nelements = 0;

while 1
  line = fgetl(mesh_file);
  [tok, val] = strtok(line);

  if(strcmp(tok, 'NumNodes') == 1)
    nnodes = str2num(val);
    continue;
  end

  if(strcmp(tok, 'Num3DElements') == 1)
    nelements = str2num(val);
    break;
  end
end


if(nnodes == 0 || nelements == 0)
  disp('could not read number of elements and nodes');
  fclose(mesh_file);
  exit;
end

disp(['found ', num2str(nnodes), ' nodes and ', num2str(nelements), ' elements']);

nodes=zeros(nnodes, 4);
elements=zeros(nelements, 9);

bDone = 0;
count = 0;

while bDone == 0
  line = fgetl(mesh_file);
  
  % Skip empty line
  if isempty(line)
    continue;

  % Break at end of file
  elseif strcmp(line, '[Node BC]') == 1
    disp('mesh reading complete');
    bDone = 1;
    break;
  

  % ... the node section ...
  elseif strcmp(line, '[Nodes]') == 1
    disp('reading nodes');
    line = fgetl(mesh_file);
    count = 0;
    while count < nnodes
      line = fgetl(mesh_file);

      count = count + 1;
      val = sscanf(line, '%d %g %g %g', 4);
      nodes(count, :) = val';
    end
 
  % ... the element section ...
  elseif strcmp(line, '[3D Elements]') == 1
    disp('reading elements');
    line = fgetl(mesh_file);
    line = fgetl(mesh_file);
    count = 0;
    while count < nelements
      line = fgetl(mesh_file);
      line = fgetl(mesh_file);

      count = count + 1;
      val = sscanf(line, '%d %d %d %d %d %d %d %d', 8);
      elements(count, :) = [count val'];
    end
    bDone = 1;
  end
  
end

fclose(mesh_file);
