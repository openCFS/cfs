function [ nodes, elements ] = readmesh( filename )

mesh_file = fopen(filename, 'r');

nnodes = 0;
nelements = 0;

while 1
  line = fgetl(mesh_file);
  text = sscanf(line, '%s', 1);

  if(strcmp(text, 'nodes') == 1)
    line = fgetl(mesh_file);
    val = sscanf(line, '%g', 1);
    nnodes = val;
    continue;
  end

  if(strcmp(text, 'elements') == 1)
    line = fgetl(mesh_file);
    val = sscanf(line, '%g', 1);
    nelements = val;
    break;
  end
end

if(nnodes == 0 || nelements == 0)
  disp('could not read number of elements and nodes');
  fclose(mesh_file);
  exit;
end

nodes=zeros(nnodes, 4);
elements=zeros(nelements, 9);

bDone = 0;

count = 0;

while bDone == 0
  line = fgetl(mesh_file);
  
  text = sscanf(line, '%s', 1);
  
  % Skip empty line
  if length(text) == 0
    continue;
  % Break at end of file
  elseif strcmp(text, 'end elements') == 1
    disp('mesh reading complete')
    break;
  % ... the element section ...
  elseif strcmp(text,'Coordinates') == 1
        while 1
          line = fgetl(mesh_file);

          text = sscanf(line, '%s', 1);

          % Skip empty line
          if length(text) == 0
            continue;
          % Break?
          elseif strcmp(text, 'end') == 1
              disp('all coordinates read')
              break;
          else
            count = count + 1;
            val = sscanf(line, '%d %g %g %g', 4);
            nodes(count, :) = val';
          end
        end
 
  % ... the node section ...
  elseif strcmp(text, 'Elements') == 1
    count = 0;
        while 1
          line = fgetl(mesh_file);
          text = sscanf(line, '%s', 1);

          % Skip empty line
          if length(text) == 0
            continue;
          % Break?
          elseif strcmp(text, 'end') == 1
              bDone = 1;
              break;
          else
              count = count + 1;
              val = sscanf(line, '%d %d %d %d %d %d %d %d %d', 9);
              elements(count, :) = val';
          end
        end
  end
  
end

fclose(mesh_file);

