function [ nodes, elements ] = readmesh( filename )

mesh_file = fopen(filename, 'r');

nodes=[];
elements=[];

nelements = 0;
nnodes = 0;

bDone = 0;

while bDone == 0
  line = fgetl(mesh_file);
  
  text = sscanf(line, '%s', 1);
  
  % Skip empty line
  if length(text) == 0
    continue;
  % Break?
  elseif strcmp(text, 'end elements') == 1
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
              break;
          else
              nnodes = nnodes + 1;
              val = sscanf(line, '%d %g %g %g', 9);
              if length(nodes) == 0
                  nodes = val';
              else
                  nodes = [nodes;val'];
              end
          end
        end
 
  % ... the node section ...
  elseif strcmp(text, 'Elements') == 1
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
              nelements = nelements + 1;
              val = sscanf(line, '%d %d %d %d %d %d %d %d %d', 9);
              if length(elements) == 0
                  elements = val';
              else
                  elements = [elements;val'];
              end
          end
        end
  end
  
end

fclose(mesh_file);