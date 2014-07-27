function [ nodes, edges ] = readgeo( filename )

geo_file = fopen(filename, 'r');

nodes=[];
edges=[];

line = fgetl(geo_file);
line = fgetl(geo_file);
line = fgetl(geo_file);
line = fgetl(geo_file);
line = fgetl(geo_file);
line = fgetl(geo_file);
line = fgetl(geo_file);
  

bDone = 0;
bNodes = 0;
bEdges = 0;

while bDone ~= 1
    line = fgetl(geo_file);
    text = sscanf(line, '%s', 1);
 
    if strcmp(text,'1') == 1 && bEdges == 0
        bNodes = 1;
        val = sscanf(line, '%d %d', 2);
        nodeidx = val(2);
        line = fgetl(geo_file);
        val = sscanf(line, '%g %g %g', 3);
        if length(nodes) == 0
            nodes = [nodeidx val(1) val(2) val(3)];
        else
            nodes = [nodes; nodeidx val(1) val(2) val(3)];
        end
    elseif strcmp(text,'2') == 1
        bEdges = 1;
        val = sscanf(line, '%d %d', 2);
        edgeidx = val(2);
        line = fgetl(geo_file);
        val = sscanf(line, '%d %d', 2);
        if length(edges) == 0
            edges = [edgeidx val(1) val(2)];
        else
            edges = [edges; edgeidx val(1) val(2)];
        end
    elseif bEdges == 1
        bDone = 1;
    end
end
  

fclose(geo_file);
