function elems = plotmesh(arg1,showNumbers,arg2)
% PLOTMESH plots an ANSYS mesh
%
%   nodes = plotmesh(file) plots the nodes given in file and returns the
%   defining nodes.
%   nodes = plotmesh(nodes) plots the nodes given in nodes and returns
%   them.
%   nodes = plotmesh(~,showNumbers) plots the nodes from above and their 
%   corresponding numbers. The function returns the defining nodes.
%   nodes = plotmesh(~,~,bcs) plots the nodes from above and their 
%   corresponding numbers. Additionally the labeled nodes are coloured.
%   The function returns the defining nodes.
%
%   Example:
%   plotmesh('./myMesh.mesh',1)
%

scale = 3;

if nargin < 2
    showNumbers = 0;
end

if ischar(arg1)
    [nodes,elems,bcs] = Homogenization.readmesh(arg1);
    [filepath,filename,ext] = fileparts(arg1);
    ismesh = strcmp(ext,'.mesh');
else
    nodes = arg1;
    bcs = arg2;
    ismesh = 1;
end

if ~ismesh % density file
    density = nodes;
    [nodes,elems,bcs] = Homogenization.readmesh([filepath,'/',filename,'.mesh']);
else
    density = ones(size(elems,1),1);
end

bcs = bcs(abs(bcs(:,2)-1)>1e-2,:);

% node names: 1 border, 2 controlpoints, 3 bottom, 4 top, 5 left, 6 right
normalNodes = setdiff(1:size(nodes,1),bcs(:,1));
controlNodes = bcs(abs(bcs(:,2)-2)<1e-2,1);
bottom = bcs(abs(bcs(:,2)-3)<1e-2,1);
top = bcs(abs(bcs(:,2)-4)<1e-2,1);
left = bcs(abs(bcs(:,2)-5)<1e-2,1);
right = bcs(abs(bcs(:,2)-6)<1e-2,1);

nodesToPlot = 1:min(size(nodes,1),1e100);
% nodesToPlot = 1:16640;
% nodesToPlot = 1:size(nodes,1);

figure;
hold on
if isempty(elems)
    normalNodesToPlot = normalNodes(normalNodes<=max(nodesToPlot));
    controlNodesToPlot = controlNodes(controlNodes<=max(nodesToPlot));
    bottomNodesToPlot = bottom(bottom<=max(nodesToPlot));
    topNodesToPlot = top(top<=max(nodesToPlot));
    leftNodesToPlot = left(left<=max(nodesToPlot));
    rightNodesToPlot = right(right<=max(nodesToPlot));
    scatter(nodes(normalNodesToPlot,1),nodes(normalNodesToPlot,2),scale*10,'MarkerEdgeColor',[0 0 0]);
    scatter(nodes(leftNodesToPlot,1),nodes(leftNodesToPlot,2),scale*10,'MarkerEdgeColor',[1 .2 1]);
    scatter(nodes(rightNodesToPlot,1),nodes(rightNodesToPlot,2),scale*10,'MarkerEdgeColor',[1 .8 0]);
    scatter(nodes(bottomNodesToPlot,1),nodes(bottomNodesToPlot,2),scale*10,'MarkerEdgeColor',[0 .8 0]);
    scatter(nodes(topNodesToPlot,1),nodes(topNodesToPlot,2),scale*10,'MarkerEdgeColor',[.2 .2 1]);
    scatter(nodes(controlNodesToPlot,1),nodes(controlNodesToPlot,2),scale*15,'MarkerEdgeColor',[1 0 0]);
else
    p = patch('Faces',elems,'Vertices',nodes);
    cdata = reshape(flipud(density)',size(density,1)*size(density,2),1);
    set(p,'FaceColor','flat','FaceVertexCData',cdata,'CDataMapping','scaled');
end

if showNumbers
    labels = strcat(num2str((1:2:2*nodesToPlot(end))'),',',num2str((2:2:2*nodesToPlot(end))'));
    text(nodes(nodesToPlot,1),nodes(nodesToPlot,2),labels);
end
xmin = min(nodes(:,1));
xmax = max(nodes(:,1));
ymin = min(nodes(:,2));
ymax = max(nodes(:,2));
axis equal
axis tight
xlim([xmin-(xmax-xmin)/100*scale,xmax+(xmax-xmin)/100*scale]);
ylim([ymin-(ymax-ymin)/100*scale,ymax+(ymax-ymin)/100*scale]);
hold off
