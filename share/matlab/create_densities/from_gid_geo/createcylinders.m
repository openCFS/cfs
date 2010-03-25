function cyl = createcylinders( nodes, edges, radius)
disp(['createcylinders with radius = ', num2str(radius)]);

nelements = size(edges,1);

cyl = zeros(nelements, 7);

for i=1:nelements
    
    p = nodes(edges(i,2),2:4);
    q = nodes(edges(i,3),2:4);
    
    cyl(i,:) = [p q radius];
end

