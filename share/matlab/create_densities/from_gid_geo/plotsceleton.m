function plotsceleton( nodes, edges)

nelements = size(edges,1);

figure;
hold on;

for i=1:nelements
    
    p = nodes(edges(i,2),2:4);
    q = nodes(edges(i,3),2:4);
    
    plot3([p(1) q(1)], [p(2) q(2)], [p(3) q(3)], 'k-*');
    
end

axis equal;
