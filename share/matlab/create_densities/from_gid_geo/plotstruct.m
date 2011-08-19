function plotstruct( nodes, elements, density)

nelements = size(elements,1);

figure;
hold on;

for i=1:nelements
    
    centre_x = .125*(sum(nodes(elements(i,2:9), 2)));
    centre_y = .125*(sum(nodes(elements(i,2:9), 3)));
    centre_z = .125*(sum(nodes(elements(i,2:9), 4)));
    
    p = [centre_x centre_y centre_z];
    
    if density(i) > .5
        plot3(centre_x, centre_y, centre_z, 'g*');
    else
        plot3(centre_x, centre_y, centre_z, 'b*');
    end
    
end