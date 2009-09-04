function [ density] = createdens( nodes, elements, bVoid, circles, cylindres)

nelements = size(elements,1);

if nargin < 4
    ncirc = 0;
else
    ncirc = size(circles, 1);
end

if nargin < 5
    ncyl = 0;
else
    ncyl = size(cylindres, 1);
end

if bVoid == 1
    density = ones(nelements, 1);
else
    density = .0001*ones(nelements, 1);
end
    
for i=1:nelements
    
    centre_x = .125*(sum(nodes(elements(i,2:9), 2)));
    centre_y = .125*(sum(nodes(elements(i,2:9), 3)));
    centre_z = .125*(sum(nodes(elements(i,2:9), 4)));
    
    p = [centre_x centre_y centre_z];
    
    bWhite = 0;
    
    for j=1:ncirc
       centre = circles(j,1:3);
       radius = circles(j,4);
        
       if(norm (centre-p) < radius+1.0e-4)
           bWhite = 1;
           break;
       end
    end

    if bWhite == 0
        for j=1:ncyl
           centre1 = cylindres(j,1:3);
           centre2 = cylindres(j,4:6);
           radius = cylindres(j,7);

           lambda = RelativeCoordinate (centre1, centre2, p);
           
           if lambda < 0 || lambda > 1
               continue;
           end
           
           centre = (1-lambda)*centre1 + (lambda)*centre2;
           
           if(norm (centre-p) < radius+1.0e-4)
               bWhite = 1;
               break;
           end
        end
    end
    
    if bWhite == 1
        if bVoid == 1
            density(i) = 1.0e-4;
        else
            density(i) = 1.0;
        end
    end
    
end