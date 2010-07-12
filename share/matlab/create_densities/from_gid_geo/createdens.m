function [ density] = createdens( nodes, elements, bVoid, circles, cylinders )
disp('create density information');

nelements = size(elements,1);

if nargin < 4
    ncirc = 0;
else
    ncirc = size(circles, 1);
end

if nargin < 5
    ncyl = 0;
else
    ncyl = size(cylinders, 1);
end

if bVoid == 1
    density = ones(nelements, 1);
else
    density = .00001*ones(nelements, 1);
end

if( ncirc > 0 )
  for j=1:ncirc
    center = circles(j,1:3);
    radius = circles(j,4) + 1.0e-4;
          
    for i=1:nelements
      % nodes is the node number and then three coords
      % elements is the element number and then eight nodes
      center_x = .125*(sum(nodes(elements(i,2:9), 2)));
      center_y = .125*(sum(nodes(elements(i,2:9), 3)));
      center_z = .125*(sum(nodes(elements(i,2:9), 4)));
        
      p = [center_x center_y center_z];
        
      bWhite = 0;
      
      if(norm (center-p) < radius)
        bWhite = 1;
      end

      if bWhite == 1
          if bVoid == 1
              density(i) = 1.0e-5;
          else
              density(i) = 1.0;
          end
      end
      
  end
  end
end % end of circles

    
if( ncyl > 0 )
  fprintf('found %d cylinders\n', ncyl);
  cylinders2 = zeros(ncyl, 5);

  % walk over all cylinders
  % precalculate values for later
  for j=1:ncyl
    center1 = cylinders(j,1:3);
    center2 = cylinders(j,4:6);
    radius = cylinders(j,7);

    c1minc2 = center1 - center2;
    c1minc2norm = radius/norm(c1minc2);
    c1minc2norm2 = 1./norm(c1minc2)^2;

    cylinders2(j, :) = [c1minc2 c1minc2norm2 c1minc2norm];
  end

  % walk over all elements
  for i=1:nelements
    bWhite = 0;
            
    % nodes is the node number and then three coords
    % elements is the element number and then eight nodes
    center_x = .125*(sum(nodes(elements(i,2:9), 2)));
    center_y = .125*(sum(nodes(elements(i,2:9), 3)));
    center_z = .125*(sum(nodes(elements(i,2:9), 4)));
    
    n = [center_x center_y center_z];
      
    % walk over all cylinders
    for j=1:ncyl
      center1 = cylinders(j,1:3);
      center2 = cylinders(j,4:6);
      radius = cylinders(j,7) + 1.0e-4;

      lambda = dot(center1 - n, cylinders2(j, 1:3)) * cylinders2(j, 4);

      if lambda <  -cylinders2(j, 5) || lambda > 1+cylinders2(j, 5)
        continue;
      end

      if lambda < 0 || lambda > 1
        % set to 0 or 1 for the caps
        if lambda < 0
          lambda = 0;
        else
          lambda = 1;
        end
        % put continue here if you want to ignore the caps
        % continue;
      end

      % d is now the projection of the element center
      % to the line which defines the cylinder
      d = (1-lambda)*center1 + lambda*center2;
            
      % calculate the distance of the element center to the
      % projection onto the line
      if(norm (d-n) < radius)
        bWhite = 1;
        break;
      end

    end % for over cylinders

    if bWhite == 1
      if bVoid == 1
        density(i) = 1.0e-5;
      else
        density(i) = 1.0;
      end
    end

  end % for over elements
      
end % if ncyl > 0
