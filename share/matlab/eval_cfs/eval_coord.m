function [eval_dim, eval_coord_r_index, eval_coord_r] = eval_coord(eval_type,dimension,coords,node_coords,node_connectivity)
%EVAL_COORD evaluate the nearest existing grid-point and extract all nodes
%necessary
switch eval_type
    case 'point'
        %evaluate a single point
        eval_dim = 1;
        if length(coords) == 3
            eval_coords = coords;
        else
            error('Wrong coordinate specification!')
        end
        eval_plane = 0;
    case 'xy'
        %evaluate a cut through the domain (only applicable for 3D meshes!)
        eval_dim = 2;
        if dimension == 3
            if length(coords) == 3
                eval_coords = coords(3);
            elseif length(coords) == 1
                eval_coords = coords;
            else
                error('Wrong coordinate specification!')
            end
        elseif dimension == 2
            warning('A 2D-cut through a 2D-mesh is not meaningful, displaying the whole plane...')
            eval_coords = 0;
        else
            error('Faulty dimension!')
        end
        eval_plane = 1;
    case 'yz'
        %evaluate a cut through the domain (only applicable for 3D meshes!)
        eval_dim = 2;
        if dimension == 3
            if length(coords) == 3
                eval_coords = coords(1);
            elseif length(coords) == 1
                eval_coords = coords;
            else
                error('Wrong coordinate specification!')
            end
        elseif dimension == 2
            warning('A 2D-cut through a 2D-mesh is not meaningful, displaying the whole plane...')
            eval_coords = 0;
        else
            error('Faulty dimension!')
        end
        eval_plane = 2;
    case 'zx'
        %evaluate a cut through the domain (only applicable for 3D meshes!)
        eval_dim = 2;
        if dimension == 3
            if length(coords) == 3
                eval_coords = coords(2);
            elseif length(coords) == 1
                eval_coords = coords;
            else
                error('Wrong coordinate specification!')
            end
        elseif dimension == 2
            warning('A 2D-cut through a 2D-mesh is not meaningful, displaying the whole plane...')
            eval_coords = 0;
        else
            error('Faulty dimension!')
        end
        eval_plane = 3;
end

switch eval_plane
    case 0
        %find the gridpoint nearest to the specified point
        eval_coord_r_index = knnsearch(node_coords,eval_coords);
        eval_coord_r = node_coords(eval_coord_r_index,:);
        disp('New evaluation-point:')
        disp(['x: ' num2str(eval_coord_r(1))])
        disp(['y: ' num2str(eval_coord_r(2))])
        disp(['z: ' num2str(eval_coord_r(3))])
        disp('------------------------------------------------------------')
    case 1
        %xy-plane
        %find the gridpoint nearest to the defined one which enables us
        %to define a plane, where every point is exactly lying on it
        eval_coord_r_index_main = knnsearch(node_coords(:,3),eval_coords);
%       --TEST CODE START--
%             search for one of the elements containing the node
%             [element_index, ~] = find(node_connectivity==eval_coord_r_index,1);
%             compute the maximum distance normal to the plane to get the
%             element-size
%             [~, elements] = size(node_connectivity);
%             z_diff = zeros(elements,1);
%             node_list_element = node_connectivity(element_index,:);
%             for uu=1:elements
%                 z_diff(uu) = node_coords(node_list_element(uu),3)-...
%                     node_coords(eval_coord_r_index,3);
%             end
%             h_element = max(abs(z_diff));
%       --TEST CODE END--
        eval_coord_r_main = node_coords(eval_coord_r_index_main,:);
        eval_coord_r_index = find(node_coords(:,3)==eval_coord_r_main(1,3));
        eval_coord_r = node_coords(eval_coord_r_index,:);
        disp('New evaluation-coordinate:')
        disp(['z: ' num2str(eval_coord_r_main(3))])
        disp('------------------------------------------------------------')
    case 2
        %yz plane
        %find the gridpoint nearest to the defined one which enables us
        %to define a plane, where every point is exactly lying on it
        eval_coord_r_index_main = knnsearch(node_coords(:,1),eval_coords);
        eval_coord_r_main = node_coords(eval_coord_r_index_main,:);
        eval_coord_r_index = find(node_coords(:,1)==eval_coord_r_main(1,1));
        eval_coord_r = node_coords(eval_coord_r_index,:);
        disp('New evaluation-coordinate:')
        disp(['x: ' num2str(eval_coord_r_main(1))])
        disp('------------------------------------------------------------')
    case 3
        %zx plane
        %find the gridpoint nearest to the defined one which enables us
        %to define a plane, where every point is exactly lying on it
        eval_coord_r_index_main = knnsearch(node_coords(:,2),eval_coords);
        eval_coord_r_main = node_coords(eval_coord_r_index_main,:);
        eval_coord_r_index = find(node_coords(:,2)==eval_coord_r_main(1,2));
        eval_coord_r = node_coords(eval_coord_r_index,:);
        disp('New evaluation-coordinate:')
        disp(['y: ' num2str(eval_coord_r_main(2))])
        disp('------------------------------------------------------------')
end
end