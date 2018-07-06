function [varargout] = scatter_to_grid(varargin)
%SCATTER_TO_GRID (V1.0) this function takes a cell array consisting of scattered data and
%interpolates it on a grid fitting the boundary
switch nargin
    case 1
        val_array = varargin{1};
        elements = 50;
    case 2
        val_array = varargin{1};
        elements = varargin{2};
end
    
%surpress warning
lastwarn('');
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');

if sum(isnan(val_array{1}(:,3)))>0
    %xy-plane
    grid_dim_id = 1;
    boundary_indeces = boundary(val_array{1}(:,1),val_array{1}(:,2));
    min_x = min(val_array{1,1}(:,1));
    min_y = min(val_array{1,1}(:,2));
    max_x = max(val_array{1,1}(:,1));
    max_y = max(val_array{1,1}(:,2));
    step_x = (max_x-min_x)/(elements-1);
    step_y = (max_y-min_y)/(elements-1);
    x_points = min_x:step_x:max_x;
    y_points = min_y:step_y:max_y;
    grid_vec = zeros(length(x_points)*length(y_points),2);
    for ii=1:length(x_points)
        for uu=1:length(y_points)
            grid_vec((ii-1)*length(x_points)+uu,1) = x_points(ii);
            grid_vec((ii-1)*length(x_points)+uu,2) = y_points(uu);
        end
    end
    in_indeces = inpolygon(grid_vec(:,1),grid_vec(:,2),val_array{1}(boundary_indeces,1),...
        val_array{1}(boundary_indeces,2));
    in_grid_vec = NaN(length(in_indeces),2);
    for ii=1:length(in_indeces)
        if in_indeces(ii)==1
            in_grid_vec(ii,1) = grid_vec(ii,1);
            in_grid_vec(ii,2) = grid_vec(ii,2);
        end
    end

    varargout{1} =  in_grid_vec(:,1);
    varargout{2} =  in_grid_vec(:,2);
    varargout{3} = NaN(length(in_indeces),1);
    
    reverseStr = '';
    length_val_array = length(val_array);
    grid_val_array = cell(length_val_array,1);
    grid_val_matrix_array = cell(length_val_array,1);
    for ii=1:length_val_array
        %create interpolation object
        int_obj = scatteredInterpolant(val_array{ii}(:,1),val_array{ii}(:,2),val_array{ii}(:,4));
        grid_val_vec = int_obj(in_grid_vec(:,1),in_grid_vec(:,2));
        sub_val_array = NaN(length(grid_val_vec),4);
        sub_val_array(:,1) = in_grid_vec(:,1);
        sub_val_array(:,2) = in_grid_vec(:,2);
        sub_val_array(:,4) = grid_val_vec;
        grid_val_array{ii} = sub_val_array;
        
        val_matrix = NaN(length(y_points),length(x_points));
        counter = 1;
        for uu=1:length(x_points)
            for oo=1:length(y_points)
                if in_indeces((uu-1)*length(x_points)+oo) == 1
                    % CHANGE
                    val_matrix(oo,uu) = grid_val_vec(counter);
                end
                counter = counter + 1;
            end
        end
        grid_val_matrix_array{ii} = val_matrix;
        
        percentDone = 100 * ii / length_val_array;
        msg = sprintf('Interpolating on grid: %3.1f', percentDone);
        fprintf([reverseStr, msg, '%%']);
        reverseStr = repmat(sprintf('\b'), 1, length(msg)+1); 
    end
    varargout{4} = grid_val_array;
    varargout{5} = grid_val_matrix_array;
    varargout{6} = grid_dim_id;
elseif sum(isnan(val_array{1}(:,1)))>0
    %yz-plane  
    grid_dim_id = 2;
    boundary_indeces = boundary(val_array{1}(:,2),val_array{1}(:,3));
    min_y = min(val_array{1,1}(:,2));
    min_z = min(val_array{1,1}(:,3));
    max_y = max(val_array{1,1}(:,2));
    max_z = max(val_array{1,1}(:,3));
    step_y = (max_y-min_y)/(elements-1);
    step_z = (max_z-min_z)/(elements-1);
    y_points = min_y:step_y:max_y;
    z_points = min_z:step_z:max_z;
    grid_vec = zeros(length(y_points)*length(z_points),2);
    for ii=1:length(y_points)
        for uu=1:length(z_points)
            grid_vec((ii-1)*length(y_points)+uu,1) = y_points(ii);
            grid_vec((ii-1)*length(y_points)+uu,2) = z_points(uu);
        end
    end
    in_indeces = inpolygon(grid_vec(:,1),grid_vec(:,2),val_array{1}(boundary_indeces,2),...
        val_array{1}(boundary_indeces,3));
    in_grid_vec = NaN(length(in_indeces),2);
    for ii=1:length(in_indeces)
        if in_indeces(ii)==1
            in_grid_vec(ii,1) = grid_vec(ii,1);
            in_grid_vec(ii,2) = grid_vec(ii,2);
        end
    end
    
    varargout{1} = NaN(length(in_indeces),1);
    varargout{2} =  in_grid_vec(:,1);
    varargout{3} =  in_grid_vec(:,2);
    
    reverseStr = '';
    length_val_array = length(val_array);
    grid_val_array = cell(length_val_array,1);
    grid_val_matrix_array = cell(length_val_array,1);
    for ii=1:length_val_array
        %create interpolation object
        int_obj = scatteredInterpolant(val_array{ii}(:,2),val_array{ii}(:,3),val_array{ii}(:,4));
        grid_val_vec = int_obj(in_grid_vec(:,1),in_grid_vec(:,2));
        sub_val_array = NaN(length(grid_val_vec),4);
        sub_val_array(:,2) = in_grid_vec(:,1);
        sub_val_array(:,3) = in_grid_vec(:,2);
        sub_val_array(:,4) = grid_val_vec;
        grid_val_array{ii} = sub_val_array;
        val_matrix = NaN(length(z_points),length(y_points));
        counter = 1;
        for uu=1:length(y_points)
            for oo=1:length(z_points)
                if in_indeces((uu-1)*length(y_points)+oo) == 1
                    % CHANGE
                    val_matrix(oo,uu) = grid_val_vec(counter);
                end
                counter = counter + 1;
            end
        end
        grid_val_matrix_array{ii} = val_matrix;
        
        percentDone = 100 * ii / length_val_array;
        msg = sprintf('Interpolating on grid: %3.1f', percentDone);
        fprintf([reverseStr, msg, '%%']);
        reverseStr = repmat(sprintf('\b'), 1, length(msg)+1); 
    end
    varargout{4} = grid_val_array;
    varargout{5} = grid_val_matrix_array;
    varargout{6} = grid_dim_id;
elseif sum(isnan(val_array{1}(:,2)))>0
    %zx-plane  
    grid_dim_id = 3;
    boundary_indeces = boundary(val_array{1}(:,1),val_array{1}(:,3),1);
    min_x = min(val_array{1,1}(:,1));
    min_z = min(val_array{1,1}(:,3));
    max_x = max(val_array{1,1}(:,1));
    max_z = max(val_array{1,1}(:,3));
    step_x = (max_x-min_x)/(elements-1);
    step_z = (max_z-min_z)/(elements-1);
    x_points = min_x:step_x:max_x;
    z_points = min_z:step_z:max_z;
    grid_vec = zeros(length(x_points)*length(z_points),2);
    for ii=1:length(x_points)
        for uu=1:length(z_points)
            grid_vec((ii-1)*length(x_points)+uu,1) = x_points(ii);
            grid_vec((ii-1)*length(x_points)+uu,2) = z_points(uu);
        end
    end
    in_indeces = inpolygon(grid_vec(:,1),grid_vec(:,2),val_array{1}(boundary_indeces,1),...
        val_array{1}(boundary_indeces,3));
    in_grid_vec = NaN(length(in_indeces),2);
    for ii=1:length(in_indeces)
        if in_indeces(ii)==1
            in_grid_vec(ii,1) = grid_vec(ii,1);
            in_grid_vec(ii,2) = grid_vec(ii,2);
        end
    end
    
    varargout{1} = grid_vec(:,1);
    varargout{2} = NaN(length(in_indeces),1);
    varargout{3} = grid_vec(:,2);
    
    reverseStr = '';
    length_val_array = length(val_array);
    grid_val_array = cell(length_val_array,1);
    grid_val_matrix_array = cell(length_val_array,1);
    for ii=1:length_val_array
        %create interpolation object
        int_obj = scatteredInterpolant(val_array{ii}(:,1),val_array{ii}(:,3),val_array{ii}(:,4),'linear','none');
        grid_val_vec = int_obj(in_grid_vec(:,1),in_grid_vec(:,2));
        sub_val_array = NaN(length(grid_val_vec),4);
        sub_val_array(:,1) = in_grid_vec(:,1);
        sub_val_array(:,3) = in_grid_vec(:,2);
        sub_val_array(:,4) = grid_val_vec;
        grid_val_array{ii} = sub_val_array;
        
        val_matrix = NaN(length(z_points),length(x_points));
        counter = 1;
        for uu=1:length(x_points)
            for oo=1:length(z_points)
                if in_indeces((uu-1)*length(x_points)+oo) == 1
                    % CHANGE
                    val_matrix(oo,uu) = grid_val_vec(counter);
                end
                counter = counter + 1;
            end
        end
        grid_val_matrix_array{ii} = val_matrix;
        
        percentDone = 100 * ii / length_val_array;
        msg = sprintf('Interpolating on grid: %3.1f', percentDone);
        fprintf([reverseStr, msg, '%%']);
        reverseStr = repmat(sprintf('\b'), 1, length(msg)+1); 
    end
    varargout{4} = grid_val_array;
    varargout{5} = grid_val_matrix_array;
    varargout{6} = grid_dim_id;
else
    %3D
    
end
fprintf('\n')
[warnMsg, ~] = lastwarn;
if ~isempty(warnMsg)
    disp('The following warning appeared during the calculations:')
    disp(warnMsg)
end