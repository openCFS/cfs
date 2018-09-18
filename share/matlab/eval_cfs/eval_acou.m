function varargout = eval_acou(varargin)
%EVAL_ACOU (V1.0) this function evaluates all results considering acoustics
%and is also able to apply postprocessing to the simulation results
switch nargin
    case 0
        %return the information needed by the main program
        info_cell = cell(3,1);
        info_cell{1} = 'acouPotential';
        info_cell{2} = 'acouPressure';
        info_cell{3} = 'acouVelocity';
        varargout{1} = info_cell;
        return
    case 1
        eval_type = varargin{1};
        switch eval_type
            case 'acouPotential'
                varargout{1} = 'Acoustic potential (m²s)';
            case 'acouPressure'
                varargout{1} = 'Acoustic pressure (Pa)';
            case 'acouVelocity'
                varargout{1} = 'Particle velocity (m/s)';
        end
        return
    case 5
        %evaluate results considering the chosen postprocessing
        result_type = varargin{1};
        postprocess = varargin{2};
        time_array = varargin{3};
        val_array = varargin{4};
        rho0 = varargin{5};
    case 6
        %plot the data, not yet implemented
        
        comp_data()
        
        plot_quiver();
        plot_quiver3();
        return
end

switch postprocess
    case 'acouPressure'
        %calculate the acoustic pressure from the acoustic potential
        if strcmp(result_type,'acouPotential') == 1
            dt = time_array(2)-time_array(1);
            post_time_array = time_array(1:end-1);
            post_val_array = cell(length(val_array)-1,1);
            if length(val_array{1}) == 4
                %point
                for ii=1:length(val_array)-1
                    sub_val_array = NaN(1,4);
                    sub_val_array(:,1) = val_array{ii,1}(:,1);
                    sub_val_array(:,2) = val_array{ii,1}(:,2);
                    sub_val_array(:,3) = val_array{ii,1}(:,3);
                    sub_val_array(:,4) = rho0*(val_array{ii+1,1}(1,4)-...
                        val_array{ii,1}(1,4))/dt;
                    post_val_array{ii,1} = sub_val_array;
                end
            end
            if sum(isnan(val_array{1}(:,3)))>0
                %xy-plane
                for ii=1:length(val_array)-1
                    sub_val_array = NaN(length(val_array{ii,1}),3);
                    sub_val_array(:,1) = val_array{ii,1}(:,1);
                    sub_val_array(:,2) = val_array{ii,1}(:,2);
                    sub_val_array(:,4) = rho0*(val_array{ii+1,1}(:,4)-...
                        val_array{ii,1}(:,4))/dt;
                    post_val_array{ii,1} = sub_val_array;
                end
            elseif sum(isnan(val_array{1}(:,1)))>0
                %yz-plane
                for ii=1:length(val_array)-1
                    sub_val_array = NaN(length(val_array{ii,1}),3);
                    sub_val_array(:,2) = val_array{ii,1}(:,2);
                    sub_val_array(:,3) = val_array{ii,1}(:,3);
                    sub_val_array(:,4) = rho0*(val_array{ii+1,1}(:,4)-...
                        val_array{ii,1}(:,4))/dt;
                    post_val_array{ii,1} = sub_val_array;
                end 
            elseif sum(isnan(val_array{1}(:,2)))>0
                %zx-plane
                for ii=1:length(val_array)-1
                    sub_val_array = NaN(length(val_array{ii,1}),3);
                    sub_val_array(:,1) = val_array{ii,1}(:,1);
                    sub_val_array(:,3) = val_array{ii,1}(:,3);
                    sub_val_array(:,4) = rho0*(val_array{ii+1,1}(:,4)-...
                        val_array{ii,1}(:,4))/dt;
                    post_val_array{ii,1} = sub_val_array;
                end      
            end
        else
            warning('No matching postprocessing method found, the program will continue without any postprocessing!')
        end
        varargout{1} = post_time_array;
        varargout{2} = post_val_array;
        varargout{3} = 'Soundpressure (Pa)';
        varargout{4} = ''; %if the post_eval_dim field is left out, the eval_dim from the main file will be used
        varargout{5} = '';
    case 'acouVelocity'
        %calculate the particle velocity from the acoustic potential
        if strcmp(result_type,'acouPotential') == 1
            %construct a equally spaced grid from the arbitrary input-grid
            [x_vec, y_vec, z_vec, grid_val_array, grid_val_matrix, grid_dim_id] = ...
                scatter_to_grid(val_array);
            disp(newline)
        else
            %wrong input, can´t compute the particle velocity
        end
        length_val_array = length(val_array);
        grad_array = cell(length_val_array,1);
        %compute the gradient depending on the plane
        switch grid_dim_id
            case 1
                %xy-plane
                for ii=1:length_val_array
                    grad_cell = cell(1,3);
                    [grad_cell{1,1}, grad_cell{1,2}] = gradient(grid_val_matrix{ii});
                    grad_cell{1,1} = -grad_cell{1,1};
                    grad_cell{1,2} = -grad_cell{1,2};
                    grad_array{ii} = grad_cell;
                end
            case 2
                %yz-plane
                for ii=1:length_val_array
                    grad_cell = cell(1,3);
                    [grad_cell{1,2}, grad_cell{1,3}] = gradient(grid_val_matrix{ii});
                    grad_cell{1,2} = -grad_cell{1,2};
                    grad_cell{1,3} = -grad_cell{1,3};
                    grad_array{ii} = grad_cell;
                end
            case 3
                %zx-plane
                for ii=1:length_val_array
                    grad_cell = cell(1,3);
                    [grad_cell{1,1}, grad_cell{1,3}] = gradient(grid_val_matrix{ii});
                    grad_cell{1,1} = -grad_cell{1,1};
                    grad_cell{1,3} = -grad_cell{1,3};
                    grad_array{ii} = grad_cell;
                end
        end
        post_time_array = time_array;
        varargout{1} = post_time_array;
        varargout{2} = grad_array;
        varargout{3} = 'Particle velocity (m/s)';
        varargout{4} = '';
        switch grid_dim_id
            case 1
                cont_out_cell{1} = 'xy';
                cont_out_cell{2} = x_vec;
                cont_out_cell{3} = y_vec;
                cont_out_cell{4} = z_vec;
                varargout{5} = cont_out_cell;
            case 2
                cont_out_cell{1} = 'yz';
                cont_out_cell{2} = x_vec;
                cont_out_cell{3} = y_vec;
                cont_out_cell{4} = z_vec;
                varargout{5} = cont_out_cell;
            case 3
                cont_out_cell{1} = 'zx';
                cont_out_cell{2} = x_vec;
                cont_out_cell{3} = y_vec;
                cont_out_cell{4} = z_vec;
                varargout{5} = cont_out_cell;
        end
end

end