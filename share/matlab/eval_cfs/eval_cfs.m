function varargout = eval_cfs(f_path,varargin)
%EVAL_H5 (V1.0) Evaluate .cfs files.
%   EVAL_H5(f_path,result_type,multi_step,eval_type,coords,postprocess,return_procedure,return_procedure_post)
%   evaluates a .cfs result file from a TRANSIENT CFS++ simulation. The input
%   parameters allow the user to specify the point or plane at which the
%   file will be evaluated as well as different types of postprocessing and
%   returning the results.
%   eval_h5(f_path) reads the file and only prints the available results
%   from all multisteps
%   eval_h5(f_path,result_type,multi_step,eval_type,coords,postprocess,return_procedure)
%   eval_h5(f_path,result_type,multi_step,eval_type,coords,postprocess,return_procedure,return_procedure_post)  
%
%   arguments: (input)
%   "f_path":                   'C:\...\file.cfs'
%   "result_type":              'acouPotential'
%   "eval_type":                'point' 'xy' 'yz' 'zx'
%   "coords":                   [double double double]
%   "postprocess":              'acouPressure'
%   "return_procedure":         'plot' 'ret_array'
%   "return_procedure_post":    'plot' 'ret_array'
%
%   Currently only the first multistep is implemented, but it can be easily
%   extended to cover more than one

%% check if the file exists and extract its name
try 
    info = h5info(f_path);
    path_name = info.Filename;
    filename_index = strfind(path_name,'\');
    filename = path_name(filename_index(end)+1:end);
    disp('------------------------------------------------------------')
    disp(['Evaluating ' filename])
    disp('------------------------------------------------------------')
    %debug
    %disp(size(varargin))
catch
    error('The file does not exist or is not supported!')
end

%% load evaluation-files
% invoke all user-made applications based on the evaluation-bare-bone
% (seperate items with a semi-colon!!)
eval_list_string = ['eval_acou'; 'eval_mech'];
eval_list = {};
eval_file_list = {};
[n,~] = size(eval_list_string);
for ii=1:n
    eval_dummy = feval(eval_list_string(ii,:));
    eval_list = [eval_list; eval_dummy];
    % create a second array where the type of function is stored which can
    % create the desired result considering postprocessing
    for uu=1:length(eval_dummy)
        eval_file_list = [eval_file_list; eval_list_string(ii,:)];
    end
end

%% define input parameters
%default values (just for testing, should be overwritten with actual values
%from the input!)
multi_step = 1;
eval_type = 'acouPotential';
coords = [0 0 0];
postprocess = 'acouPressure';

[~, m_v] = size(varargin);
switch m_v
    case 0
        %print only the result-type
        %scan how many multisteps have been executed
        multi_info = h5info(f_path,'/Results/Mesh');
        multi_info = multi_info.Groups;
        disp('The following results are available:')
        for uu=1:length(multi_info)
            h5_info = h5info(f_path,['/Results/Mesh/MultiStep_' num2str(uu) '/Step_1']);
            h5_info = h5_info.Groups;
            for ii=1:length(h5_info)
                info_string = h5_info(ii).Name;
                info_index = strfind(info_string,'/');
                h5_info_eval{ii,uu} = info_string(info_index(end)+1:end);
            end
            disp(['MultiStep_' num2str(uu)])
            disp(h5_info_eval)
        end
        disp('Results which can be processed:')
        disp(eval_list)
        return
    case 1
        
    case 2
        
    case 6
        result_type = varargin{1};
        multi_step = varargin{2};
        eval_type = varargin{3};
        coords = varargin{4};
        postprocess = varargin{5};
        return_procedure = varargin{6};
        return_procedure_post = return_procedure;
    case 7
        result_type = varargin{1};
        multi_step = varargin{2};
        eval_type = varargin{3};
        coords = varargin{4};
        postprocess = varargin{5};
        return_procedure_post = varargin{7};
end

%% extract the grid information and prepare for evaluation
dimension = h5readatt(f_path,'/Mesh','Dimension');
node_coords=h5read(f_path,'/Mesh/Nodes/Coordinates');
node_coords = transpose(node_coords);

node_coords(:,1:3) = round(node_coords(:,1:3),12);

node_connectivity = h5read(f_path,'/Mesh/Elements/Connectivity');
node_connectivity = transpose(node_connectivity);

%retrieve the information from the input parameters considering type of the
%evaluation and the coordinates
[eval_dim, eval_coord_r_index, eval_coord_r] = eval_coord(eval_type,dimension,coords,node_coords,node_connectivity);


%% extract the results and process them based on the grid information

%extract the regions which contain results
disp('Extracting regions...')
res_regions_info = h5info(f_path,['/Results/Mesh/MultiStep_1/Step_1/' result_type]);
res_regions_info = res_regions_info.Groups;

disp([num2str(length(res_regions_info)) ' regions found!'])
res_regions = cell(length(res_regions_info),1);
for ii=1:length(res_regions_info)
    path = res_regions_info(ii).Name;
    region_index = strfind(path,'/');
    region = path(region_index(end)+1:end);
    res_regions{ii} = region;
    disp(region)
end

disp('------------------------------------------------------------')
disp('Evaluating results')
disp('------------------------------------------------------------')

%retrieve the unit of the source-data
eval_unit = h5read(f_path,['/Results/Mesh/MultiStep_' ...
    num2str(multi_step) '/ResultDescription/' result_type '/Unit']);
eval_unit = eval_unit{1};

switch eval_dim
    case 1
        %evaluate a single point
        
        %retrieve the number of timesteps from the .h5 file
        MS_info = h5info(f_path,'/Results/Mesh/MultiStep_1');
        time_steps = length(MS_info.Groups)-1;
        disp([num2str(time_steps) ' time-steps found!'])
        %search for the region which contains the node
        for ii=1:length(res_regions_info)
            nodes = h5read(f_path,['/Mesh/Regions/' res_regions{ii} '/Nodes']);
            for uu=1:length(nodes)
                if nodes(uu) == eval_coord_r_index
                    eval_region = res_regions{ii};
                    node_index = uu;
                    disp(['Evaluating node ' num2str(nodes(uu))])
                    break
                end
            end
            if nodes(uu) == eval_coord_r_index
                break
            end
        end
        
        disp('Extracting values of the time-steps')
        [n_nva, ~] = size(h5read(f_path,...
            ['/Results/Mesh/MultiStep_1/Step_1/' num2str(result_type)...
            '/' eval_region '/Nodes/Real']));
        sub_val_array = NaN(1,3+n_nva);
        val_array_plot = zeros(time_steps,1);
        val_array = cell(time_steps,1);
        time_array = zeros(time_steps,1);
        reverseStr = '';
        for ii=1:time_steps
            time_string = ['/Results/Mesh/MultiStep_1/Step_' num2str(ii)];
            time_array(ii,1) = h5readatt(f_path,time_string,'StepValue');
            val_string = ['/Results/Mesh/MultiStep_1/Step_' num2str(ii) '/' ...
                num2str(result_type) '/' eval_region '/Nodes/Real'];
            node_val_array = h5read(f_path,val_string);
            sub_val_array(1,1:3) = eval_coord_r;
            sub_val_array(1,4:end) = node_val_array(:,node_index);
            val_array{ii,1} = sub_val_array;
            val_array_plot(ii,1) = node_val_array(node_index);
            percentDone = 100 * ii / time_steps;
            msg = sprintf('Reading results: %3.1f', percentDone);
            fprintf([reverseStr, msg, '%%']);
            reverseStr = repmat(sprintf('\b'), 1, length(msg)+1); 
        end
        disp([newline 'Finished!'])

        %return results
        if strcmp(return_procedure,'plot') == 1
            disp('Plotting results')
            figure
            plot(time_array(:,1),val_array_plot(:,1))
        elseif strcmp(return_procedure,'ret_array') == 1
            varargout{1} = time_array;
            varargout{2} = val_array;
            varargout{3} = [result_type ' (' eval_unit ')'];
            if isempty(postprocess)==1
                return 
            end 
        else
            warning('No valid return-procedure found, the default procedure will be applied (plot)')
            figure
            plot(time_array(:,1),val_array_plot(:,1))
        end
    case 2
        %evaluate a plane
        reverseStr = '';
        %retrieve the number of timesteps from the .h5 file
        MS_info = h5info(f_path,'/Results/Mesh/MultiStep_1');
        time_steps = length(MS_info.Groups)-1;
        disp([num2str(time_steps) ' time-steps found!'])
        
        %search for the regions which contain the nodes
        disp('Extracting node-indeces and creating dependency-matrix')
        eval_region = cell(length(eval_coord_r_index),1);
        node_index = zeros(length(eval_coord_r_index),1);
        for ii=1:length(res_regions_info)
            reverseStr = '';
            nodes = h5read(f_path,['/Mesh/Regions/' res_regions{ii} '/Nodes']);
            for oo=1:length(eval_coord_r_index)
                for uu=1:length(nodes)
                    if nodes(uu) == eval_coord_r_index(oo)
                        eval_region{oo} = res_regions{ii};
                        node_index(oo) = uu;
                        break
                    end
                end
                percentDone = 100 * oo / length(eval_coord_r_index);
                msg = ['Reading indeces ', num2str(res_regions{ii})];
                msg = [msg sprintf(': %3.1f', percentDone)];
                fprintf([reverseStr, msg, '%%']);
                reverseStr = repmat(sprintf('\b'), 1, length(msg)+1); 
            end
            fprintf('\n')
        end
        
        %check for "fake" datapoints
        nan_vec = check_data(f_path,eval_coord_r_index,...
            eval_region,result_type,node_index,eval_type,eval_coord_r);
        counter = 1;
        new_vec_length = length(eval_coord_r_index)-sum(nan_vec);
        eval_coord_r_index_new = zeros(new_vec_length,1);
        eval_region_new = cell(new_vec_length,1);
        node_index_new = zeros(new_vec_length,1);
        for ii=1:length(eval_coord_r_index)
            if nan_vec(ii,1) == 0
                eval_coord_r_index_new(counter,1) = eval_coord_r_index(ii,1);
                eval_region_new{counter,1} = eval_region{ii,1};
                node_index_new(counter,1) = node_index(ii,1);
                counter = counter + 1;
            end
        end
        eval_coord_r_new = node_coords(eval_coord_r_index_new,:);
        
        
        disp('Extracting values of the time-steps')
        %define a size limit when to start executing the reading-process in
        %parallel
        size_limit = 500e+06;
        val_array = cell(time_steps,1);
        time_array = zeros(time_steps,1);
        reverseStr = '';
        file_info = dir(f_path);
        file_size = file_info.bytes;
        
        if file_size < size_limit
            %do not execute the function in a parallel-pool
            %execute the function in a parallel-pool
            disp('Requested file is within the size-limit, using standard computation')
            reverseStr = '';
            fprintf('Reading results: ');
            for ii=1:time_steps
                array_buffer = read_val_array(ii,...
                    f_path,eval_coord_r_index_new,eval_region_new,result_type,...
                    node_index_new,eval_type,eval_coord_r_new);
                time_array(ii,1) = array_buffer{1};
                val_array{ii,1} = array_buffer{2};
                percentDone = 100 * ii / time_steps;
                msg = sprintf('%3.1f', percentDone);
                fprintf([reverseStr, msg, '%%']);
                reverseStr = repmat(sprintf('\b'), 1, length(msg)+1); 
            end
            fprintf('\n');
        else
            %execute the function in a parallel-pool
            disp('Requested file exceeds size-limit, using parallel computation')
            p_pool = gcp('nocreate'); % If no pool, do not create new one.
            if isempty(p_pool)
                disp('No active parallel pool found, starting a new one')
            else
                poolsize = p_pool.NumWorkers;
                disp(['Using active parallel pool ... connected to ' ...
                    num2str(poolsize) ' workers'])
            end
            p_pool = gcp();
            for ii=1:time_steps
                func_inst(ii) = parfeval(p_pool,@read_val_array,1,ii,f_path,...
                    eval_coord_r_index_new,eval_region_new,result_type,node_index_new,...
                    eval_type,eval_coord_r_new);
            end
            % Collect the results as they become available
            reverseStr = '';
            counter = 1;
            fprintf('Reading results: ');
            for ii=1:time_steps
                [complete_index,array_out] = fetchNext(func_inst);
                time_array(complete_index,1) = array_out{1};
                val_array{complete_index,1} = array_out{2};
                percentDone = 100 * counter / time_steps;
                counter = counter + 1;
                msg = sprintf('%3.1f', percentDone);
                fprintf([reverseStr, msg, '%%']);
                reverseStr = repmat(sprintf('\b'), 1, length(msg)+1); 
            end
            fprintf('\n');
        end
        disp('Finished!')
        
        %return results
        if strcmp(return_procedure,'plot') == 1
            %disp('Preparing for plotting')
            %[plot_val_array_x, plot_val_array_y, plot_val_array_z] = ...
            %    create_plot_array(val_array);
            %plot_surf(plot_val_array_x, plot_val_array_y, plot_val_array_z);
            try
                plot_scatter(val_array);
            catch exception
                warning(getReport(exception));
            end
        elseif strcmp(return_procedure,'ret_array') == 1
            varargout{1} = time_array;
            varargout{2} = val_array;
            varargout{3} = [result_type ' (' eval_unit ')'];
            if isempty(postprocess)==1
                return 
            end
        else
            warning('No valid return-procedure found, the default procedure will be applied (ret_array)')
            varargout{1} = time_array;
            varargout{2} = val_array;
            varargout{3} = [result_type ' (' eval_unit ')'];
            return
        end
end


%% postprocess the results
disp('------------------------------------------------------------')
disp('Postprocessing results')
disp(['Evaluating ' postprocess])
rho0 = 1.2;

postprocess_index_array = strcmp(postprocess,eval_list);
postprocess_index = find(postprocess_index_array);

[post_time_array, post_val_array, val_desc, post_eval_dim, cont_out] = ...
    feval(eval_file_list{postprocess_index},...
    result_type,postprocess,time_array,val_array,rho0);
if isempty(post_eval_dim)==1
    post_eval_dim = eval_dim;
end

%plotting needs to be done in the individual eval file!!!!
if strcmp(return_procedure_post,'plot')
    if post_eval_dim == 1
        disp('Plotting results')
        figure
        plot(post_time_array(:,1),post_val_array(:,1))
        titlestring = ['CFS-simulation ' filename];
        title(titlestring, 'Interpreter', 'none')
        xlabel('Time (s)')
        ylabel(val_desc)
        set(gca,'fontsize',20)
    elseif post_eval_dim == 2
        
    end
elseif strcmp(return_procedure_post,'ret_array')
    varargout{4} = post_time_array;
    varargout{5} = post_val_array;
    varargout{6} = val_desc;
    varargout{7} = cont_out;
else
    warning('No valid return-procedure found, the default procedure will be applied (ret_array)')
    varargout{4} = post_time_array;
    varargout{5} = post_val_array;
    varargout{6} = val_desc;
    varargout{7} = cont_out;
end

disp('Finished!')
disp('------------------------------------------------------------')
end