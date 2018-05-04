function array_out = read_val_array(ii,f_path,eval_coord_r_index,...
    eval_region,result_type,node_index,eval_type,eval_coord_r)
time_string = ['/Results/Mesh/MultiStep_1/Step_' num2str(ii)];
time_array = h5readatt(f_path,time_string,'StepValue');
last_region = '';
eval_region_short=cellstr(eval_region(~cellfun('isempty',eval_region)));
%test2=cell2table(test)
%unk=unique(test2)
%unk{:,{'test'}}
[n_nva, ~] = size(h5read(f_path,...
    ['/Results/Mesh/MultiStep_1/Step_1/' num2str(result_type)...
    '/' eval_region_short{1} '/Nodes/Real']));
sub_val_array = NaN(length(eval_coord_r_index),3+n_nva);
for oo=1:length(eval_coord_r_index)
    %create an array with x- and y-coordinates including the
    %value of the respective point
    
    %check if a result region is available for the nodes
    if isempty(eval_region{oo}) == 1
        sub_val_array(oo,4:end) = NaN;
    else
        %load the required region
        if strcmp(last_region,eval_region{oo}) == 0
            val_string = ['/Results/Mesh/MultiStep_1/Step_' num2str(ii) '/' ...
                num2str(result_type) '/' eval_region{oo} '/Nodes/Real'];
            node_val_array = h5read(f_path,val_string);
            last_region = eval_region{oo};
        end
        %read the value from the file
        sub_val_array(oo,4:end) = node_val_array(:,node_index(oo));
    end

    %we need to round here because otherwise slight 
    %inacuracies destroy the uniqueness of certain values
    switch eval_type
        case 'xy' 
            %add the x- and y-coordinates
            sub_val_array(oo,1) = eval_coord_r(oo,1);
            sub_val_array(oo,2) = eval_coord_r(oo,2);
        case 'yz'
            %add the y- and z-coordinates
            sub_val_array(oo,2) = eval_coord_r(oo,2);
            sub_val_array(oo,3) = eval_coord_r(oo,3);
        case 'zx'
            %add the z- and x-coordinates
            sub_val_array(oo,3) = eval_coord_r(oo,3);
            sub_val_array(oo,1) = eval_coord_r(oo,1);
    end
    %write it to the whole array
    val_array = sub_val_array;
end
array_out = cell(1,2);
array_out{1} = time_array;
array_out{2} = val_array;
end