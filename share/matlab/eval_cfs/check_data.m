function nan_vec = check_data(f_path,eval_coord_r_index,...
    eval_region,result_type,node_index,eval_type,eval_coord_r)
%Checks if there are results for the data points, otherwise they will be
%disregarded
array_out = read_val_array(1,f_path,eval_coord_r_index,eval_region,...
    result_type,node_index,eval_type,eval_coord_r);
val_vec = array_out{2}(:,4);
nan_vec = isnan(val_vec);
end
