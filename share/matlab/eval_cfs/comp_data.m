function [comp_x, comp_y, comp_val_array_mat, plane_info] = comp_data(x_vec,y_vec,z_vec,val_array_mat)
%COMP_DATA (V1.0) checks for empty entries and returns the compressed values-array with the
%corresponding vectors
comp_val_array_mat = cell(length(val_array_mat),1);
if isempty(val_array_mat{1,1}{1,1})==1
    %yz-plane
    comp_x = y_vec;
    comp_y = z_vec;
    for ii=1:length(val_array_mat)
        sub_val_array = cell(1,2);
        sub_val_array{1,1} = val_array_mat{ii,1}{1,2};
        sub_val_array{1,2} = val_array_mat{ii,1}{1,3};
        comp_val_array_mat{ii} = sub_val_array; 
    end
    plane_info = 'yz';
elseif isempty(val_array_mat{1,1}{1,2})==1
    %zx-plane
    comp_x = x_vec;
    comp_y = z_vec;
    for ii=1:length(val_array_mat)
        sub_val_array = cell(1,2);
        sub_val_array{1,1} = val_array_mat{ii,1}{1,1};
        sub_val_array{1,2} = val_array_mat{ii,1}{1,3};
        comp_val_array_mat{ii} = sub_val_array; 
    end
    plane_info = 'zx';
elseif isempty(val_array_mat{1,1}{1,3})==1
    %xy-plane
    comp_x = x_vec;
    comp_y = y_vec;
    for ii=1:length(val_array_mat)
        sub_val_array = cell(1,2);
        sub_val_array{1,1} = val_array_mat{ii,1}{1,1};
        sub_val_array{1,2} = val_array_mat{ii,1}{1,2};
        comp_val_array_mat{ii} = sub_val_array; 
    end
    plane_info = 'xy';
end