%% eval_h5 example usage
%This example-file demonstrates how the eval_h5 package can be used in
%four different cases.
%For the examples we consider a coupled simulation
%(acoustic-electrostatic-mechanic) where the acoustic problem has been 
%solved via the 'acouPotential' formulation

%define the path of the file
%fpath = 'C:\TEST\EXAMPLE.cfs';
fpath = 'C:\Users\Dominik\Documents\MEMS\CFS\MEMS_couple\MEMS_Couple_alpha_direct.cfs';
%define an arbitrary point
p = [0 0 0];

%% Example 1: Extract the values of an arbitrary point
%we do not need the source data, therefore we just extract the already
%postprocessed data
[~, ~, ~, ex1_post_time_array, ex1_post_val_array, ex1_post_unit,~] = ...
    eval_cfs(fpath,'acouPotential',1,'point',p,'acouPressure','ret_array');
plot_point(ex1_post_time_array, ex1_post_val_array, ex1_post_unit);

%% Example 2: Extract the values of an arbitrary point
[ex2_time_array, ex2_val_array, ex2_unit] = eval_cfs(fpath,...
    'mechDisplacement',1,'point',p,'','ret_array');
plot_point(ex2_time_array,ex2_val_array,ex2_unit,'z')

%% Example 3: Extract the values of a plane and animate it
%we do not need the source data, therefore we just extract the already
%postprocessed data
[~, ~, ~, ex3_post_time_array, ex3_post_val_array, ex3_post_unit, ...
    ex3_cont_out] = eval_cfs(fpath,'acouPotential',1,'yz',p,...
    'acouPressure','ret_array');
plot_scatter(ex3_post_val_array)

%% Example 4: Extract the values of plane and use a vector plot
%we do not need the source data, therefore we just extract the already
%postprocessed data
[~, ~, ~, ex4_post_time_array, ex4_post_val_array, ex4_post_unit, ...
    ex4_cont_out] = eval_cfs(fpath,'acouPotential',1,'yz',p,...
    'acouVelocity','ret_array');
[ex4_comp_x, ex4_comp_y, ex4_comp_val_array_mat, ex4_plane_info] = ...
    comp_data(ex4_cont_out{2},ex4_cont_out{3},ex4_cont_out{4},...
    ex4_post_val_array);
plot_quiver(ex4_comp_x,ex4_comp_y,ex4_comp_val_array_mat)
