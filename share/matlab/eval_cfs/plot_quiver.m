function plot_quiver(x_plot,y_plot,val_array)
%PLOT_QUIVER (V1.0) 2D-quiver-plot 
%This function animates a quiver plot to be able to visualize a
%vector-field over time. This function can only handly a uniformly spaced
%grid, although each axis can have its own uniform spacing.

%initiate variables
time_steps = length(val_array);
max_vec = zeros(time_steps,2);
val_array_mat_scaled = cell(time_steps,2);
%extract unique coordinates for our grid
x_q = unique(x_plot);
y_q = unique(y_plot);
%extract the maximum value to be able to scale the plot
for mm=1:length(val_array)
    max_vec(mm,1) = max(max(abs(val_array{mm,1}{1,1})));
    max_vec(mm,2) = max(max(abs(val_array{mm,1}{1,2})));
end
max_val = max(max(max_vec));
[n,m] = size(val_array{1,1}{1,1});
%scale to element-size
scale_x = (max(x_plot(:))-min(x_plot(:)))/(n-1);
scale_y = (max(y_plot(:))-min(y_plot(:)))/(m-1);
for nn=1:length(val_array)
    val_array_mat_scaled{nn,1}=val_array{nn,1}{1,1}/max_val*scale_x;
    val_array_mat_scaled{nn,2}=val_array{nn,1}{1,2}/max_val*scale_y;
end
%animate the plot
figure
q = quiver(x_q,y_q,val_array_mat_scaled{1,1},val_array_mat_scaled{1,2},'AutoScale','off');
axis manual
for ii=2:time_steps
    set(q,'udata',val_array_mat_scaled{ii,1},'vdata',val_array_mat_scaled{ii,2})
    pause(0.05)
end
end