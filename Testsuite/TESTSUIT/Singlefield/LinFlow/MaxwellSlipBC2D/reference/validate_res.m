%% Evaluate Maxwell slip BC

eval_step = 15; % last step
hy = 0.04; % height of evaluation node
w = 1e-3; % width of the channel, evaluation at right boundary of chamber


%% Sequence Steps (1-hollowincompressible, 2-incompressible, 3-compressible)

for seqStep = 1:3  
    
% read in 
[tx,res_cellx,coords_cellx] = read_res('results_hdf5/MaxwellSlipBC2D_fine.cfs','y',hy,'fluidMechVelocity',seqStep);
[ty,res_celly,coords_celly] = read_res('results_hdf5/MaxwellSlipBC2D_fine.cfs','x',w,'fluidMechVelocity',seqStep);
% sort

x_coords = zeros(length(coords_cellx),1);
y_coords = zeros(length(coords_celly),1);
sorted_res_cellx = cell(size(res_cellx));
sorted_res_celly = cell(size(res_celly));

for ii=1:length(coords_cellx)
    x_coords(ii) = coords_cellx{ii}(1);
end

for ii=1:length(coords_celly)
    y_coords(ii) = coords_celly{ii}(2);
end

[sorted_x_coords, sorting_vecx] = sort(x_coords);
[sorted_y_coords, sorting_vecy] = sort(y_coords);

for ii=1:length(coords_cellx)
    sorted_res_cellx{ii} = res_cellx{sorting_vecx(ii)};
end

for ii=1:length(coords_celly)
    sorted_res_celly{ii} = res_celly{sorting_vecy(ii)};
    if sorted_y_coords(ii) == hy
        eval_node_index_y = ii;
    end        
end

% evaluate

resx = zeros(length(coords_cellx),1);
resy = zeros(length(coords_celly),1);

for ii=1:length(coords_cellx)
    dummy_vec = sorted_res_cellx{ii};
    resx(ii) = dummy_vec(eval_step);
end
for ii=1:length(coords_celly)
    dummy_vec = sorted_res_celly{ii};
    resy(ii) = dummy_vec(eval_step);
end

%figure
%plot(sorted_x_coords,resx)

sorted_x_coords_diff = sorted_x_coords(1:end-1)+diff(sorted_x_coords);
duydx = diff(resx)./diff(sorted_x_coords);

duxdy = diff(resy)./diff(sorted_y_coords);

% backwards difference at eval point
duxdy(eval_node_index_y);
% forward difference at eval point
duxdy(eval_node_index_y+1);

% store results for seqSteps
if seqStep==1
    duydx_1=duydx;
    duxdy_1=duxdy;
    resx_1=resx; % results at h=0.04
    resy_1=resy; % results at x=right boundary of channel
elseif seqStep==2
    duydx_2=duydx;
    duxdy_2=duxdy;
    resx_2=resx;
    resy_2=resy;
elseif seqStep==3
    duydx_3=duydx;
    duxdy_3=duxdy;
    resx_3=resx;
    resy_3=resy;
end
end
%figure
%plot(sorted_x_coords_diff,duydx)

C=1000;
meanFreePath=68e-9;
% the Maxwell BC prescribes -(duy/dx+dux/dy)*C*meanFreePath = u_t
%disp(['Tangential velocity based on the normal derivative (nominal value): ' num2str(-duydx(end)*meanFreePath*C) ' m/s'])
%disp(['Actual tangential velocity (actual value): ' num2str(resx(end)+resy(eval_node_index_y+1)) ' m/s'])
%disp(['Actual tangential velocity (actual value): ' num2str(resx(end)) ' m/s'])

disp(['Tangential velocity based on the normal derivative for seqStep 1 (nominal value): ' num2str(-(duydx_1(end)+duxdy_1(eval_node_index_y))*meanFreePath*C) ' m/s'])
disp(['Actual tangential velocity for seqStep 1 (actual value): ' num2str(resx_1(end)) ' m/s'])
disp(['Tangential velocity based on the normal derivative for seqStep 2 (nominal value): ' num2str(-(duydx_2(end)+duxdy_2(eval_node_index_y))*meanFreePath*C) ' m/s'])
disp(['Actual tangential velocity for seqStep 2 (actual value): ' num2str(resx_2(end)) ' m/s'])
disp(['Tangential velocity based on the normal derivative for seqStep 3 (nominal value): ' num2str(-(duydx_3(end)+duxdy_3(eval_node_index_y))*meanFreePath*C) ' m/s'])
disp(['Actual tangential velocity for seqStep 3 (actual value): ' num2str(resx_3(end)) ' m/s'])
