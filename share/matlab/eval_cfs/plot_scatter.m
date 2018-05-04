function plot_scatter(varargin)
%PLOT_SCATTER (V1.0) this function creates a scatter plot and animates it
switch nargin
    case 1
        comp_val_array = varargin{1};
        view_point = '3D';
    case 2
        comp_val_array = varargin{1};
        view_point = varargin{2};
end

time_steps = length(comp_val_array);
max_vec = zeros(time_steps,1);
for mm=1:length(comp_val_array)
    max_vec(mm,1) = max(abs(comp_val_array{mm}(:,4)));
end
max_val = max(max_vec);
figure
if sum(isnan(comp_val_array{1}(:,3)))>0
    %xy-plane
    s = scatter3(comp_val_array{1}(:,1),comp_val_array{1}(:,2),comp_val_array{1}(:,4),'filled','CData',comp_val_array{1}(:,4));
elseif sum(isnan(comp_val_array{1}(:,1)))>0
    %yz-plane
    s = scatter3(comp_val_array{1}(:,2),comp_val_array{1}(:,3),comp_val_array{1}(:,4),'filled','CData',comp_val_array{1}(:,4));
elseif sum(isnan(comp_val_array{1}(:,2)))>0
    %zx-plane
    s = scatter3(comp_val_array{1}(:,1),comp_val_array{1}(:,3),comp_val_array{1}(:,4),'filled','CData',comp_val_array{1}(:,4));
end

axis([xlim ylim -max_val max_val])
axis manual
caxis([-max_val max_val])
caxis manual
alpha(.5)

switch view_point
    case '3D'
        %standard mode
    case 'top'
        view(2)
end

for ii=2:time_steps
    s.ZData = comp_val_array{ii}(:,4);
    s.CData = comp_val_array{ii}(:,4);
    pause(0.05)
end
end