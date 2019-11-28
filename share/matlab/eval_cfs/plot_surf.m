function plot_surf(varargin)
%PLOT_SURF (V1.0) this function creates a surface plot and animates it
switch nargin
    case 3
        x_coord = varargin{1};
        y_coord = varargin{2};
        val_mat = varargin{3};
        view_point = '3D';
    case 4
        x_coord = varargin{1};
        y_coord = varargin{2};
        val_mat = varargin{3};
        view_point = varargin{4};
end

x_vec = unique(x_coord);
y_vec = unique(y_coord);
max_vec = zeros(length(val_mat),1);
for mm=1:length(val_mat)
    max_vec(mm,1) = max(max(abs(val_mat{mm})));
end
max_val = max(max_vec);
figure
s = surf(x_vec,y_vec,val_mat{1});
axis([xlim ylim -max_val max_val])
axis manual
caxis([-max_val max_val])
caxis manual

switch view_point
    case '3D'
        %standard mode
    case 'top'
        view(2)
end

for ii=2:length(val_mat)
    s.ZData = val_mat{ii};
    pause(0.05)
end
end