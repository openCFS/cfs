function plot_point(varargin)
%plot the data of a single point over time
time_array_vec = varargin{1};
val_array = varargin{2};
[~, m_v] = size(varargin);
switch m_v
    case 3
        coord = 4;
    case 4
        coord_string = varargin{4};
        switch coord_string
            case 'x'
                coord = 4;
            case 'y'
                coord = 5;
            case 'z'
                coord = 6;
        end
end
val_array_vec = zeros(length(val_array),1);
for ii=1:length(val_array)
    val_array_vec(ii,1) = val_array{ii}(1,coord);
end
figure
plot(time_array_vec,val_array_vec)
xlabel('Time (s)')
ylabel(varargin{3})
end