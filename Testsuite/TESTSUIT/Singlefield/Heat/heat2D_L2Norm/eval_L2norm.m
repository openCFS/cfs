%% L2 norm verification

% here we calculate the corresponding norms of the problem
% for the tests using the relative L2 norm we expect to get either inf (for
% the case when the reference is zero) or 1 (if the excitation is zero) for
% all time steps

n_elem = 5;
d = 0.1;
amp = 1;

%% Full/over integration

% analytical solution of \int_\Omega (amp*y/d)^2 d\Omega

sol = @(y) (amp^2*y^3/(3*d^2))*d;

squared_norm_full = sol(d);

norm_full = sqrt(squared_norm_full);

disp(['Norm obtained with correct integration order: ' num2str(norm_full)])


%% Under integration

y_node = 0:d/n_elem:d;

y_mid = y_node(1:end-1)+d/n_elem/2;

disp_func = @(y) amp*y/d;


% calc approximation

squared_norm = 0;
for ii=1:length(y_mid)
    squared_norm = squared_norm + disp_func(y_mid(ii))^2*d*d/n_elem;
end
norm = sqrt(squared_norm);

disp(['Norm obtained with too low integration order (order 1): ' num2str(norm)])

