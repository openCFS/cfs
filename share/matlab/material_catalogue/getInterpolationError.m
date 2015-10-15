% function getInterpolationError(interpolatedcataloguefile, exactcataloguefile)

level = 5;
tensoridx = 5;

interpolationcataloguefile = sprintf('catalogues/detailed_stats_presets3D_L%d',level);
interpolatedcataloguefile = sprintf('catalogues/detailed_stats_presets3D_L%d_interpolant',level);
exactcataloguefile = sprintf('catalogues/detailed_stats_presets3D_%d',2^level-1);

%-------------------------------------------------------------------------%

% Get interpolated data from file
intdata = dlmread(interpolatedcataloguefile,'\t',0,0);

% Get real data from file
realdata = dlmread(exactcataloguefile,'\t',1,0);
realdata = realdata( realdata(:,3)==intdata(1,3), :);

% Get interpolation data from file
interpolationdata = dlmread(interpolationcataloguefile,'\t',1,0);
interpolationdatapoints = interpolationdata( interpolationdata(:,3)==intdata(1,3) ,:);

[nintpoints,m1] = size(intdata);
[nrealpoints,m2] = size(realdata);

if nintpoints ~= nrealpoints
    error('Catalogues do not match.');
end

m = min(m1,m2);

% absolute error
aerr = abs( intdata(:,1:m) - realdata(:,1:m) );
aerrnorm = norm(aerr(:,tensoridx));
aerrmax = max(abs(aerr(:,tensoridx)));

if any(aerr(:,1:3) > 1e-6)
    error('Catalogues do not match.');
end

% relative error
relerr = aerr./realdata(:,1:m);
relerrnorm = norm(relerr(:,tensoridx));
relerrmax = max(abs(relerr(:,tensoridx)));

fprintf('\n')
fprintf('level:       %d\n', level)
fprintf('tensorindex: %d\n', tensoridx)
fprintf('p[2] = %f\n', intdata(1,3))
fprintf('\n')
fprintf('maximal absolute error: %f\n', aerrmax)
fprintf('norm of absolute error: %f\n', aerrnorm)
fprintf('\n')
fprintf('maximal relative error: %f\n', relerrmax)
fprintf('norm of relative error: %f\n', relerrnorm)
fprintf('\n')

% plot data
subplot(1,2,1)
title( sprintf('Tensorindex %d for p[2]=%f', tensoridx, intdata(1,3)) )
hold on;
scatter3(realdata(:,1),realdata(:,2),realdata(:,tensoridx),10,[0 .8 1],'*');
scatter3(realdata(:,1),realdata(:,2),intdata(:,tensoridx),15,'k');
scatter3(interpolationdatapoints(:,1),interpolationdatapoints(:,2),interpolationdatapoints(:,tensoridx),'r*');
hold off;
legend('realdata', 'intdata','supp points','Location','SouthOutside')
grid on
view(-30,20)

% plot relative error
subplot(1,2,2)
title('relative error')
hold on;
scatter3(realdata(:,1),realdata(:,2),relerr(:,tensoridx),10,relerr(:,tensoridx))
zlim([0,10])
colormap('Cool')
scatter(interpolationdatapoints(:,1),interpolationdatapoints(:,2),'r*');
hold off;
legend('reldiff','supp points','Location','SouthOutside')
grid on
view(-30,20)