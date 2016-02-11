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
scatter3(realdata(:,1),realdata(:,2),realdata(:,tensoridx),5,[0 .8 1],'*');
scatter3(realdata(:,1),realdata(:,2),intdata(:,tensoridx),10,'k');
scatter3(interpolationdatapoints(:,1),interpolationdatapoints(:,2),interpolationdatapoints(:,tensoridx),'r*');
hold off;
view(-30,20)
legh1 = legend(gca,'realdata', 'intdata','supp points','Location','SouthOutside','Orientation','horizontal');
pos1 = get(legh1,'Position');
set(legh1,'Position',[pos1(1) 0.02 pos1(3) pos1(4)]);
grid on

% plot relative error
subplot(1,2,2)
title('relative error')
hold on;
scatter3(realdata(:,1),realdata(:,2),relerr(:,tensoridx),10,relerr(:,tensoridx))
colormap('Cool')
scatter(interpolationdatapoints(:,1),interpolationdatapoints(:,2),'r*');
hold off;
view(-30,20)
legh2 = legend(gca,'reldiff','supp points','Location','SouthOutside','Orientation','horizontal');
pos2 = get(legh2,'Position');
set(legh2,'Position',[pos2(1) 0.02 pos2(3) pos2(4)]);
grid on

axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.29,0,...
    sprintf('maximal absolute error: %f   norm of absolute error: %f', aerrmax, aerrnorm),...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(0.75,0,...
    sprintf('maximal relative error: %f   norm of relative error: %f', relerrmax, relerrnorm),...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')

h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition', [1 1 28 19]);
print(h, '-dpdf', sprintf('relerr_3D_L%d_TI%d_p2_0%d.pdf', level, tensoridx, intdata(1,3)*1e6));
