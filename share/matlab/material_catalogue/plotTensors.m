tensoridx = 1;
dimstoplot = [1 2];
plotwithinterpolatedvalues = 0;
onlyboundary = 0;
% cataloguefile = 'catalogues/detailed_stats_presets3D_65_bndry';
cataloguefile = 'catalogues/detailed_stats_presets2D_11';
% cataloguefile2 = 'catalogues/detailed_stats_presets3D_L9';
% intfile = 'interpolatedValues_128';

Rows = [
% 16/32
% 0.0
% 0.125
% 0.25
% 0.5
% 0.75
% 0.875
% 1.0
];

Rows = linspace(0.01,0.99,9)';
Rows = .5;

% Rows = [
% 1/np 1/np
% 2/128 2/128
% 0.0	0.0
% 0.0	0.5
% 0.0	1.0
% 0.5	0.0
% 0.5	0.5
% 0.5	1.0
% 1.0	0.0
% 1.0	0.5
% 1.0	1.0
% 0 0
% 0.05 0.05
% 0.1 0.1
% 0.25 0.25
% 0.25 0.5
% 0.25 0.75
% 0.5 0.25
% 0.5 0.5
% 0.5 0.75
% 0.75 0.25
% 0.75 0.5
% 0.75 0.75
% ];

% Rows = [
% 0.0	0.0	0.0
% 0.0	0.0	0.5
% 0.0	0.0	1.0
% 0.0	0.5	0.0
% 0.0	0.5	0.5
% 0.0	0.5	1.0
% 0.0	1.0	0.0
% 0.0	1.0	0.5
% 0.0	1.0	1.0
% 0.5	0.0	0.0
% 0.5	0.0	0.5
% 0.5	0.0	1.0
% 0.5	0.5	0.0
% 0.5	0.5	0.5
% 0.5	0.5	1.0
% 0.5	1.0	0.0
% 0.5	1.0	0.5
% 0.5	1.0	1.0
% 1.0	0.0	0.0
% 1.0	0.0	0.5
% 1.0	0.0	1.0
% 1.0	0.5	0.0
% 1.0	0.5	0.5
% 1.0	0.5	1.0
% 1.0	1.0	0.0
% 1.0	1.0	0.5
% 1.0	1.0	1.0
% ];


% Get data from file
if ~exist('data','var')
    fid = fopen(cataloguefile);
    tmp = textscan(fid,'%dD\t%d\t%d');
    if isempty(tmp{3})
        fseek(fid,0,-1);
        tmp = textscan(fid,'%dD\tL%d\t%d');
    end
    dim = double(tmp{1});
    nPoints = double(tmp{3});
    fclose(fid);
    data = dlmread(cataloguefile,'\t',1,0);
else
    dim = size(data,2)-7;
end
if plotwithinterpolatedvalues
    if ~exist('intVals','var')
        intVals = load(intfile);
    end
    s = strsplit(intfile,'_');
    np = str2double(s{2});
end

idcs = setdiff(1:dim,dimstoplot);

% Rows = unique(data(:,idcs),'rows');

% index = load('testing.txt')';
index = 1:size(Rows,1);

% Get full material tensor
[~,idx] = ismember(data(:,1:dim),ones(1,dim),'rows');
fullmaterial = data(idx|zeros(size(data,1),1),:);

if exist('cataloguefile2','var')
    data2 = dlmread(cataloguefile2,'\t',1,0);
end

figure
for i = index
    % Take hyperplane
    [ism,idx] = ismembertol(data(:,idcs),Rows(i,:),1e-10,'ByRows',true);
    B = data(logical(idx),:);
    if exist('cataloguefile2','var')
        data3 = data2(data2(:,idcs)==Rows(i,:),:);
    end
    
%     scatter(B(:,dim),B(:,tensoridx+4));

    % Get only values at boundary
    if onlyboundary
        fac = zeros(size(B,1),1);
        prod = ones(size(B,1),1);
        for j = 1:length(dimstoplot)
            fac = B(:,dimstoplot(j)) - ones(size(B,1),1);
            prod = prod .* B(:,dimstoplot(j)) .* fac;
        end
        boundaryidcs = abs(prod) < eps;
        B = B(boundaryidcs,:);
    end
    
    assert(~isempty(B),sprintf('No points to plot. Choose other row values.'))
        
    
    % Remove full material
%     B = B( abs(B(:,tensoridx+4) - fullmaterial(tensoridx+4)) > 1e-10,:);

    if plotwithinterpolatedvalues
        [ismIV,idxIV] = ismember(intVals(:,idcs),Rows(i,:),'rows');
        if ~( sum(ismIV))
            continue
        end
        C = intVals(logical(idxIV),:);
        nPointsPerDim = sqrt(size(C,1))-1;
    else
        nPointsPerDim = size(data,1)^(1/dim)+1;
    end
    
    % Plot data
    if length(dimstoplot) == 2
        if plotwithinterpolatedvalues
            assert(abs(mod(C(1,dimstoplot(1))*nPointsPerDim,1)) < 1e-3)
            S = sparse(round(C(:,dimstoplot(1))*nPointsPerDim+1),round(C(:,dimstoplot(2))*nPointsPerDim+1),C(:,tensoridx+dim));
            [X,Y] = meshgrid(unique(C(:,dimstoplot(1))), unique(C(:,dimstoplot(2))));
            surf(X,Y,S')
            shading interp
            hold on
            if exist('cataloguefile2','var')
                scatter3(data3(:,dimstoplot(1)),data3(:,dimstoplot(2)),data3(:,tensoridx+dim),'r')
            end
        else
            if mod(nPointsPerDim,1) == 0
                assert(abs(mod(B(1,dimstoplot(1))*nPointsPerDim,1)-1) < 1e-3,sprintf('%f',B(1,dimstoplot(1))*nPointsPerDim))
                S = sparse(round(B(:,dimstoplot(1))*nPointsPerDim),round(B(:,dimstoplot(2))*nPointsPerDim),B(:,tensoridx+dim));
                [X,Y] = meshgrid(unique(B(:,dimstoplot(1))), unique(B(:,dimstoplot(2))));
                surf(X,Y,S')
            end
        end
        colorbar
%         if plotwithinterpolatedvalues
%             hold on;
%             scatter3(B(:,dimstoplot(1)),B(:,dimstoplot(2)),B(:,tensoridx+4),25,B(:,tensoridx+dim))
            scatter3(B(:,dimstoplot(1)),B(:,dimstoplot(2)),B(:,tensoridx+dim),10,'k*')
%             hold off;
%         end
        xlim([-0.01 1])
        ylim([-0.01 1])
        zlim([-0.01 1.2])
        xlabel( sprintf('Dimension %d', dimstoplot(1)) );
        ylabel( sprintf('Dimension %d', dimstoplot(2)) );
        zlabel( sprintf('Tensor entry %d', tensoridx) );
        title( sprintf('Other dimension(s) %s', sprintf('%.8f ', Rows(i,:))) );
        grid minor;
        hold off;
        input('Press enter for next hyperplane');
        continue;
    elseif length(dimstoplot) == 3
        scatter3(B(:,dimstoplot(1)),B(:,dimstoplot(2)),B(:,dimstoplot(3)),25*B(:,tensoridx+4)+1,B(:,tensoridx+4))
        xlim([-0.01 1])
        ylim([-0.01 1])
        zlim([-0.01 1])
        xlabel( sprintf('Dimension %d', dimstoplot(1)) );
        ylabel( sprintf('Dimension %d', dimstoplot(2)) );
        zlabel( sprintf('Dimension %d', dimstoplot(3)) );
        title( sprintf('Other dimension(s) %s', sprintf('%.8f ', Rows(i,:))) );
        grid minor;
        hold off;
        waitforbuttonpress;
        continue;
    end

    reply = input('Print Index to file? Y/N [N]:','s');
    if isempty(reply)
        continue;
    end
    disp(i)

    xindex = input('xIndex:','s');
    C = strsplit(xindex, {' ', '\t', ',', ', '});
    xidx = zeros(1,length(C));
    for j = 1:length(C)
     xidx(j) = str2double(C{j});
    end
%     yindex = input('yIndex:','s');
%     C = strsplit(yindex, {' ', '\t', ',', ', '});
%     yidx = zeros(1,length(C));
%     for j = 1:length(C)
%         yidx(j) = str2double(C{j});
%     end
    level = input('level:');
    
    [~,idx] = ismember(B(:,dimstoplot),xidx' / 2^level,'rows');
    C = B(idx|zeros(size(B,1),1),:);

%     % figure;
%     scatter(C(:,1),C(:,tensoridx+4),25,C(:,tensoridx+4))

    fid = fopen('testing.txt','at');
    fprintf(fid,'%d\n',i);
    fclose(fid);
    
    fid = fopen('wrongpoints.txt','at');
%     %fprintf(fid,'%.10f\t%.10f\t%.10f\t%.10f\n', C(C(:,tensoridx+4)>eps,1:4)');
    fprintf(fid,'%.10f\t%.10f\t%.10f\t%.10f\n', C(:,1:4)');
    fclose(fid);
end
