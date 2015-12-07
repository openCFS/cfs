tensoridx = 4;
dimstocheck = [2 3];
cataloguefile = 'catalogues/detailed_stats_presets3D_127';

threshold = 0.1;

% Get data from file
if ~exist('data','var')
    fid = fopen(cataloguefile);
    tmp = textscan(fid,'%dD\t%d\t%d');
    if isempty(tmp{3})
        fseek(fid,0,-1);
        tmp = textscan(fid,'%dD\tL%d\t%d');
    end
    dim = tmp{1};
    nPoints = tmp{3};
    fclose(fid);
    data = dlmread(cataloguefile,'\t',1,0);
else
    dim = size(data,2)-7;
end

idcs = setdiff(1:dim,dimstocheck);

Rows = unique(data(:,idcs),'rows');

index = 1:size(Rows,1);
nPointsperdim = index(end)+1;

wrongPoints = zeros(nPoints,3);
pidx = 1;

% Edge detection filter
filterkernel1 = -ones(3); filterkernel1(2,2) = 8;

for i = index
    p = [];
    % Take hyperplane
    [ism,idx] = ismember(data(:,idcs),Rows(i,:),'rows');
    B = data(logical(idx),:);
    % Put values into 2-D matrix
    B(:,1:3) = int32( B(:,1:3)*nPointsperdim );
    C = zeros(index(end), index(end));
    C( sub2ind(size(C), B(:,dimstocheck(1)), B(:,dimstocheck(2)) ) ) = B(:,tensoridx+3);
    % Filter for edge (peak) detection 
    D = abs(filter2(filterkernel1,C));
    
    % Take values above threshold
    E = D > threshold;
    % Sum up 5x5 parts of matrix using a simple sum filter
    filterkernel2 = ones(5);
    F = abs(filter2(filterkernel2,E));
    % Get indices of peaks (sum == 25)
    [idxi,idxj] = find(abs(F - 25) < 1e-10);
    p = [p; idxi+1, idxj+1 , i*ones(length(idxi),1)];
    for ii = -2:2
        for jj = -2:2
            E(idxi+ii,idxj+jj) = 0;
        end
    end
    % Get indices of peaks (sum == 9)
    filterkernel2 = ones(3);
    F = abs(filter2(filterkernel2,E));
    [idxi,idxj] = find(abs(F - 9) < 1e-10);
    p = [p; idxi, idxj , i*ones(length(idxi),1)];
    for ii = -1:1
        for jj = -1:1
            E(idxi+ii,idxj+jj) = 0;
        end
    end
    % Get indices of peaks (sum == 1)
    [idxi,idxj] = find( E(2:end-1,2:end-1) );
    p = [p; idxi+1, idxj+1 , i*ones(length(idxi),1)];
    wrongPoints(pidx:pidx+size(p,1)-1,:) = p(:,[dimstocheck idcs]);
    pidx = pidx + size(p,1);
    wrongPoints = wrongPoints(wrongPoints(:,1) > 0,:);

%     subplot(2,2,1)
%     surf(C)
%     zlim([0 1.2])
%     subplot(2,2,2)
%     surf(D(2:end-1,2:end-1))
%     zlim([0 1.2])
%     subplot(2,2,3)
%     surf(+E)
%     subplot(2,2,4)
%     surf(F)
end
