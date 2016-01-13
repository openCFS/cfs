function callCatalogueGeneration(gridfile)
% CALLCATALOGUEGENERATION is an interface for callMatlab.sh to the
% Homogenization package.

if ischar(gridfile)
    tmp = strsplit(gridfile,{'/','\'});
    threadID = str2double(tmp{end});
    if isnan(threadID)
        threadID = tmp{end};
    end
else
    threadID = 1;
end

Homogenization.getInterpolationValues(gridfile, threadID)
