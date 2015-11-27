function callCatalogueGeneration(gridfile)
% CALLCATALOGUEGENERATION is an interface for callMatlab.sh to the
% Homogenization package.

tmp = strsplit(gridfile,{'/','\'});
threadID = str2double(tmp{end});
if isnan(threadID)
    threadID = tmp{end};
end

Homogenization.getInterpolationValues(gridfile, threadID)
