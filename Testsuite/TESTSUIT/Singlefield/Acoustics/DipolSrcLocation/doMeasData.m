
fileNames = {'history/forwardSrc500-acouPressure-node-172-MicPos.hist', ...
             'history/forwardSrc500-acouPressure-node-204-MicPos.hist', ...
             'history/forwardSrc500-acouPressure-node-216-MicPos.hist', ...
             'history/forwardSrc500-acouPressure-node-277-MicPos.hist', ...
             'history/forwardSrc500-acouPressure-node-411-MicPos.hist', ...
             'history/forwardSrc500-acouPressure-node-446-MicPos.hist'};

outputName = 'testMeas.txt';
         
%do data steering         
[dummy numFiles] = size(fileNames);
measData = zeros(numFiles,3);
%nodeNumbers = zeros(numFiles,1);

for i=1:numFiles
    [data, nodeNr] = getDataFromSingleHist(fileNames(i));
    measData(i,:)  = [nodeNr; data(2); data(3) ];   
    %%nodeNumbers(i) = nodeNr;
end

%save to outout file
fileID = fopen(outputName,'w');
for i=1:numFiles
    fprintf(fileID,'%6.0f %12.8f %12.8f\n',measData(i,:));
end
fclose(fileID);

