
function [data, nodeNr] = getDataFromSingleHist(filenameStr)

%filename = 'history/forwardSrc500-ms1-acouPressure-node-1144-MicPos.hist';
filename = char(filenameStr);

mmfile = fopen(filename,'r');
if ( mmfile == -1 )
    disp(filename);
    error('File not found');
end;

%get node number
commentline = fgets(mmfile);
lineStr =  sscanf(commentline,'%c');
pos=strfind(lineStr,'#');
nodeStr = lineStr(pos(2)+1:length(lineStr));
nodeNr = round(str2num(nodeStr));

while length(commentline) > 0 & commentline(1) == '#',
    commentline = fgets(mmfile);
end
data = sscanf(commentline,'%lf%lf%lf');

%close file
fclose(mmfile);

end