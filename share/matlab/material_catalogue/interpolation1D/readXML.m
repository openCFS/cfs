function [coeffs, params] = readXML(inputXML)

if exist(inputXML, 'file') ~= 2
    throw( MException('readXML:FileNotFound','File %s not found.',inputXML) );
end

doc = xmlread(inputXML);

paramnames = {
'param1';
'param2';
'param3'
};

coeffnames = {
'coeff11';
'coeff12';
'coeff13';
'coeff14';
'coeff15';
'coeff16';
'coeff22';
'coeff23';
'coeff24';
'coeff25';
'coeff26';
'coeff33';
'coeff34';
'coeff35';
'coeff36';
'coeff44';
'coeff45';
'coeff46';
'coeff55';
'coeff56';
'coeff66';
'vol';
'microloadfactor'
};

npoints = 0;
nintervals = 0;
ncoefficients = 0;

% read sizes of params
for i = 1:length(paramnames)
    list1 = doc.getElementsByTagName(paramnames{i});
    if list1.getLength < 1
        continue
    end
    list2 = list1.item(0).getElementsByTagName('matrix');
    npoints = str2double(list2.item(0).getAttribute('dim1'));
    break
end

% read size of coeff
for i = 1:length(coeffnames)
    list1 = doc.getElementsByTagName(coeffnames{i});
    if list1.getLength < 1
        continue
    end
    list2 = list1.item(0).getElementsByTagName('matrix');
    nintervals = str2double(list2.item(0).getAttribute('dim1'));
    ncoefficients = str2double(list2.item(0).getAttribute('dim2'));
    break
end

assert(npoints == nintervals+1);

% read params
ptemp = zeros(npoints, length(paramnames));
for i = 1:length(paramnames)
    list1 = doc.getElementsByTagName(paramnames{i});
    if list1.getLength > 1
        warning('MATLAB:readXML','File contains multiple %s', paramnames{i});
    end
    if list1.getLength < 1
%        warning('MATLAB:readXML','File contains no %s',name);
        continue
    end
    list2 = list1.item(0).getElementsByTagName('real');
    charArray = list2.item(0).getTextContent;
    ptemp(:,i) = str2num(charArray);
end

% read coeffs
ctemp = zeros(nintervals, ncoefficients, length(coeffnames));
for i = 1:length(coeffnames)
    name = coeffnames{i};
    list1 = doc.getElementsByTagName(name);
    if list1.getLength > 1
        warning('MATLAB:readXML','File contains multiple %s', name);
    end
    if list1.getLength < 1
%        warning('MATLAB:readXML','File contains no %s',name);
        continue
    end
    list2 = list1.item(0).getElementsByTagName('real');
    charArray = list2.item(0).getTextContent;
    ctemp(:,:,i) = str2num(charArray);
end

params = struct(...
'param1', ptemp(:,1),...
'param2', ptemp(:,1),...
'param3', ptemp(:,1)...
);

coeffs = struct(...
'coeff11', ctemp(:,:,1),...
'coeff12', ctemp(:,:,2),...
'coeff13', ctemp(:,:,3),...
'coeff14', ctemp(:,:,4),...
'coeff15', ctemp(:,:,5),...
'coeff16', ctemp(:,:,6),...
'coeff22', ctemp(:,:,7),...
'coeff23', ctemp(:,:,8),...
'coeff24', ctemp(:,:,9),...
'coeff25', ctemp(:,:,10),...
'coeff26', ctemp(:,:,11),...
'coeff33', ctemp(:,:,12),...
'coeff34', ctemp(:,:,13),...
'coeff35', ctemp(:,:,14),...
'coeff36', ctemp(:,:,15),...
'coeff44', ctemp(:,:,16),...
'coeff45', ctemp(:,:,17),...
'coeff46', ctemp(:,:,18),...
'coeff55', ctemp(:,:,19),...
'coeff56', ctemp(:,:,20),...
'coeff66', ctemp(:,:,21),...
'coeffvol', ctemp(:,:,22),...
'microloadfactor', ctemp(:,:,23)...
);
