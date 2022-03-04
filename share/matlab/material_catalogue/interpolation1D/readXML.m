function coeffs = readXML(inputXML)

if exist(inputXML, 'file') ~= 2
    throw( MException('readXML:FileNotFound','File %s not found.',inputXML) );
end

doc = xmlread(inputXML);

names = {
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


% read sizes
for i = 1:length(names)
    name = names{i};
    list1 = doc.getElementsByTagName(name);
    if list1.getLength < 1
        continue
    end
    list2 = list1.item(0).getElementsByTagName('matrix');
    nintervals = str2double(list2.item(0).getAttribute('dim1'));
    ncoefficients = str2double(list2.item(0).getAttribute('dim2'));
    break
end

% read coeffs
temp = zeros(nintervals, ncoefficients, length(names));
for i = 1:length(names)
    name = names{i};
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
    temp(:,:,i) = str2num(charArray);
end

coeffs = struct(...
'coeff11', temp(:,:,1),...
'coeff12', temp(:,:,2),...
'coeff13', temp(:,:,3),...
'coeff14', temp(:,:,4),...
'coeff15', temp(:,:,5),...
'coeff16', temp(:,:,6),...
'coeff22', temp(:,:,7),...
'coeff23', temp(:,:,8),...
'coeff24', temp(:,:,9),...
'coeff25', temp(:,:,10),...
'coeff26', temp(:,:,11),...
'coeff33', temp(:,:,12),...
'coeff34', temp(:,:,13),...
'coeff35', temp(:,:,14),...
'coeff36', temp(:,:,15),...
'coeff44', temp(:,:,16),...
'coeff45', temp(:,:,17),...
'coeff46', temp(:,:,18),...
'coeff55', temp(:,:,19),...
'coeff56', temp(:,:,20),...
'coeff66', temp(:,:,21),...
'coeffvol', temp(:,:,22),...
'microloadfactor', temp(:,:,23)...
);

