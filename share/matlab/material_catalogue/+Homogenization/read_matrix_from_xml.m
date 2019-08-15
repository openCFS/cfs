function [ matrix ] = read_matrix_from_xml(infoxml)
% READ_MATRIX_FROM_XML reads the elasticity tensor from a xml file
%
%   matrix = read_matrix_from_xml(infoxml) gets the element
%   '//homogenizedTensor/ ... /real' from infoxml and stores
%   its text content into matrix
%
%   Example:
%   Eh = read_matrix_from_xml('./myTestrun.info.xml')
%

if exist(infoxml, 'file') ~= 2
    throw( MException('read_matrix_from_xml:FileNotFound','File %s not found.',infoxml) );
end

doc = xmlread(infoxml);

% Get Element '//homogenizedTensor/ ... /real'
list1 = doc.getElementsByTagName('homogenizedTensor');
if list1.getLength > 1
    warning('MATLAB:ReadMatrixFromXML','File contains multiple homogenized tensors');
end
if list1.getLength < 1
    warning('MATLAB:ReadMatrixFromXML','File contains no homogenized tensor');
end
list2 = list1.item(0).getElementsByTagName('real');
charArray = list2.item(0).getTextContent;
matrix = str2num(charArray);

end
