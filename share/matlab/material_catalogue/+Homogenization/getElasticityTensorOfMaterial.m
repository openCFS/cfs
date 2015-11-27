function [ C ] = getElasticityTensorOfMaterial(file)
% GETELASTICITYTENSOR returns the (isotropic) elasticity tensor for the
% material given in file.
%   
%   Example:
%   C = getElasticityTensor('./inv_tensor.xml')
%

doc = xmlread(file);

% Get Element '//domain/ ... /region'
list1 = doc.getElementsByTagName('domain');
list2 = list1.item(0).getElementsByTagName('region');

% Get material type
attributes = list2.item(0).getAttributes;
for i = 0:attributes.getLength-1
    if attributes.item(i).getName.equals('material')
        material = char(attributes.item(i).getValue);
        break;
    end
end

switch material
    case 'Steel'
        E = 2.00E+11;
        nu = 0.29;
    case '99lines'
        E = 1;
        nu = 0.3;
    case '99mod'
        E = 1.00E+12;
        nu = 0.3;
    case 'sigmund'
        E = 5;
        nu = 0.35;
    case 'sigmund_void'
        E = 0.00086;
        nu = 0.23076;
end

C = E / (1-nu^2) * [ 1 , nu , 0 ; nu , 1 , 0 ; 0 , 0 , (1-nu)/2 ];
