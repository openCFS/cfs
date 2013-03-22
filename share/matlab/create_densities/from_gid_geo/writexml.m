function writexml( filename, elements, density )
disp('writing xml file');

xml_file = fopen(filename, 'w');

fprintf(xml_file, '<?xml version="1.0"?>\n<cfsErsatzMaterial writeCounter="16">\n <header>\n  <design constant="false" initial="0.5" name="density" region="all" scale="false"/>\n  <transferFunction application="mech" design="density" param="1.0" type="simp"/>\n </header>\n<set id="struct">\n');


nelements = size(elements,1);

for i=1:nelements
   fprintf(xml_file, '<element nr="%d" type="density" design="%g"/>\n', elements(i,1), density(i));
end

fprintf(xml_file, '</set>\n</cfsErsatzMaterial>');

fclose(xml_file);
