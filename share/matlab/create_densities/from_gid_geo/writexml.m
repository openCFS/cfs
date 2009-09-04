function writexml( filename, elements, density )

xml_file = fopen(filename, 'w');

nelements = size(elements,1);

for i=1:nelements
   fprintf(xml_file, '<element nr="%d" type="density" design="%g"/>\n', elements(i,1), density(i));
end


fclose(xml_file);