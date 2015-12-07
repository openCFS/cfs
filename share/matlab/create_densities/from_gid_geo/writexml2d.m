function writexml2d( filename, elements, density )

xml_file = fopen(filename, 'w');

fprintf(xml_file, '<?xml version="1.0"?>\n<cfsErsatzMaterial writeCounter="16">\n <header>\n  <design constant="false" initial="0.5" name="density" region="all" scale="false"/>\n  <transferFunction application="mech" design="density" param="1.0" type="simp"/>\n </header>\n<set id="struct">\n');


nelements = size(elements,1);

projdens=zeros(3600,1);
for i=1:3600
  projdens(i)=projdens(i)+0.001;
end

for s=1:60
  for i=1:3600
    ind=i + 3600*(s-1);
    if density(ind) > projdens(i)
      projdens(i) = density(ind);
    end
  end
end

for i=1:3600
   fprintf(xml_file, '<element nr="%d" type="density" design="%g"/>\n', elements(i,1), projdens(i));
end

fprintf(xml_file, '</set>\n</cfsErsatzMaterial>');

fclose(xml_file);
