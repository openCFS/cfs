The code in this directory allows to generate (optimization) densitiy.xml files
from arbitrary geometries exported from GiD. As GiD can read some stuff this is
quite generic. It has been written by Michael Stingl.

----------------------

Hallo Fabian,

ich haab die Matlab Software gesten abend fetiggemacht.
Wir können jetzt mesh files und geo files (für die graphen)  in Matlab einlesen
und mittels einfacher Geometriedefinitionen density files erzeugen.
Ich hab mich mittels einfacher Visualisierung auch schon teils von der korrekten Funktionsweise überzeugt.

Das Kugelbeispiel hab ich angehängt.

Anbei schicke ich nochmal die eingelesene Geometrie von Jan. Ist doch korrekt oder?
Im file jan_box.txt  ist die Boundingbox für den Graphen, den ich von dir bekommen habe enthalten.

Sobald ich ein zugeh. Mesh habe, kann ich sofort das density file erzeugen.

Ich füge auch mal die Matlab codes bei.

Aufruf (Kugel):

list=mshlist(); % erzeugt Liste mit vorhandenen Mesh files
[nn,ee]=readmesh(list{4}); % liest Knoten (nn) und Elemene(ee) ein
circls = [.5 .5 .5 .2]; % Definiert Kugel um (.5,.5,.5) mit Radius .2
ddd=createdens(nn,ee,circls); % erzeugt density Daten
writexml('cube_20.xml',ee,ddd); % schreibt zugehöriges XML file
plotstruct(nn,ee,ddd); % Visualisiert Kugel ...


Einlesen von Geometriedaten (Jan):

listg=geolist;
[nn1,edges]=readgeo(listg{1});
createcylinders(nn1,edges,radius);

...

Gruss,
Michael

