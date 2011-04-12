% 'msh/unit-cube-60.msh'
% 'msh/unit-cube-80.msh'
mshfile = '~/meshes/cube-30.mesh';

% 'geo/jan-cube.geo'
% 'geo/jan-diamond.geo'
geofile = 'geo/jan-diamond.geo';

disp(datestr(now));

disp(['reading ', mshfile]);
[nn,ee]=readgidmesh(mshfile);

disp(['reading ', geofile]);
[nn1,ee1]=readgeo(geofile);


disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.09);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.09.density.xml',ee,ret);

if 0 == 1
disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.07);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.07.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.11);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.11.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.13);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.13.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.15);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.15.density.xml',ee,ret);


disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.06);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.06.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.08);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.08.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.10);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.10.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.12);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.12.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.14);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.14.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.16);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.16.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.18);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.18.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.20);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.20.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.22);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.22.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.24);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.24.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.26);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.26.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.28);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.28.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.30);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-diamond30-0.30.density.xml',ee,ret);
end
