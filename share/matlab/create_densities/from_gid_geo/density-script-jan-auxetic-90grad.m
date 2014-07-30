mshfile = 'msh/struct-80.mesh';

% 'geo/auxetic.geo'
geofile = 'geo/auxetic-90grad.geo'

disp(datestr(now));

disp(['reading ', mshfile]);
[nn,ee]=readgidmesh(mshfile);

disp(['reading ', geofile]);
[nn1,ee1]=readgeo(geofile);

if 0 == 1

  disp(datestr(now));
  cyl=createcylinders(nn1,ee1,0.45);
  ret=createdens( nn, ee, 0, [], cyl);
  writexml('jan-auxetic80-90grad-0.45.density.xml',ee,ret);

end

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.55);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-auxetic80-90grad-0.55.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.65);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-auxetic80-90grad-0.65.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.75);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-auxetic80-90grad-0.75.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.85);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-auxetic80-90grad-0.85.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,0.95);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-auxetic80-90grad-0.95.density.xml',ee,ret);

disp(datestr(now));
cyl=createcylinders(nn1,ee1,1.05);
ret=createdens( nn, ee, 0, [], cyl);
writexml('jan-auxetic80-90grad-1.05.density.xml',ee,ret);

if 1 == 1

  disp(datestr(now));
  cyl=createcylinders(nn1,ee1,0.40);
  ret=createdens( nn, ee, 0, [], cyl);
  writexml('jan-auxetic80-90grad-0.40.density.xml',ee,ret);

  disp(datestr(now));
  cyl=createcylinders(nn1,ee1,0.50);
  ret=createdens( nn, ee, 0, [], cyl);
  writexml('jan-auxetic80-90grad-0.50.density.xml',ee,ret);

  disp(datestr(now));
  cyl=createcylinders(nn1,ee1,0.60);
  ret=createdens( nn, ee, 0, [], cyl);
  writexml('jan-auxetic80-90grad-0.60.density.xml',ee,ret);

  disp(datestr(now));
  cyl=createcylinders(nn1,ee1,0.70);
  ret=createdens( nn, ee, 0, [], cyl);
  writexml('jan-auxetic80-90grad-0.70.density.xml',ee,ret);

  disp(datestr(now));
  cyl=createcylinders(nn1,ee1,0.80);
  ret=createdens( nn, ee, 0, [], cyl);
  writexml('jan-auxetic80-90grad-0.80.density.xml',ee,ret);

  disp(datestr(now));
  cyl=createcylinders(nn1,ee1,0.90);
  ret=createdens( nn, ee, 0, [], cyl);
  writexml('jan-auxetic80-90grad-0.90.density.xml',ee,ret);

  disp(datestr(now));
  cyl=createcylinders(nn1,ee1,1.0);
  ret=createdens( nn, ee, 0, [], cyl);
  writexml('jan-auxetic80-90grad-1.00.density.xml',ee,ret);

end



if 1 == 1

  disp(datestr(now));
  cyl=createcylinders(nn1,ee1,1.10);
  ret=createdens( nn, ee, 0, [], cyl);
  writexml('jan-auxetic80-90grad-1.10.density.xml',ee,ret);

  disp(datestr(now));
  cyl=createcylinders(nn1,ee1,1.15);
  ret=createdens( nn, ee, 0, [], cyl);
  writexml('jan-auxetic80-90grad-1.15.density.xml',ee,ret);

  disp(datestr(now));
  cyl=createcylinders(nn1,ee1,1.20);
  ret=createdens( nn, ee, 0, [], cyl);
  writexml('jan-auxetic80-90grad-1.20.density.xml',ee,ret);

  disp(datestr(now));
  cyl=createcylinders(nn1,ee1,1.25);
  ret=createdens( nn, ee, 0, [], cyl);
  writexml('jan-auxetic80-90grad-1.25.density.xml',ee,ret);

  disp(datestr(now));
  cyl=createcylinders(nn1,ee1,1.30);
  ret=createdens( nn, ee, 0, [], cyl);
  writexml('jan-auxetic80-90grad-1.30.density.xml',ee,ret);

end
