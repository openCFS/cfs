function pp= computebox(nn1)

minx=min(nn1(:,2));
maxx=max(nn1(:,2));
miny=min(nn1(:,2));
miny=min(nn1(:,3));
maxy=max(nn1(:,3));
minz=min(nn1(:,4));
maxz=max(nn1(:,4));
pp(1,:) = [minx, miny, minz];
pp(2,:) = [maxx, miny, minz];
pp(3,:) = [minx, maxy, minz];
pp(4,:) = [maxx, maxy, minz];
pp(5,:) = [minx, miny, maxz];
pp(6,:) = [maxx, miny, maxz];
pp(7,:) = [minx, maxy, maxz];
pp(8,:) = [maxx, maxy, maxz];
pp=pp';
plot3([pp(1,1) pp(1,2) pp(1,3) pp(1,4) pp(1,5) pp(1,6) pp(1,7) pp(1,8)], [pp(2,1) pp(2,2) pp(2,3) pp(2,4) pp(2,5) pp(2,6) pp(2,7) pp(2,8)], [pp(3,1) pp(3,2) pp(3,3) pp(3,4) pp(3,5) pp(3,6) pp(3,7) pp(3,8)], '*');
