load Eh11
load Eh12
load Eh13
load Eh22
load Eh23
load Eh33

n = size(Eh11,3);
a = zeros(1,n);
b = zeros(1,n);
c = zeros(1,n);
d = zeros(1,n);
e = zeros(1,n);
f = zeros(1,n);

for i=1:n
    a(i) = Eh11(3,3,i);
    b(i) = Eh12(3,3,i);
    c(i) = Eh13(3,3,i);
    d(i) = Eh22(3,3,i);
    e(i) = Eh23(3,3,i);
    f(i) = Eh33(3,3,i);
end

% Filter
iN = 1;
a = filter(ones(1,iN)/iN,1,a);
b = filter(ones(1,iN)/iN,1,b);
c = filter(ones(1,iN)/iN,1,c);
d = filter(ones(1,iN)/iN,1,d);
e = filter(ones(1,iN)/iN,1,e);
f = filter(ones(1,iN)/iN,1,f);

range = -pi/2:pi/(n-1):pi/2;
figure; hold on;
plot(range,a,'b-');
plot(range,b,'r-');
plot(range,c,'g-');
plot(range,d,'y-');
plot(range,e,'k-');
plot(range,f,'m-');
axis tight
legend('Eh11','Eh12','Eh13','Eh22','Eh23','Eh33')
hold off
