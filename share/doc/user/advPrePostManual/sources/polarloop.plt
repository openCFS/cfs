a = a+1
set trange [-pi:pi]
plot a*0.2*cos(t),a*0.2*sin(t) lt 5

set trange [-1:1]
plot t*cos(a*pi/4),t*sin(a*pi/4) lt 4

if(a<5) reread
