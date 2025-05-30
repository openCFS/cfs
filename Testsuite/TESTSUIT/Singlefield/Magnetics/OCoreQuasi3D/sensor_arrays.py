from params import *;
# python snippet to write circle coordinates in appropriate format 
# to define CFS sensorArray output.

from numpy import sin,cos,arange, array,pi,zeros_like,savetxt,sqrt;


dr = arange(0.0,r_domain,5e-3); 
P = array([dr*cos(pi/4),dr*sin(pi/4),zeros_like(dr)]);
savetxt('sensor_arrays.csv',P.T,delimiter=',');
