function [vol] = calc3Dvol(s1,s2,s3)

vol = s1^2+s2^2+s3^2;
if s1 >= s2 && s1 >= s3
   dvol = -s1*s3^2-s1*s2^2; 
elseif s2 >= s1 && s2 >= s3
   dvol = -s2*s1^2-s2*s3^2;   
else
   dvol = -s3*s1^2-s3*s2^2;
end
vol = vol + dvol;

end
