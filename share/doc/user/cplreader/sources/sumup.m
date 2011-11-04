f = fopen('out.bin', "r");
% Specify here the length of the vector from the HDF5 file.
sizeVec = 10;
vec = fread (f, sizeVec, "double", 0);

function s = sumUp(m,n)
s = 0;
for i=1:n
  s = s + m(i);
end
endfunction

sumUp(vec,sizeVec)
