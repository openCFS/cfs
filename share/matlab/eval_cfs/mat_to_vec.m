function vec = mat_to_vec(mat)
%MAT_TO_VEC (V1.0) create a vector from a matrix
[n,m] = size(mat);
vec = zeros(n*m,1);
for ii=1:n
    for uu=1:m
        % CHANGE
        vec((ii-1)*n+uu) = mat(uu,ii);
    end
end
end