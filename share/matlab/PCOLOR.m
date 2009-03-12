function ch = PCOLOR( A )

% PCOLOR: Make pseudocolor plot of a matrix (returns handle for colorbar)

[m,n] = size(A);
AA = zeros(m+1,n+1);
AA(1:m,1:n) = A;
pcolor(AA);
axis square;
axis ij;
ch = colorbar;
