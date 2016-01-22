function [ file ] = matrixToMeshAndDensity(A,filename,shearingAngle)
% MATRIXTOMESHORDENSITY  -  Generates eihter a (sparse) ANSYS mesh file or
% a density file and a full ANSYS mesh file out of a given density matrix.
% 
% This is due to the fact, that we cannot apply periodic boundary
% conditions, if there are no nodes at the boundary. In this case the
% version with a density file is used.
%
% @param:
%       A               density matrix (entry has to be 0 for no material)
%       meshfile        name of generated .mesh file
%       shearingAngle   shearing angle of mesh in (-pi/2,pi/2)
%                       (optional, default = 0)
%

if nargin < 3
    shearingAngle = 0.0;
end

[n,m] = size(A);

% Check if matrix is periodic
msg = 'Matrix has no periodic boundary';
assert( all(A(1,:) == A(n,:)), msg )
assert( all(A(:,1) == A(:,m)), msg )

file = [filename,'.dens'];
Homogenization.matrixToMesh(ones(n,m),[filename,'.mesh'],shearingAngle);
Homogenization.matrixToDensity(A,file);
