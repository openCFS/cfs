%% for one dens and one alpha
function [H] = calc_H(H0,dens,alpha,p)
  Q = get_rotation_matrix(2,alpha,0,0);
  H = dens.^p * mtimesx(Q,mtimesx(H0,Q'));
end