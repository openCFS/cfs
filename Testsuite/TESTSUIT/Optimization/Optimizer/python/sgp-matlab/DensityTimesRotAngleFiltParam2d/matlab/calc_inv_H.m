%% for given H(rho), compute H(rho,alpha)^(-1) for all combinations of rho and alpha
% in: H - size 3 x 3 x n_elems^2
%     alpha - size n_elems x 1
function Hi = calc_inv_H(H,alpha)
  na = size(alpha,1);
  % 2D: 3 x 3: 3D: 6 x 6
  rows = size(H,1);
  Hi = zeros(rows,rows,na,na);
  
  for n = 1:na
    Qinv = get_inv_rot_mat(alpha(n));
    for r = 1:na
      Hi(:,:,r,n) = Qinv'*H(:,:,r)*Qinv;
    end
  end
  
  Hi = reshape(Hi, [3 3 na*na]);
  Hi = permute(Hi, [3 1 2]);
end