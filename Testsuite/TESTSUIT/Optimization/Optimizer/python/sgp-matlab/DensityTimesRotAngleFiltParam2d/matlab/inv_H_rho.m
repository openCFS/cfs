function [H] = inv_H_rho(H0_inv, rho, p)
  % multiply core stiffness matrix with each entry of vector (dens^p)
  fact = 1.0./(rho.^p);
      
  %(\rho^p*stiff)^(-1) = (\rho^p)^(-1)*stiff^(-1)
  H = bsxfun(@times,reshape(fact,1,1,numel(fact)),H0_inv);
end