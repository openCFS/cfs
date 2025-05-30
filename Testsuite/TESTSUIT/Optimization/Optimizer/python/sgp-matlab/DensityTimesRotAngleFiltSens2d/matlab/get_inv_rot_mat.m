function [Q]=get_inv_rot_mat(alpha)
  Q = zeros(3,3);
  % computed with symbolic toolbox
  Q(1,1) = 1 - sin(alpha)^2;
  Q(1,2) = sin(alpha)^2;
  Q(1,3) = sin(2*alpha);
  
  Q(2,1) = sin(alpha)^2;
  Q(2,2) = 1 - sin(alpha)^2;
  Q(2,3) = -sin(2*alpha);
  
  Q(3,1) = -sin(2*alpha)/2;
  Q(3,2) = sin(2*alpha)/2;
  Q(3,3) = cos(2*alpha);
end