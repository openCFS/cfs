function [Q] = get_rotation_matrix(dim,alpha,beta,gamma)
  if nargin == 2
    beta = 0;
    gamma = 0;
  end
  
  % for numerical comparisons
  eps = 1e-10;
  assert(dim == 2 || dim == 3);
  %assert(-eps < alpha < 2*pi+eps);
  %assert(-eps < beta < 2*pi+eps);
  %assert(-eps < gamma < 2*pi+eps);
  
  if dim == 2
    % output is a cell of cells
    Q = arrayfun(@get_rotation_matrix_3x3_2d, alpha, 'un',0);
    % convert to array of size: 3 x 3 x nelems
    %fprintf(1,'size Q: %d %d; size alpha: %d %d', size(Q), size(alpha))
    Q = reshape(cell2mat(Q),3,3,size(Q,1));
  else  
    Q = get_rotation_matrix_6x6(alpha,beta,gamma);
  end
end


function [Q_inv] = get_rotation_matrix_inv(dim,alpha,beta,gamma)
  if nargin == 2
    beta = 0;
    gamma = 0;
  end
  
  % for numerical comparisons
  eps = 1e-10;
  assert(dim == 2 || dim == 3);
  %assert(-eps < alpha < 2*pi+eps);
  %assert(-eps < beta < 2*pi+eps);
  %assert(-eps < gamma < 2*pi+eps);
  
  if dim == 2
    Q_inv = arrayfun(@inv_rotation_matrix_3x3_2d, alpha, 'un',0)';
    % convert to array of size: 3 x 3 x nelems
    Q_inv = reshape(cell2mat(Q_inv),3,3,size(Q_inv,2));
    
    % Q * Q_inv = I (identity)
    assert(max(Q(1)*Q_inv(1)-eye(3,3),[],'all') < eps);
  else  
    disp('inverted rot matrix not implemented for 3D!');
    assert(false);
  end
end

function R = get_rot_2x2(angle)
  R = zeros(2,2);

  R(1,1) =  cos(angle);
  R(1,2) =  -sin(angle);
  R(2,1) =  sin(angle);
  R(2,2) =  cos(angle);
end

function Q = get_rotation_matrix_3x3_2d(alpha)
  R = get_rot_2x2(alpha);

  Q = zeros(3,3);

  Q(1,1) = R(1,1)*R(1,1);
  Q(1,2) = R(1,2)*R(1,2);
  Q(1,3) = 2.0*R(1,1)*R(1,2);

  Q(2,1) = R(2,1)*R(2,1);
  Q(2,2) = R(2,2)*R(2,2);
  Q(2,3) = 2.0*R(2,1)*R(2,2);

  Q(3,1) = R(1,1)*R(2,1);
  Q(3,2) = R(1,2)*R(2,2);
  Q(3,3) = R(1,1)*R(2,2) + R(1,2)*R(2,1);
end

function Q = inv_rotation_matrix_3x3_2d(alpha)
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

function Q = get_rotation_matrix_6x6(alpha,beta,gamma)
  R = zeros(3,3);
  R(1,1) =  cos(beta) * cos(gamma);
  R(1,2) = -cos(beta) * sin(gamma);
  R(1,3) =  sin(beta);
  R(2,1) =  cos(alpha)*sin(gamma) + sin(alpha)*sin(beta)*cos(gamma);
  R(2,2) =  cos(alpha)*cos(gamma) - sin(alpha)*sin(beta)*sin(gamma);
  R(2,3) = -sin(alpha)*cos(beta);
  R(3,1) =  sin(alpha)*sin(gamma) - cos(alpha)*sin(beta)*cos(gamma);
  R(3,2) =  sin(alpha)*cos(gamma) + cos(alpha)*sin(beta)*sin(gamma);
  R(3,3) =  cos(alpha)*cos(beta);
  
  Q = zeros(6,6);
  
  Q(1,1) = R(1,1) * R(1,1);
  Q(1,2) = R(1,2) * R(1,2);
  Q(1,3) = R(1,3) * R(1,3);
  Q(1,4) = 2.0 * R(1,2) * R(1,3);
  Q(1,5) = 2.0 * R(1,1) * R(1,3);
  Q(1,6) = 2.0 * R(1,1) * R(1,2);

  Q(2,1) = R(2,1) * R(2,1);
  Q(2,2) = R(2,2) * R(2,2);
  Q(2,3) = R(2,3) * R(2,3);
  Q(2,4) = 2.0 * R(2,2) * R(2,3);
  Q(2,5) = 2.0 * R(2,1) * R(2,3);
  Q(2,6) = 2.0 * R(2,1) * R(2,2);

  Q(3,1) = R(3,1) * R(3,1);
  Q(3,2) = R(3,2) * R(3,2);
  Q(3,3) = R(3,3) * R(3,3);
  Q(3,4) = 2.0 * R(3,2) * R(3,3);
  Q(3,5) = 2.0 * R(3,1) * R(3,3);
  Q(3,6) = 2.0 * R(3,1) * R(3,2);

  Q(4,1) = R(2,1) * R(3,1);
  Q(4,2) = R(2,2) * R(3,2);
  Q(4,3) = R(2,3) * R(3,3);
  Q(4,4) = R(2,2) * R(3,3) + R(2,3) * R(3,2);
  Q(4,5) = R(2,1) * R(3,3) + R(2,3) * R(3,1);
  Q(4,6) = R(2,1) * R(3,2) + R(2,2) * R(3,1);

  Q(5,1) = R(1,1) * R(3,1);
  Q(5,2) = R(1,2) * R(3,2);
  Q(5,3) = R(1,3) * R(3,3);
  Q(5,4) = R(1,2) * R(3,3) + R(1,3) * R(3,2);
  Q(5,5) = R(1,1) * R(3,3) + R(1,3) * R(3,1);
  Q(5,6) = R(1,1) * R(3,2) + R(1,2) * R(3,1);

  Q(6,1) = R(1,1) * R(2,1);
  Q(6,2) = R(1,2) * R(2,2);
  Q(6,3) = R(1,3) * R(2,3);
  Q(6,4) = R(1,2) * R(2,3) + R(1,3) * R(2,2);
  Q(6,5) = R(1,1) * R(2,3) + R(1,3) * R(2,1);
  Q(6,6) = R(1,1) * R(2,2) + R(1,2) * R(2,1);
end
