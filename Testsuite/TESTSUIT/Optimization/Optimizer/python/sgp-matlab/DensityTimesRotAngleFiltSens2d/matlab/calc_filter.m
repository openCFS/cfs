function [F] = calc_filter(nx, ny, nz, radius)
  nelems = nx * ny * nz;
  if isequal(radius,0)
    F = speye(nelems,nelems);
    return;
  end

  %% prepare filter
  iH  = ones(nelems*(2*(ceil(radius)-1)+1)^2,1);
  jH  = ones(size(iH));
  sH  = zeros(size(iH));
  k = 0;

  for k1 = 1:nz
    for i1 = 1:nx
      for j1 = 1:ny
        e1 = (k1-1)*nx*ny + (i1-1)*ny+j1;
        for k2 = max(k1-(ceil(radius)-1),1):min(k1+(ceil(radius)-1), nz)
          for i2 = max(i1-(ceil(radius)-1),1):min(i1+(ceil(radius)-1), nx)
            for j2 = max(j1-(ceil(radius)-1),1):min(j1+(ceil(radius)-1), ny)
              e2 = (k2-1)*nx*ny + (i2-1)*ny+j2;
              k = k+1;
              iH(k) = e1;
              jH(k) = e2;
              sH(k) = max(0,radius-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
            end
          end
        end
      end
    end
  end

  %% write filter to sparse matrix
  H    = sparse(iH,jH,sH);
  tmp  = spdiags(1./sum(H,2),0,double(nelems),double(nelems));
  F = H*tmp;
end