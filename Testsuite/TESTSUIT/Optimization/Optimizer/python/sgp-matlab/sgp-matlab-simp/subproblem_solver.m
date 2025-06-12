classdef subproblem_solver < handle
  properties (Access = private)
    xl = 1e9;
    xu = 1;
    
    % dimension of the problem
    dim = 2;
    % bounds of mesh
    nx = 0;
    ny = 0;
    nz = 0;
    % number of design elements
    n_elems = 0;
    % number of grid levels
    levels = 1;
    % number of parameter samples per grid
    samples = 5;
    % lower asymptote
    L = 0;
    % penalty for regularization
    p_filt = 0;
    % we assume density filtering
    % used for exact approximation of filter
    F = -1;
    dAFii = 0;
    dAF = 0;
    
    %% assume SIMP
    % core stiffness tensor and its inverse
    core_stiffness = 0;
    core_stiffness_inv = 0;
    % exponent for power law
    p = 1;
    % store all density values on finest grid level + x^2 in last entry
    x_all = 0;
    % precomputed arrays
    H_all = 0;
    % inverse of H_all
    H_all_i = 0;
  end
  
  methods
    function this = subproblem_solver(nx, ny, nz, xl, xu, levels, samples, L, p_filt, rad, stiffness, p)
      disp("entering MATLAB")
      addpath(genpath('mtimesx-master'));
      
      assert(samples >= 0);
      assert(levels >= 0);
      assert(rad >= 0);
      assert(nx > 0 && ny > 0 && nz > 0);
      
      if nz > 1
        this.dim = 3;
      end
      
      this.nx = nx; this.ny = double(ny); this.nz = double(nz);
      
      this.n_elems = nx * ny * nz;
      
      this.xl = xl;
      this.xu = xu;
      this.levels = levels;
      this.samples = samples;
      this.L = L;
      this.p_filt = p_filt;
      this.core_stiffness = stiffness;
      % need inverse (assume L=0) for subproblems
      this.core_stiffness_inv = stiffness^(-1);
      this.p = p;
      fprintf("MATLAB:  xl: %2.9f, xu: %2.6f, levels: %2.6f, samples: %2.6f, p_filt: %2.6f, rad: %2.6f, p: %2.6f\n", this.xl, this.xu, this.levels, this.samples, this.p_filt, rad, this.p);
      
      % all density values on finest grid
      x = (this.xl:(this.xu-this.xl)/(this.samples^this.levels):this.xu)';
      H = this.eval_stiffness_tensor(x);
      Hi = this.eval_inv_stiffness_tensor(x);
      % store only non-zero values
      % order: H_11, H_22, H_33, H_12, H_13, H_23, vol, penbw=0
      % don't store vol 
      vol = x;
      if this.dim == 2
        this.H_all = [squeeze(H(1,1,:)) squeeze(H(2,2,:)) squeeze(H(3,3,:)) squeeze(H(1,2,:)) squeeze(H(1,3,:)) squeeze(H(2,3,:))];
        this.H_all_i = [squeeze(Hi(1,1,:)) squeeze(Hi(2,2,:)) squeeze(Hi(3,3,:)) squeeze(Hi(1,2,:)) squeeze(Hi(1,3,:)) squeeze(Hi(2,3,:)) vol(:) zeros(size(vol))];
      else
        this.H_all = [squeeze(H(1,1,:)) squeeze(H(2,2,:)) squeeze(H(3,3,:)) squeeze(H(4,4,:)) squeeze(H(5,5,:)) squeeze(H(6,6,:)) squeeze(H(1,2,:)) squeeze(H(1,3,:)) squeeze(H(2,3,:))];
        this.H_all_i = [squeeze(Hi(1,1,:)) squeeze(Hi(2,2,:)) squeeze(Hi(3,3,:)) squeeze(Hi(4,4,:)) squeeze(Hi(5,5,:)) squeeze(H(6,6,:)) squeeze(Hi(1,2,:)) squeeze(Hi(1,3,:)) squeeze(Hi(2,3,:)) vol(:) zeros(size(vol))];
      end
      
      % need this structure for evaluation of f_merit of model
      this.x_all = [x(:) x.^2];
      
      % compute filter matrix
      this.F = this.calc_filter(rad);
     
      tmp = speye(this.n_elems);
      % Q
      this.dAF = (tmp-this.F)'*(tmp-this.F); % bei Uebergabe noch halbieren ...
      % diag(Q)
      this.dAFii = diag(this.dAF);
    end
    
    function F = calc_filter(this, radius)
      % check if filter was already computed
      if this.F ~=-1
        % matlab engine does not support sparse matrices as return type to python
        F = full(this.F);
        return;
      end
      
      if isequal(radius,0)
        F = speye(this.n_elems,this.n_elems);
        return;
      end
      
      %% prepare filter
      iH  = ones(this.n_elems*(2*(ceil(radius)-1)+1)^2,1);
      jH  = ones(size(iH));
      sH  = zeros(size(iH));
      k = 0;
      
      for k1 = 1:this.nz
        for i1 = 1:this.nx
          for j1 = 1:this.ny
            e1 = (k1-1)*this.nx*this.ny + (i1-1)*this.ny+j1;
            for k2 = max(k1-(ceil(radius)-1),1):min(k1+(ceil(radius)-1), this.nz)
              for i2 = max(i1-(ceil(radius)-1),1):min(i1+(ceil(radius)-1), this.nx)
                for j2 = max(j1-(ceil(radius)-1),1):min(j1+(ceil(radius)-1), this.ny)
                  e2 = (k2-1)*this.nx*this.ny + (i2-1)*this.ny+j2;
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
      tmp  = spdiags(1./sum(H,2),0,double(this.n_elems),double(this.n_elems));
      F = H*tmp;
    end
    
    %% returns interpolated stiffness tensor/matrix
    % dens: density value
    
    function [C] = eval_stiffness_tensor(this, dens)
      % multiply core stiffness matrix with each entry of vector (dens^p)
      % in 2d: C is a 3 x 3 x n_elems^2 matrix
      fact = dens.^this.p;
      
      C = bsxfun(@times,reshape(fact,1,1,numel(fact)),this.core_stiffness);
    end
    
    function [C] = eval_inv_stiffness_tensor(this, dens)
      % multiply core stiffness matrix with each entry of vector (dens^p)
      % in 2d: C is a 3 x 3 x n_elems^2 matrix
      fact = 1.0./(dens.^this.p);
      
      % (\rho^p*stiff)^(-1) = (\rho^p)^(-1)*stiff^(-1)
      C = bsxfun(@times,reshape(fact,1,1,numel(fact)),this.core_stiffness_inv);
    end
    
    
    function [x_new,f_min] = solve(this,x,dfH,lambda)
      %fprintf(1,"MATLAB solve(): size of H= %d %d %d; size of dfH= %d %d %d\n",size(H),size(dfH));
      f_min = 0;
      
      % permute H and dfH to conform with input of mtimesx:
      % orig size of H: n_elems x 3 x 3  -> input of mtimesx: 3 x 3 x n_elems
      % dA = df(H^(k))/dH
      dA = permute(dfH, [2,3,1]);
      
      x_old = x(1:this.n_elems);
      %disp(x_old)
      % assume L=0 as this was also assumed in (H^(k)-L)^(-1)
      % dLA = H^(k) - L
      dLA = this.eval_stiffness_tensor(x_old);
      %fprintf(1,"size dA= %d %d %d; size dLA = %d %d %d \n",size(dA),size(dLA));
      dAL = mtimesx(dA,dLA);
      
      % -(H^(k)-L) * dfH * (H^(k)-L)
      BBglobA = -mtimesx(dLA,dAL);
      
      n_all_samples = this.samples^this.levels;
      m_grid_bisect = this.samples;
      grid_bisect = 0:m_grid_bisect;
      
      %fprintf("MATLAB size this.dAF=%d %d; size x=%d %d; size this.dAFii=%d %d\n",size(this.dAF),size(x),size(this.dAFii));
      % Q*x_bar - diag(Q)*x_bar
      % Q = (I - AF)' * (I-AF) with filter matrix AF
      dAFiR1 = this.dAF*x_old-x_old.*this.dAFii;
      
      jstart = ones(this.n_elems,1);
      
      for q=1:this.levels
        ddx = this.samples^(this.levels-q);
        jjall = bsxfun(@plus,jstart,ddx*grid_bisect);
        
        Ftmpall = zeros(this.n_elems,m_grid_bisect+1);
        for i=1:this.n_elems
          idx = jjall(i,:);
          HALLi = this.H_all_i(idx,:);
          
          if this.dim == 3
            % assume orthotropic tensor
            BBBi = [BBglobA(1,1,i); BBglobA(2,2,i); BBglobA(3,3,i); ...
              BBglobA(4,4,i); BBglobA(5,5,i); BBglobA(6,6,i); ...
              BBglobA(1,2,i); BBglobA(1,3,i); BBglobA(2,3,i); ...
              lambda; 0 ];
          else % dim = 2
            BBBi = [BBglobA(1,1,i); BBglobA(2,2,i); BBglobA(3,3,i); ...
              BBglobA(1,2,i); BBglobA(1,3,i); BBglobA(2,3,i); ...
              lambda; 0 ];
          end
          
          if this.p_filt > 0
            XALL = this.x_all(idx,:);
            dAFiRi = [dAFiR1(i); 0.5*this.dAFii(i) ];
            Ftmpall(i,:) = HALLi*BBBi + this.p_filt*XALL*dAFiRi;
          else
            Ftmpall(i,:) = HALLi*BBBi;
          end
          
        end
        
        % get minimum for each row of Ftmpall
        [f_min,jmin] = min(Ftmpall,[],2);
        
        %disp(Ftmpall)
        %disp(min(jmin))
        
        % compute indices of found min values on finest grid
        jval_tmp = jstart + ddx*grid_bisect(jmin)';
        
        % set lower interval bound for next level
        % special cases: we're already at lower or upper design value bound
        i_jstart = jmin-1;   
        i_jstart(i_jstart < 1) = 1;
        i_jstart(i_jstart > m_grid_bisect-1) = m_grid_bisect -1;
        
        % valid for next level
        jstart = jstart + ddx*grid_bisect(i_jstart)';
        
        % e.g. jmin = 2
        % next level: left = 1 and right = 3
        % -> factor 2
        m_grid_bisect = 2*this.samples;
        grid_bisect = 0:m_grid_bisect;
        
      end
      
      %disp(min(jval_tmp))
      x_new = this.x_all(jval_tmp,1)';
            
    end % function
    
  end % class methods
  
end % class

