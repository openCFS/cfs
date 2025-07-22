classdef subproblem_solver < handle
  properties (Access = private)
    % lower + upper bound for density variable
    xl = 1e9;
    xu = 1;
    
    % lower + upper bound for angle variable
    al = 0;
    au = 2*pi;
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
    % filter radius
    rad_filt = 0;
    % filter sensitivity
    filt_sens = 1;
    
    %% assume SIMP
    % core stiffness tensor H0 and its inverse
    H0 = 0;
    H0_inv = 0;
    % exponent for power law
    p = 1;
    % store all density values on finest grid level + x^2 in last entry
    x_all = 0;
    % store all angle values on finest grid level
    a_all = 0;
    % need this for regularization term
    % store all design parameter combinations + squared sum
    design_all = 0;
    % precomputed arrays
    % rotated core stiffness tensors - density not included yet!
    H_all = 0;
    % inverse of H_all
    H_all_i = 0;
    vol = 0;
  end
  
  methods
    function this = subproblem_solver(nx, ny, nz, xl, xu, al, au, levels, samples, L, filt_sens, p_filt, rad, stiffness, p)
      disp("entering MATLAB")
      addpath(genpath('mtimesx-master'));
      
      assert(samples >= 0);
      assert(levels >= 0);
      assert(rad >= 0);
      assert(nx > 0 && ny > 0 && nz > 0);
      
      if nz > 1
        this.dim = 3;
      end
      
      this.filt_sens = filt_sens;
      
      this.nx = nx; this.ny = double(ny); this.nz = double(nz);
      
      this.n_elems = nx * ny * nz;
      
      % lower + upper bound for density variable
      this.xl = xl;
      this.xu = xu;
      % lower + upper bound for angle variable
      this.al = al;
      this.au = au;
      this.levels = levels;
      this.samples = samples;
      this.L = L;
      this.p_filt = p_filt;
      this.H0 = stiffness;
      % need inverse (assume L=0) for subproblems
      this.H0_inv = inv(stiffness);
      this.p = p;
      this.rad_filt = rad;
      fprintf("MATLAB:  xl: %2.9f, xu: %2.6f, al: %2.9f, au: %2.6f, levels: %2.6f, samples: %2.6f, p_filt: %2.6f, rad: %2.6f, p: %2.6f\n", this.xl, this.xu, this.al, this.au, this.levels, this.samples, this.p_filt, rad, this.p);
      if this.filt_sens == 1
        fprintf("using sensitivity filtering")
      end
      
      % all density values on finest grid
      n_all_samples = (this.samples)^this.levels;
      x = (this.xl:(this.xu-this.xl)/(n_all_samples):this.xu)';
      this.x_all = x;
      % all angle values on finest grid
      alpha = (this.al:(this.au-this.al)/(n_all_samples):this.au)';
      this.a_all = alpha;
      
      % need this structure for evaluation of f_merit of model
      % store all possible combinations of dens and alpha (column major!)
      % dim: #{samples on finest grid}^2 x 1
      tmp_x = repmat(x,size(x,1),1);
      % dim: #{samples on finest grid}^2 x 1
      tmp_alpha = reshape(repmat(alpha,1,size(alpha,1))',[size(tmp_x,1) 1]);
      % dim: #{samples on finest grid}^2 x 1
      square = tmp_x.^2 + tmp_alpha.^2;
      this.design_all = [tmp_x(:) tmp_alpha(:)  square(:)];
      %this.design_all = [tmp_x(:) tmp_x.^2];
      
      % compute stiffness tensors with power law
      % H_rho_inv is a 3 x 3 x n_elems matrix
      H_rho_inv = inv_H_rho(this.H0_inv,x,this.p);
      assert(size(H_rho_inv,1) == 3);
      assert(size(H_rho_inv,2) == 3);
      assert(size(H_rho_inv,3) == size(x,1));
      
      % for each density, apply inverse rotation matrices on H_rho_inv
      % H(rho,alpha) = Q(alpha)^(-T) H(rho)^(-1) Q(alpha)^(-1)
      % Hi has dim {samples on finest grid}^2 x 3 x 3
      Hi = calc_inv_H(H_rho_inv,alpha);
      
      % store only non-zero values
      % order: H_11, H_22, H_33, H_12, H_13, H_23, vol, penbw=0
      vol = repmat(x,1,size(alpha,1));
      this.vol = vol;
      if this.dim == 2
        % for all x and angles: [H11 H22 H33 H12 H13 H23]
        % -> size of H_all: #{sample on finest grid} x 6*#{sample on finest grid}
        % squeeze(H_rot(1,1,:,:)) has dimensions: #{sample on finest grid} x *#{sample on finest grid}
        % -> 1st entry: density, 2nd entry: angle
        %this.H_all = [squeeze(H_rot(:,1,1)) squeeze(H_rot(:,2,2)) squeeze(H_rot(:,3,3)) squeeze(H_rot(:,1,2)) squeeze(H_rot(:,1,3)) squeeze(H_rot(:,2,3))];
        %this.H_all = H_rot;
        % inverse matrix and volume are independent on the angle
        % -> size of H_all: #{sample on finest grid} x 8
        this.H_all_i = [Hi(:,1,1) Hi(:,2,2) Hi(:,3,3) Hi(:,2,3) Hi(:,1,3) Hi(:,1,2) vol(:) zeros(size(vol(:)))];
      else
        %this.H_all = [squeeze(H_rot(:,1,1)) squeeze(H_rot(:,2,2)) squeeze(H_rot(:,3,3)) squeeze(H_rot(:,4,4)) squeeze(H_rot(:,5,5)) squeeze(H_rot(:,6,6)) squeeze(H_rot(:,1,2)) squeeze(H_rot(:,1,3)) squeeze(H_rot(:,2,3))];
        this.H_all_i = [Hi(:,1,1) Hi(:,2,2) Hi(:,3,3) Hi(:,4,4) Hi(:,5,5) Hi(:,6,6) Hi(:,2,3) Hi(:,1,3) Hi(:,1,2) vol(:) zeros(size(vol(:)))];
      end
      
      % compute filter matrix
      this.F = this.get_filter(rad);
      %this.F = mmread('/home/bvu/TEST/cfs-python-matlab/sgp_rot/filter_matrix.mtx');
      
      id = speye(this.n_elems);
      this.dAF = (id-this.F)'*(id-this.F); % bei Uebergabe noch halbieren ...
      this.dAFii = diag(this.dAF);
    end
    
    function F = get_filter(this, rad)
      % check if filter was already computed
      if this.F ~=-1
        % matlab engine does not support sparse matrices as return type to python
        F = full(this.F);
      else
        F = calc_filter(this.nx, this.ny, this.nz, rad);
      end
    end
    
    
    %% returns rotated stiffness tensor/matrix for all combinations of dens and alpha
    % alpha: assume 2d - rotation about one axis only
    % dens: pseudo density
    function [C] = eval_rot_tens(this, dens, alpha)
      assert(all(size(dens) == size(alpha)));
      n = size(dens,1);
      % 2D: 3 x 3: 3D: 6 x 6
      r = (this.dim-1)*3;
      C = zeros(n,n,r,r);
      
      for i = 1:n
        fact = dens(i)^this.p;
        for j = 1:n
          % \rho^p * Q*D*Q^T
          C(i,j,:,:) = this.calc_rot_tensor(fact,alpha(j));
        end
      end
    end
    
    %% returns material tensor for given pair of dens and alpha
    function [C] = eval_tensor(this, dens, alpha)
      assert(all(size(dens) == size(alpha)));
      % output is a cell of cells
      C = arrayfun(@(x,y) calc_H(this.H0,x,y,this.p),dens,alpha,'un',false);
      % convert to array of size: 3 x 3 x nelems
      %fprintf(1,'size Q: %d %d; size alpha: %d %d', size(Q), size(alpha))
      C = reshape(cell2mat(C)',3,3,size(C,1));
      
%       C = zeros(r,r,n);
%       for i = 1:n
%         C(:,:,i) = calc_H(this.H0,dens(i),alpha(i),this.p);
%       end
    end
    
   % function [C] = eval_inv_stiffness_tensor(this, dens, alpha)
    function [C_inv] = eval_inv_stiffness_tensor(this,C)
      % in 2d: C is a 3 x 3 x #rho*# phi 
      assert(ndims(C) == 3);
      n = size(C,1);
      C_inv = zeros(size(C));
      for i = 1:n
          C_inv(i,:,:) = inv(squeeze(C(i,:,:)));
      end
    end
    
    function [x_new,f_min] = solve(this,x,dfH,lambda)
      % x is a long vector: [\rho_1, \rho_2, \rho_2, ...., \phi_1, \phi_2, \phi_3, ...]
      x = reshape(x,this.n_elems,2);
      
      %fprintf(1,"MATLAB solve(): size of x=%d %d; size of dfH= %d %d %d\n", size(x), size(dfH));
      f_min = 0;
      
      % permute H and dfH to conform with input of mtimesx:
      % orig size of H: n_elems x 3 x 3  -> input of mtimesx: 3 x 3 x n_elems
      % dA = df(H^(k))/dH
      dA = permute(dfH, [2,3,1]);
      
      Hbar = this.eval_tensor(x(:,1),x(:,2));
      
      if this.rad_filt > 0. && this.filt_sens
        % sensitivity filtering
        trE = zeros(this.n_elems,1);
        for e = 1:this.n_elems
          trE(e) = trace(Hbar(:,:,e));
        end
        % reorganize derivative ...
        dff = reshape(dA,3,3,this.nx,this.ny);
        trE1 = reshape(trE,this.nx,this.ny)/trace(this.H0);
        dff = filterE(this.nx,this.ny,this.rad_filt,trE1,dff);
        dA = reshape(dff,3,3,this.n_elems);
      end
      
      %% assume L=0 
      % dLA = Q*(H^(k) - L)*Q'
      dLA = Hbar;
      
      %fprintf(1," size dA: %d %d %d; size dLA: %d %d %d\n", size(dA), size(dLA));
      dAL = mtimesx(dA,dLA);
      
      % -Q*(H^(k) - L)*Q' * dfH * Q*(H^(k) - L)*Q'
      BBglobA = -mtimesx(dLA,dAL);
      
      %fprintf("MATLAB size this.dAF=%d %d; size x=%d %d; size this.dAFii=%d %d\n",size(this.dAF),size(x),size(this.dAFii));
      % density values: x(:,1)
      % angle values: x(:,2)
      dAFiR1 = this.dAF*x(:,1)-x(:,1).*this.dAFii;
      dAFiR2 = this.dAF*x(:,2)-x(:,2).*this.dAFii;
      
      % number of points on finest grid in one dimension
      npoints_all = (this.samples)^this.levels+1;
      
      m = this.n_elems;
      jstart = ones(m,1); kstart = ones(m,1);
      
      % number of samples on current level
      m_grid_bisect = this.samples;
      grid_bisect = 0:m_grid_bisect;
      
      %fid = fopen("fast_comp.txt", "w");
      
      for q=1:this.levels
        ddx = this.samples^(this.levels-q);
        % indices for angle values on current level
        jjall = bsxfun(@plus,jstart,ddx*grid_bisect);
        % indices for density values on current level
        kkall = bsxfun(@plus,kstart,ddx*grid_bisect);
        
        Ftmpall = zeros(this.n_elems,(m_grid_bisect+1)^2);
        for i=1:m
          idx_d = jjall(i,:)';
          idx_a = kkall(i,:)';
          % generate pairs of indices for dens and angle
          % see https://stackoverflow.com/questions/7446946/how-to-generate-all-pairs-from-two-vectors-in-matlab-using-vectorised-code
          [aa dd] = meshgrid(idx_a,idx_d);
          %pairs_i = [dd(:) aa(:)];
          % vary all densities, the angles
          idx = sub2ind([npoints_all npoints_all], dd(:), aa(:));
          HALLi = this.H_all_i(idx,:);
          if this.dim == 3
            % store only non-zero values
            % we need factor 2 for off-diagonal entries
            % as Frobenius inner product is calculate with half of the matrix (symmetry)
            % assume orthotropic tensor
            BBBi = [BBglobA(1,1,i); BBglobA(2,2,i); BBglobA(3,3,i); ...
              BBglobA(4,4,i); BBglobA(5,5,i); BBglobA(6,6,i); ...
              2*BBglobA(2,3,i); 2*BBglobA(1,3,i); 2*BBglobA(1,2,i); ...
              lambda; 0 ];
          else % dim = 2
            BBBi = [BBglobA(1,1,i); BBglobA(2,2,i); BBglobA(3,3,i); ...
              2*BBglobA(2,3,i); 2*BBglobA(1,3,i); 2*BBglobA(1,2,i); ...
              lambda; 0 ];
          end
          
          XALL = this.design_all(idx,:);
          dAFiRi = [dAFiR1(i); dAFiR2(i); 0.5*this.dAFii(i) ];
          Ftmpall(i,:) = HALLi*BBBi;
          if this.filt_sens ~= 1
            Ftmpall(i,:) = Ftmpall(i,:) + this.p_filt*[XALL*dAFiRi]';
          end
          
        end %element
        
        [f_min,idxmin] = min(Ftmpall,[],2);
        
        [jmin,kmin]=ind2sub([m_grid_bisect+1,m_grid_bisect+1],idxmin);
        
        % compute indices of found min values on finest grid
        jval_tmp = jstart + ddx*grid_bisect(jmin)';
        kval_tmp = kstart + ddx*grid_bisect(kmin)';
        
        % for debugging
%         for i = 1:m
%           sol = [this.x_all(jval_tmp) this.a_all(kval_tmp)];
%           fprintf(1,"level= %d element=%d obj_min=%f rho=%f theta=%f\n",q,i,f_min(i),sol(i,1),sol(i,2));
%           fprintf(fid,"level= %d element=%d obj_min=%f rho=%f theta=%f\n",q,i,f_min(i),sol(i,1),sol(i,2));
%         end
%         fprintf(1,"\n");
%         fprintf(fid,"\n");

        %% set lower interval bound for next level
        % special cases: we're already at lower or upper design value bound
        i_jstart = jmin-1;   
        i_jstart(i_jstart < 1) = 1;
        i_jstart(i_jstart > m_grid_bisect-1) = m_grid_bisect -1;
        
        % same for k
        i_kstart = kmin-1;   % Achtung: Index, wie er aus min funktion (nach ind2sub) rauskommt
        i_kstart(i_kstart < 1) = 1;
        i_kstart(i_kstart > m_grid_bisect-1) = m_grid_bisect -1;
        
        % valid for next level
        jstart = jstart + ddx*grid_bisect(i_jstart)';
        kstart = kstart + ddx*grid_bisect(i_kstart)';
        
        % e.g. jmin = 2
        % next level: left = 1 and right = 3
        % -> factor 2
        if q == 1
          m_grid_bisect = 2*this.samples;
          grid_bisect = 0:m_grid_bisect;
        end
        
      end % level
      
      x_new = [this.x_all(jval_tmp)' this.a_all(kval_tmp)'];
      % debugging
%       tmp = this.x_all(jval_tmp);
%       fmt = ['xnew : [', repmat('%g, ', 1, numel(tmp)-1), '%g]\n'];
%       fprintf(1, fmt, tmp);
    end % function
    
  end % class methods
  
end % class

