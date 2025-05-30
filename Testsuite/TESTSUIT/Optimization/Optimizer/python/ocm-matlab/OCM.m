%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef OCM < handle
  properties (Access = private)
    xl = -1e6;
    xu = 99999;
    move = 0.2;
    damp = 0.5;
  end
  
  methods
    %% set ocm parameters via constructor
    function this = OCM(xl, xu, move, damp)
      this.xl = xl;
      this.xu = xu;
      this.move = move;
      this.damp;
      fprintf("MATLAB: set  xl: %2.6f, xu: %2.6f, move: %2.6f, damp: %2.6f \n", this.xl, this.xu, this.move, this.damp) 
    end
    
    %% does one ocm update without checking
    % x: current design  
    % lambda: volume constraint multiplier
    % returns xnew: updated design
    function [xnew] = oc_update(this, x, dc, dv, lambda)
      %fprintf("dc: %2.6f, dv: %2.6f, lambda: %2.6f",dc,dv,lambda)
      xnew = max(this.xl,max(x-this.move, min(this.xu, min(x+this.move,x.*(-dc./dv/lambda).^this.damp) ) ) );
    end
  end
end

