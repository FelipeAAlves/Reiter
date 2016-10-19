function [ser,cov_matrix] = simul_sims(G1, Impact, T_or_Shocks, Sigma, H)
  global file_suffix;

  nx = size(G1,1);
  if (nargin<5 | isempty(H))
    H = 1;
    nh = nx;
  else
    nh = size(H,1);
  end


  nz = size(Impact,2);
  C = chol(Sigma)';

  if (prod(size(T_or_Shocks))==1)  % is a scalar
    T = T_or_Shocks;
    shocks = (C*randn(nz,T));
  else
    shocks = C*T_or_Shocks;
    T = size(shocks,2);
  end



    ser = zeros(nh,T+1);
    x   = zeros(nx,1);
    for t = 2 : T+1
        x = G1 * x + Impact * shocks(:,t-1);
        ser(:,t) = H * x;
    end
    ser = ser';   % each column a variable

  if (nargout>1)
    % covmatrix of output series satisfies
    %   cov_matrix = G1*cov_matrix*G1' + Impact*Sigma;
    cov_matrix = H*lyapunov_symm(G1,Impact*Sigma*Impact')*H';
  end

  return;

end
