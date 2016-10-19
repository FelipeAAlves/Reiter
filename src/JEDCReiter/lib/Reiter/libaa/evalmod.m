% inputs A,B,C: "exact" solution
% A2,B2,C2:  cell array of approximate solutions
%   note: C and C2 give the variables for which accuracy is computed
% H: for each approximation, as in paper (row-form)
% T(1): length of IR; 
% T(2): length max-norm calculation
function [diffIR, maxerr_trans, stderr_inf,maxerr_long] = evalmod(A,B,C,A2,B2,C2,H,T,Sigma,discount)

global MP;  %Added 2/17/12 by AGM

scaleshock = 1;
irEX = ir_sims(A,B,T(1),C,scaleshock);

nm = length(A2);
nc = size(C,1);
nz = size(B,2);
diffIR = NaN(nm,nc,nz);
maxerr_trans = NaN(nm,1);
maxerr_long = NaN(nm,1);
stderr_inf = NaN(nm,nc);

for imom=1:nm
  irDSF = ir_sims(A2{imom},B2{imom},T(1),C2{imom},scaleshock);
  

  for iz = 1:nz
    for ic = 1:nc
      diffIR(imom,ic,iz) = distance(irEX{iz}(:,ic),irDSF{iz}(:,ic),1);
      % diffIR(imom,ic,iz) = norm(irEX{iz}(:,ic)-irDSF{iz}(:,ic))/scaleshock;
    end
  end

  if(nargout>1)
    d1 = C*A;
    d2 = A2{imom}*H{imom};

    for i=1:T(2)
      res = d1 - C2{imom}*d2;
      ns(i) = norm(res);  %AGM (2/25) this is what takes for ever!
      d1 = d1*A;
      d2 = A2{imom}*d2;

    end

    % choose maximum of norms:
    disc = discount.^(0:1:length(ns)-1);
    maxerr_trans(imom) = max(disc.*ns);
    % chose sum of square of norms:
    % maxerr_trans(imom) = sum(ns.^2);
  end
end

if(nargout>2)
  if(length(A)<=2000)
    [S1,S2,S12,covmodel,covappr,maxerror] = doubling_joint(A,B,A2,B2,Sigma,30,C,C2,H);
    disc = discount.^(2.^(0:1:size(maxerror,1)-1)-1)';
    maxerr_long(1:size(maxerror,2)) = max(repmat(disc,1,size(maxerror,2)).*maxerror)';
    for i=1:length(S2)
      if(any(covappr{i}<0))
        warning('NaN in RMSE computation');
        ser = simul_sims(A,B,shocks(nz),Sigma,C);
        ser2 = simul_sims(A2{i},B2{i},shocks(nz),Sigma,C2{i});
        var1 = var(ser);
        vardiff = var(ser-ser2);
        stderr_inf(i,:) = sqrt(vardiff./var1);
      else
        stderr_inf(i,:) = sqrt(diag(covappr{i})./ diag(covmodel));
      end
    end
  else
    ser = simul_sims(A,B,shocks(nz),Sigma,C);
    for i=1:length(A2)
      ser2 = simul_sims(A2{i},B2{i},shocks(nz),Sigma,C2{i});
      var1 = var(ser);
      vardiff = var(ser-ser2);
      stderr_inf(i,:) = sqrt(vardiff./var1);
    end
  end
end

function ser = shocks(nz)
  s = rng;  % save old state  [updated to rng 3/7/12 by AGM]
  rng(2937414);  %to draw always same shocks, for compatibility
  ser = randn(nz,100000);            % 2 -> MP.nz.  changed 2/17/12 by AGM
  rng(s);  % restore old state [updated to rng 3/7/12 by AGM]