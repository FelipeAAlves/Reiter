%%% Description
%       Program to check existence and uniqueness in linear RE models
%%% INPUT:
%   g0,g1,psi,pi,div: as defined in Sims' paper
%
% OUTPUT:
%   eu: stable solution exists iff eu(1)==1
%       solution is unique iff eu(2)==1
%   eigvals: generalized eigenvalues
%
% Michael Reiter, November2005
function [eu,eigvals,G1,impact] = checkeu(g0,g1,psi,pi,div)

  [Lambda Omega q z]=qz(g0,g1);  %5th argument not needed!
  clear g0 g1;
  % Note that q'*Lambda*z' = g0 etc.,
  %  NOT  q*Lambda*z = g0 (as the Matlab help says!)

  realsmall = 1e-10; %much stricter than Sims
  dLambda = abs(diag(Lambda));
  dOmega = abs(diag(Omega));
  if(any(dLambda<realsmall & dOmega<realsmall))
    % taken from Sims:
    disp('Coincident zeros.  Indeterminacy and/or nonexistence.')
    eu=[-2;-2];
    return;
  end

  warning off;
  ev = abs(diag(Omega)./diag(Lambda));
  warning on;
  evavs = flipud(sort(ev));

  if(str2num(version('-release'))<14) % use Sims' file
    warning('I use Sims ordering algorithm, because ordqz not available');
    [Lambda Omega q z]=qzdiv(div,Lambda,Omega,q,z);
  else
    select = ev<div;
    [Lambda Omega q z]=ordqz(Lambda,Omega,q,z,select);
  end

  dLambda = abs(diag(Lambda));
  dOmega = abs(diag(Omega));
  dLambda = max(dLambda,realsmall); %to avoid dividing by 0;
  n = size(Lambda,1);
  n_unstable = sum(dLambda<=realsmall | dOmega./dLambda>div);
  n_stable = n-n_unstable;

  iStable = 1:n-n_unstable;
  iUnstable = n-n_unstable+1:n;

  q1=q(iStable,:);
  q2=q(iUnstable,:);
  q1pi=q1*pi;
  q2pi=q2*pi;
  q2psi=q2*psi;

  iExist = rank(q2pi) == rank([q2psi q2pi]);
  clear q2psi;
  iUnique = rank(q*pi) == rank(q2pi);
  clear q;

  Phi = q1pi/q2pi;  %Sims, equ. (42)
  clear q1pi q2pi;

  eu = [iExist; iUnique];
  warning off;
  ev = diag(Omega)./diag(Lambda);
  % sort according to absolute value:
  ev = sortrows([ev abs(ev)],2);
  eigvals = ev(:,1);
  warning on;


  z1 = z(:,iStable);
  %  z2 = z(:,iUnstable); % NOT NEEDED
  L11 = Lambda(iStable,iStable);
  L12 = Lambda(iStable,iUnstable);
  L22 = Lambda(iUnstable,iUnstable);
  clear Lambda;
  O11 = Omega(iStable,iStable);
  O12 = Omega(iStable,iUnstable);
  O22 = Omega(iUnstable,iUnstable);
  clear Omega;

  %Sims, equ. (45)
  L11inv = inv(L11);
  clear L11;
  aux = [O11 (O12-Phi*O22)]*z';
  clear O11 O12 O22;
  aux2 = z1*L11inv;
  clear z1;
  G1 = real(aux2*aux);
  clear aux aux2;
  aux = [L11inv -L11inv*(L12-Phi*L22);zeros(n_unstable,n_stable) eye(n_unstable)];
  clear L11inv L12 L22;
  H = z*aux;
  clear z;
  impact = real(H*[q1-Phi*q2;zeros(n_unstable,size(psi,1))]*psi);
