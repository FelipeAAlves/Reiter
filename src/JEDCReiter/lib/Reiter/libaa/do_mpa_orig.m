% Function to solve Linear Rational Expectations Model by 
%   state and policy reduction,
%   with accuracy up to machine precision  (almost-exact, MP-exact)
% Outputs:
%   AReduc: state and policy transition matrix of reduced model
%             (conforms to G1 in Sims' package)
%   BReduc: impact matrix of reduced model
%             (conforms to Impact in Sims' package)
%   AState: implied state transition matrix for original (big) model
%   BState: implied impact matrix for big states
%   ADec: implied decision matrix for original (big) model (this is D_X is
%   Reiter's notation)
%   BDec: implied impact matrix for big decisions (this is D_E in Reiter's
%   notation)
%   G:  subspace in which decision vector lives
%   CEX: structure with results of conditional expectations approach
%   H: state aggregation matrix
%   C: space of state variables that enter Euler equation directly
function [AReduc,BReduc,AState,BState,ADec,BDec,G,CEX,H,C] = do_mpa_orig(Env,horizon)

    disp('entering do_mpa');

  splitjac;  %generates matrices T, D etc. from filename

  nStates = size(T,1);

  iEs0 = find(sum(abs(Es0),1)>1e-14);  %1 added to sum to
                                       %accomodate 1d forward
                                       %looking var.  AGM 2/18/12
  iEs1 = find(sum(abs(Es1),1)>1e-14);  % as above
  
  %assert(length(iEs0)+length(iEs1)<=100);  %write decisions as functions of few state variables, not of distribution!!
  iEs = union(iEs0,iEs1);

  Ess = full([Es0(:,iEs);Es1(:,iEs)])';
  oe = orth(Ess);
  nStates = length(iVBW);
  nch = size(oe,2);
  C = zeros(nStates,nch);
  C(iEs,:) = oe;
  if(length(T)<1500)
    CG = gramian(full(T)',C,0);
  else
    CG = gramian(T',C,horizon);  %see page 20
  end
  [CEX.u,CEX.s] = us(CG);  %SVD of Q = CG
  ii = find(CEX.s>CEX.s(1)*1e-16);  
  H = CEX.u(:,ii)';  %create H corresponding to non-zero singular values

  nMom = size(H,1);
  % H = eye(size(T));
  hatT = full(H*(T*H'));  %see eq 59
  nDSF = length(hatT);
  maxabseig = max(abs(eig(hatT)));
  if(maxabseig>=1)
    fprintf(1,'maximal eigenvalue of hatT = %e\n',maxabseig);
    error('unstable reduced model');
  end

  
  Em0 = Es0*H';  %eq 60
  Em1 = Es1*H';

  hatA = hatT;  %page 22 bullet 1
  Ed10 = Ed1\Ed0;
  nIter = 4;
  for i=1:nIter
    eea = (Em1 + Em0*hatA)*hatA;    
    Dm1 = sylvschur(Ed1,Ed0,hatA,-eea);   %page 22 bullet 2 (eq 75)
    opts.disp=0;
    % [vv,dd] = eigs(Ed10,10,'LM',opts);
    vv = [];
    G = orth([Dm1 real(vv)]);  %page 22 bullet 3
    N = G';

    nDec = size(G,2);
    nTot = nMom+nDec;

    % Sims: -Ss0 - Sd0 = Ss1 + Se0
    % here:  I - (-Ss0\Sd0)  = (-Ss0\Ss1) + (-Ss0\Se0)
    % c = zeros(size(g0,1),1);
    if(0)
      g0 = [eye(nDSF) -H*(D*G);
            -N*Em0   -N*Ed0*G];
      g1 = [hatT zeros(nDSF,nDec);
            N*Em1     N*Ed1*G];
      Psi = [H*Se;N*Ee0];
    else

      %prep for page 22 bullet 4
      HF = H*Se0;
      nz = size(HF,2);
      sigT = doubling(hatT,quadform(HF,eye(nz)),30);
      HD = Dm1'*G;
      sigHD = sigT*HD;
      sigD = HD'*sigHD;
      sigD = (sigD+sigD')/2;
      
      sigDreg = sigD + 1e-10*max(abs(eig(sigD)))*eye(size(sigD));
      
      VV = G'*Ed10*(Dm1*sigT*Dm1')*G;
      hatF = VV/sigDreg;

      if(1)
        g0 = [eye(nDSF) -H*(D*G);
              -N*(Ed1\Em0)   -N*Ed10*G];
      else
        g0 = [eye(nDSF) -H*(D*G);
              -N*(Ed1\Em0)   -hatF];
      end
      g1 = [hatT zeros(nDSF,nDec);
            N*(Ed1\Em1)     N*G];
      Psi = [H*Se;N*(Ed1\Ee0)];
    end
    Pi = [zeros(nDSF,nDec);eye(nDec)];
    div = 1+1e-10;

    % g0*y(t) = g1*y(t-1) + Psi*z
    % LONG RUN: dy/dz = inv(g0-g1)*Psi;
    compstat = (g0-g1)\Psi;

    %page 22 bullet 4 (gensys)
    [eu,eigvals,AReduc,BReduc] = checkeu(g0,g1,Psi,Pi,div);
    
    
    abseig = flipud(sort(abs(eigvals)));
    if(any(eu==0))
      fprintf(1,'eu = %d,%d\n',eu(1),eu(2));
    end
    if(i>1)
      ii=1:nDSF;
      Creduc = H'\C;
      fprintf(1,'distance between G1: %e; w.r.t. C: %e\n', ...
              distance(AReduc(ii,ii),AReducLast(ii,ii),1),...
              distance_h(AReduc(ii,ii),AReducLast(ii,ii),Creduc',1000));
    end
    
    
    if(i==nIter)
      %ADec = G*AReduc(nMom+1:nTot,1:nMom)*H;
      ADec = op_concat(G,AReduc(nMom+1:nTot,1:nMom)*H);
      BDec = G*BReduc(nMom+1:nTot,:);
      % BState = op_sum(Se,op_concat(D,BDec));
      BState = Se + D*BDec;  %In Reiter's notation this is F + D * D_E
      AState = op_sum(T,op_concat(D,ADec)); %Reiter notation: T + D * D_X

      % perhaps iterate on that one:
      % [hatA,hatB,AState,BState] = updatesol(Env,H,AState,BState);
      % OR:
      % [AReduc,BReduc,AState,BState] = updatesol(Env,H,AState,BState);

      %G1 = blockdiag([AState;ADec],zeros(0,length(MP2.iVarDec)));
      %B = [BState;BDec];
    end
    
    %change hat A and go back to page 22 bullet 2
    hatA = AReduc(1:nDSF,1:nDSF);
    AReducLast = AReduc;
  end

  H = H';  % write in column form