% global structure of parameter values:
global MP;
neps = 1;  % permanent productivity states
sigtau = 0.01;  % stdev tax shock
freq = 4;  % frequency
initMP(neps,freq,sigtau);

h1short = make_h(1);
HSumK = [h1short;zeros(n-length(h1short),1)];

% do  state space reduction:
if(MP.neps==1)  % compute results for state space reduction:

  cholSigmaA = doubling_chol(A,B*chol(getsigma)',30);
  [ObservGram,PCA,BalRed] = reduc3(A,B,cholSigmaA,getsigma, HSumK,0);%

  scaleshock = 0.01;
  nMoms = [1:10 12 15 20:5:40 50:10:80 100];
  nOrders = nMoms;
  %take HSumK, not Hbasic, because hpoly:
  CEAstr = gramianspace(A',HSumK,MP.freq*250);
  CEA = CEAstr.u;
  [mom_used,momtypes, diffIR, maxerr_trans, stderr_inf,maxerr_long] = ...
      ssapprox(A,B,getsigma,Tmax,scaleshock,HSumK,cholSigmaA,...
               Env.iExog,...
               range(PCA),...%               range(MPA3),...
               CEA,...
               range(BalRed),nMoms,nOrders,@hpoly,MP.beta);
  save(filename('resapprox'),'mom_used','momtypes','diffIR', ...
       'maxerr_trans','stderr_inf','maxerr_long');

  if(strcmp(momtypes{4},'ApplSpec'))
    momtypes{4} = 'Moments';
  end
  graphreduc(stderr_inf, diffIR, maxerr_trans, mom_used,momtypes, {'Z','\tau'},[filename('reduc') '.eps']);
end