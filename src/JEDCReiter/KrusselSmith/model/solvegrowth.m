% Matlab file for_ paper 
% "Approximate and Almost-Exact Aggregation in Dynamic Stochastic Heterogeneous-Agent Models"
% Michael Reiter, IHS, December 2009
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% Program to solve growth model with almost-exact aggregation
%
% global structure of parameter values:
global MP;
neps = 1;  % permanent productivity states; neps=1 means only i.i.d shocks
           % neps=31 is the big calibration in the paper
sigtau = 0.01;  % stdev tax shock
freq = 4;  % frequency  (4 is quarterly etc.)
initMP(neps,freq,sigtau);

% COMPUTE STEADY STATE:
[Par,Kequ,D] = hetsolve;
% COMPUTE DSGE
[Par,Kequ,D,G1,impact,abseig,Env] = hetsolve(Kequ,Par);
  

% index of static and backward looking variables:
iS = Env.iVarBWS;
% number of exogenous shocks:
nz = MP.nz;
n = length(iS);
% variable to be predicted: (aggregate capital)
Hpred = full(eyeii(n,Env.iVarStatic));  % aggr. capital is a
                                        % variable, with index Env.iVarStatic



% do almost (machine-precision) exact aggregation:
[G1_MPA,Impact_MPA,A_Dis,B_Dis,ADec_Dis,BDec_Dis,G,dummy_notused,H_MPA] = do_mpa_orig(Env,MP.freq*250);
% accuracy check 1: compare small model to induced big model:
ii = 1:size(H_MPA,2);
A_MPA = G1_MPA(ii,ii);
B_MPA = Impact_MPA(ii,:);
C_MPA = coord_h(H_MPA,Hpred);


% ACCURACY CHECKS
% check against disaggregate solution of reduced model (suffix _Dis)
%   output: - difference in impulse response
%           - max. error  in transition
%           - RMSE of infinite-horizon forecast
% length of IR etc. in accuracy check
% Tmax = MP.freq*250*ones(2,1);  
% [errAggr.diffIR_Dis,errAggr.trans_Dis,errAggr.RMSE_Dis] ...
%     = evalmod(A_Dis,B_Dis,Hpred',{A_MPA},{B_MPA},{C_MPA'},{H_MPA'},Tmax,MP.SigmaEps,MP.beta);
% if(MP.neps==1)
%   % also check against exact solution (suffix _Ex)
%   A = G1(iS,iS);
%   B = impact(iS,:);
%   [errAggr.diffIR_Ex,errAggr.trans_Ex,errAggr.RMSE_Ex] ...
%       = evalmod(A,B,Hpred',{A_MPA},{B_MPA},{C_MPA'},{H_MPA'},Tmax,MP.SigmaEps,MP.beta);
% end

