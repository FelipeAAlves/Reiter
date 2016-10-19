% falves notes
% -------------
%
% X      :  point at which to take derivatives.
% iDeriv :  index of which variables in X are to be differentiated with respect to.
%           other parts of X are assumed constant.
%
% nmax tells adjacob how many derivatives to take at a time.  I believe this has to do with
% memory management.
%

clc; clear all;
X = 0.5*ones(2,1);
iDeriv = 1:2;
nmax = 2;

J1 = adjacob(@test_equ1,X,iDeriv,nmax);
J2 = adjacob(@test_equ2,X,iDeriv,nmax);
