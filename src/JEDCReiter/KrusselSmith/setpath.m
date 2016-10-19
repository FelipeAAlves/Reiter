

%%% Setting paths:

% set here the directory into which appraggr.tar.gz is installed:
p_ = pwd;

% libraries:
% addpath([p_ '/lib/MirandaFackler']);                                      % PC STERN falves: already on path!!
% addpath('/Users/Felipe/Documents/MATLAB/Toolbox/COMPECON/CEtools')        % Mac falves: have to add

addpath([p_ '/../lib/Reiter']);
% addpath([p_ '/../lib/Reiter/Reiter_original']);

addpath([p_ '/../lib/Reiter/libaa']);
addpath([p_ '/../lib/Reiter/libm']);
addpath([p_ '/../lib/Reiter/rrqr']);

% Model
addpath([p_ '/model']);

% to allow loading of gramian objects
% a = gramian(0.5,0.5,100,1e-14); clear a;
