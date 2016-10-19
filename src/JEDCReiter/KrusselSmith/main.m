%
% File  : main.m
% Auhtor: Felipe Alves
%
% Based on
% 1. "Approximate and Almost-Exact Aggregation in Dynamic Stochastic Heterogeneous-Agent Models"
%     Michael Reiter, IHS, December 2009
%
% 2. "A Toolbox for Solving and Estimating Heterogeneous Agent Macro Models"
%     Thomas Winberry.
%
%    Program to solve Krussel-Smith Model by Reiter & Winebrry method

clc;
clearvars
clear global

%%%  OPTION  %%%
plt_figure = 1;

AggregateUncertainty = 1;

Reiter   = 1;   % does Reiter  Linearization Step
Winberry = 0;   % does Wiberry Steady States Step
% .....................................................................................

% set the Matlab path
setpath;

% set the baseline parameters
initModelParam;

% SteadyState
solveStSt;

%% Linearize
if AggregateUncertainty
    if Reiter
    linearizeHistogram;
    end
    if Winberry
        linearizePolynomial;
    end
%     IRFFigures;
    SimulFigures;
end
