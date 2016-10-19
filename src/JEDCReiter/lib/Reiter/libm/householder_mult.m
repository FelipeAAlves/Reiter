function [A] = householder_mult(A,v);

% input:   A - matrix
%          v - Householder vector 
% output:  A - transformed matrix 
%              A = A - \(\frac{\text{2}}{\text{v'v}}(\text{vv'})\)A 

% Author     : Stefan Hüeber
% Date       : May, 8, 2003
% Institution: University of Stuttgart,
%              Institut for Applied Analysis and Numerical Mathematics,
%              High Performance Scientific Computing
% Version    : 1.0

vv=v'*v;
%if(vv==0)
%    A=[];
%    return;
%end
beta = -2/(vv);
w = v'*A;          % w is a line vector
A = A + beta*v*w;
