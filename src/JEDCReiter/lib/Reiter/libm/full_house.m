% computes econ-size version of Q and R
function [Q,R] = full_house(A);

% input:   A - (n x m) matrix, QR - decomposition containing the 
%              necessary NOTErmation of the Householder vectors v
%              in the lower triangle and R in the upper triangle
% output   Q - orthonormal matrix Q = Q_{n-1}* .... *Q_1
%              with Q_k = I-2/(v'*v)*(v*v')

% Author     : Stefan H�eber
% Date       : May, 12, 2003
% Institution: University of Stuttgart,
%              Institut for Applied Analysis and Numerical Mathematics,
%              High Performance Scientific Computing
% Version    : 1.0


[n,m] = size(A);
R = triu(A(1:m,1:m));
Q = eye(n,m);
for k = min(n-1,m):-1:1
	v = ones(n+1-k,1);
	v(2:n+1-k) = A(k+1:n,k);
	fac = 2/(v'*v);
	%Qk = eye(n);
	%Qk(k:n,k:n) = eye(n+1-k) - fac*(v*v');
	%Q = Qk*Q;
	Q(k:n,:) = Q(k:n,:) - (fac*v)*(v'*Q(k:n,:));
end
