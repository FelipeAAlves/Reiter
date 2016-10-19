function [A] = qr_house(A);

% input:   A - (n x m) matrix 
% output:  A - matrix containing the necessary NOTErmation of the Householder
%              vectors v in the lower triangle and R in the upper triangle

% Author     : Stefan H�eber
% Date       : May, 8, 2003
% Institution: University of Stuttgart,
%              Institut for Applied Analysis and Numerical Mathematics,
%              High Performance Scientific Computing
% Version    : 1.0

[n,m] = size(A);
for k = 1:min(n-1,m)
	v(k:n,1) = householder(A(k:n,k));
	A(k:n,k:m) = householder_mult(A(k:n,k:m),v(k:n,1));
	A(k+1:n,k) = v(k+1:n,1);
end