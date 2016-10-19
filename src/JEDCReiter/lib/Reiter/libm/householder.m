function [v] = householder(a);

% input:   a - vector with a \neq 0
% output:  v - Householder vector of x such that v_1 = 1

% Author     : Stefan Hüeber
% Date       : May, 8, 2003
% Institution: University of Stuttgart,
%              Institut for Applied Analysis and Numerical Mathematics,
%              High Performance Scientific Computing
% Version    : 1.0

n = length(a);
v = a;
if (a(1) >= 0) beta = a(1) + norm(a);
else beta = a(1) - norm(a);
end
if(beta==0)
  v=[];
else
  v(2:n) = 1/beta * v(2:n);
  if(~all(isfinite(v)))
    v=[];
  else
    v(1) = 1;
  end
end
