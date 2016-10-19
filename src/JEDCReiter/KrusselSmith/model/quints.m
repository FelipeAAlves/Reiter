function [qu,shares] = quints(x0,Dmat)
[x,p] = ptd_merge(x0,Dmat);

rqu = linspace(0,1,6)';
ranges = [[rqu(1:5);0.9;0.95;0.99] [rqu(2:6);0.95;0.99;1]];
for i=1:8
  [x2,p2] = ptd_percrange(x,p,ranges(i,:));
  p2 = p2;
  qu(i) = dot(x2,p2);
end

X = sum(qu(1:5));  % sum up over the quintiles
shares = qu./X;

