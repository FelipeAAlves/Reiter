function [dist,err] = distance_h(A1,A2,C,Tmax)
d1 = C*A1;
d2 = C*A2;
for i=1:Tmax
  err(i) = norm(d1-d2);
  d1 = d1*A1;
  d2 = d2*A2;
end
dist = max(err);

