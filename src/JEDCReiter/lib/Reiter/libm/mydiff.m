% serves as default for_ some classes:
function d = mydiff(m)
  n = size(m,1);
  d = m(2:n,:) - m(1:n-1,:);
