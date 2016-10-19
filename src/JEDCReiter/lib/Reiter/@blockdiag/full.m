function y = full(x)
  nn = size(x);
  y = zeros(nn);
  % y = initsize(zeros(nn,x.data{:}));
  i1 = 0;
  i2 = 0;
  for i=1:length(x.data)
    [nr,nc] = size(x.data{i});
    y(i1+1:i1+nr,i2+1:i2+nc) = full(x.data{i});
    i1 = i1+nr;
    i2 = i2+nc;
  end