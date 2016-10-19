function xt = ctranspose(x)
  xt = x;
  for i=1:length(x.data)
    xt.data{i} = xt.data{i}';
  end