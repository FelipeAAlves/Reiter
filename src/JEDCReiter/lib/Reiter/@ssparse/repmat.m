function r = repmat(m,nr,nc)
  c = cell(nr,1);
  for i=1:nr
    c{i} = m;
  end
  m2 = vertcat(c{:});
  c = cell(nc,1);
  for i=1:nc
    c{i} = m2;
  end
  r = horzcat(c{:});

  
  