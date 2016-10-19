function s = sum(A)
  s = zeros(1,size(A,2));
  nc = 0;
  for i=1:length(A.data)
    ni = size(A.data{i},2);
    s(nc+1:nc+ni) = sum(A.data{i});
    nc = nc + ni;
  end