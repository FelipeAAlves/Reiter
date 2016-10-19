function is = dissparse(x)
  % assert(~issparse(x.d));  %WHY THIS???
  is = issparse(x.d) | isa(x.d,'ssparse');
