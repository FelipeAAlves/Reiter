function y = apply2cell(fcn,x,x2)
y = cell(size(x));
for i=1:prod(size(x));
  if(nargin<=2)
    y{i} = fcn(x{i});
  else
    y{i} = fcn(x{i},x2{i});
  end
end