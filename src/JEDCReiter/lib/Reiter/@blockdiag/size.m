function [n,varargout] = size(x,varargin)
  nr = 0;
  nc = 0;
  for j=1:length(x.data)
    nr = nr + size(x.data{j},1);
    nc = nc + size(x.data{j},2);
  end
  nn = [nr nc];

  if(nargout==0)
    n = nn;
  else
    varargout = cell(nargout-1,1);
    [n,varargout{:}] = handle_size(nn,varargin{:});
  end
