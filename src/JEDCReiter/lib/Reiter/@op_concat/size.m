function [n,varargout] = size(S,varargin)

  na = size(S.a);
  nb = size(S.b);
  nn = [na(1) nb(2)];

  if(nargout==0)
    n = nn;
  else
    varargout = cell(nargout-1,1);
    [n,varargout{:}] = handle_size(nn,varargin{:});
  end

  