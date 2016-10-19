function [n,varargout] = size(S,varargin)
  nn = size(S.a);
    
  if(nargout==0)
    n = nn;
  else
    varargout = cell(nargout-1,1);
    [n,varargout{:}] = handle_size(nn,varargin{:});
  end

  