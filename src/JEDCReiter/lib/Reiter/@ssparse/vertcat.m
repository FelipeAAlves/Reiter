% Matlab file to implement class "ssparse",
% Michael Reiter, Universitat Pompeu Fabra, April 2007
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
% 
function xOut=vertcat(varargin)
  nv = length(varargin);

  iR=[];iC=[];V=[];nr=0;nc=0;
  for i=1:nv
    if(~isempty(varargin{i}))
      n = size(varargin{i});
      if(nc==0)
	nc=n(2);
      else
	if(nc~=n(2))
	  error('incompatible matrices in vertcat');
	end
      end
      [ir,ic,v] = find(varargin{i});
      ir = ir+nr;
      iR = [iR;ir];
      iC = [iC;ic];
      V = [V;v];
      nr = nr+n(1);
    end
  end
  xOut = ssparse(iR,iC,V,nr,nc);