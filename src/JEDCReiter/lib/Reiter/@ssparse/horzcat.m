% Matlab file to implement class "ssparse",
% Michael Reiter, Universitat Pompeu Fabra, April 2007
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
% 
function xOut=horzcat(varargin)
  nv = length(varargin);

  iR=[];iC=[];V=[];nr=0;nc=0;
  for i=1:nv
    if(~isempty(varargin{i}))
      n = size(varargin{i});
      if(nr==0)
	nr=n(1);
      else
	if(nr~=n(1))
	  error('incompatible matrices in horzcat');
	end
      end
      [ir,ic,v] = find(varargin{i});
      ic = ic+nc;
      iR = [iR;ir];
      iC = [iC;ic];
      V = [V;v];
      nc = nc+n(2);
    end
  end
  xOut = ssparse(iR,iC,V,nr,nc);