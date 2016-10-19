% Matlab file for_ paper
% "Approximate and Almost-Exact Aggregation in Dynamic Stochastic Heterogeneous-Agent Models"
% Michael Reiter, IHS, December 2009
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
% .....................................................................................

% Matlab file to implement class "deriv1",
% forward mode of automatic differentation, first derivatives, possibly ssparse

%%% Usage:
%   1) x = deriv1(val): val is column vector
%          x is vector of independent variables, Jacobian is the identity
%   2) x = deriv1(val,iIndep): val is column vector, iIndep index of variables
%          x is vector where iIndep are independent variables, rest are constants
%   3) x = deriv1(val,0,np): val any matrix, np positive integer
%          x has zero Jacobian, assuming np independent variables
%   4) x = deriv1(val,[],J): val any matrix, J the Jacobian
%          x is matrix of values with Jacobian J
classdef deriv1
    properties
        v
        d
    end
    methods
        function obj= deriv1(val, iIndep, deriv)

            if isa(val,'deriv') || isa(val,'struct')
                obj.v = val.v;
                obj.d = val.d;
                return;
            end

            nSparse = 1;
            if (nargin==0 | nargin>3)
                error('usage: deriv1(Value,Derivative)');
            end
            obj.v=val;

            if nargin==3
                if (isscalar(iIndep) && iIndep==0)
                    np = deriv;
                    if (np>=nSparse)
                        obj.d = ssparse(numel(val),np);
                    else
                        obj.d = zeros(numel(val),np);
                    end
                else
                    obj.d=ssparse(deriv);
                end
            % only one argument
            else
                m=size(val);
                if (length(m)~=2 | m(2)~=1)
                    error('single argument to deriv1 must be column vector (independent variables)');
                end
                if (nargin==1)
                    iIndep = 1:m(1);
                end

                nIndep = length(iIndep);
                if (nIndep>=nSparse | issparse(val))
                    obj.d=ssparse(iIndep,1:nIndep,ones(nIndep,1),m(1),nIndep);
                else
                    obj.d=zeros(m(1),nIndep);
                    obj.d(iIndep,:)=eye(nIndep);
                end
            end
        end %deriv1
% other methods
%         xOut = abs(xIn)
%         xOut = acos(xIn)
%         xOut = acosh(xIn)
%         b = all(m)
%         xOut = and(x1,x2)
%         b = any(m)
%         xOut = asinh(xIn)
%         xOut = atan(xIn)
%         xOut = atanh(xIn)
%         xOut = cos(xIn)
%         xOut = cosh(xIn)
%         xOut = transpose(xIn)
%         d = diff(m)
%         display(s)
%         dissparse(x)
%         e = end(obj,k,n)
%         xOut=eq(x1,x2)
%         xOut=exp(xIn)
    end
end
