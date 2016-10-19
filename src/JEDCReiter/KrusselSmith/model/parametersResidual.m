function [value, vDerivs] = parametersResidual(vParameters,aGridMoments)
%%% Description
%       Computes the objective function and Jacobian for computing the parameters of the distribution, given
%       moments to match
%
%%% INPUT
%       (1) vParameters:  candidate parameters of distribution (to be optimized over)
%       (2) aGridMoments: (nAssetsQuad x nMoments) grid of centralized moments that enter on pdf
%
%%% OUTPUT
%       (1) value: value of objective function
%       (2) vDerivs: Jacobian
%

global MP

% Value of function
value = MP.QuadWeights' * exp( aGridMoments * vParameters );

% Derivatives of function
if nargout > 1
    % OLD expression
    % vDerivs = sum(repmat(MP.QuadWeights,[1 MP.nMoments]) .* ...
    % aGridMoments .* ...
    % repmat(exp( aGridMoments * vParameters ),...
    %     [1 MP.nMoments]),1)';

    vDerivs = zeros(MP.nMoments,1);
    for iMoment = 1 : MP.nMoments
        vDerivs(iMoment) = MP.QuadWeights' * ( exp( aGridMoments * vParameters ) .* aGridMoments(:,iMoment) );
    end
end
