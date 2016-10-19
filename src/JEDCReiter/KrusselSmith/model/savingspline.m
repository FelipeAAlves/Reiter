%%% Description 
%
%
%%% INPUT:
%   Par:        finite parameterization of the savings function
%
%%% OUTPUT:
%   x:          knot points of spline
%   y:          saving at knot points
%   Slope:      derivative of saving at knot points
%
%   The output define a spline approximation for the savings function
function S = savingspline(Par)

    global MP;

    % if (nargin==2)
    %     Par = [X(1);Par(2:end)];
    %     knotXi = X(2:end)-X(1);
    % else
    %     knotXi = MP.knotXi;
    % end

    knotXi = MP.knotXi;

    S.xcrit = Par(1);

    S.x = [ S.xcrit ; S.xcrit + knotXi; S.xcrit+1e8 ];  %add end point 1e8/
    S.y = [ 0; Par(2:end); 0 ];

    n = length(S.y);
    S.y(n) = S.y(n-1) + (S.x(n)-S.x(n-1))*(S.y(n-1)-S.y(n-2))/(S.x(n-1)-S.x(n-2));

    % S.iSlope = diff(S.x) ./ diff(S.y);
    S.Slope = diff(S.y) ./ diff(S.x);
