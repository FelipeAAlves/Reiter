
%%% Description
%       Computes the objects needed for the linearized system

function Env = dojacob(funcname, X, Env, njac, varargin)

    if (nargin<4)
        njac = 1000;
    end

    % check adjacob fnc to understand what is going on!!
    Env.jac0 = sparse( adjacob( funcname, X, Env.iVarX.x   , njac, varargin{:}) );
    Env.jac1 = sparse( adjacob( funcname, X, Env.iVarX.xlag, njac, varargin{:}) );
    Env.jace = sparse( adjacob( funcname, X, Env.iVarX.eps , njac, varargin{:}) );
    Env.jact = sparse( adjacob( funcname, X, Env.iVarX.eta , njac, varargin{:}) );

%  assert(nnz(Env.jac1(Env.iEquBWS,Env.iVarDec))==0);
    if (isfield(Env,'iVarJump'))
        assert(nnz(Env.jac1(:,Env.iVarJump))==0);
    end
