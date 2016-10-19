
%%% Description:
%       Computes the jacobian through automatic difference.
%
%%% INPUT
%       () func  :      function to be evaluated
%       () X     :      argument of fnc at which derivative is takes
%       () iDeriv:      index of which variables in X are to be diff
%       () nmax  :      maximum number of derivatives to take at EACH ITER

function J = adjacob(func, X, iDeriv, nmax, varargin)

    nX = length(X);
    if (isempty(iDeriv))
        iDeriv = 1:nX;
    end;


    %== Set the derivatives ==%
    nder = length(iDeriv);              % # derivatives
    Ind = sparse(nX,nder);              % Ind of derivatives to take, Nx variables and nder derivatives
    Ind(iDeriv,:) = speye(nder);
    ndone = 0;

    %== Preallocate J ==%
    nEqu = size( func(X, varargin{:}), 1);
    J = zeros(nEqu, nder);

    % fprintf('adjacob: %d to do\nder',nder);
    while (ndone<nder)
        if (isscalar(nmax))
            ii = ndone+1:min(nder,ndone+nmax);
        else
            %== decide how many derivatives to take ==%
            nmax0 = 1000;
            ii0 = ndone + 1:min(nder,ndone+nmax0);
            nCumul = cumsum(nmax(ii0));
            ii = ndone + (1:find(nCumul<=nmax0, 1,'last')); %last entry with 1
            assert(~isempty(ii));
        end

        %== INITIALIZE an deriv1 instance with X values ==%
        xx = deriv1(X,[],Ind(:,ii));
        cput = cputime;

        %== CALL fnc on xx ==%
        ff = func(xx,varargin{:});

        if (isa(ff,'deriv1'))
            %== Get Jacobian ==%
            j = getjac(ff);
        else
            j = zeros(size(ff,1),length(ii));
        end
        ndone = max(ii);
        J(:,ii) = j;
        % fprintf('%d done in %0.2f sec\nder',ndone,cputime-cput);
    end
