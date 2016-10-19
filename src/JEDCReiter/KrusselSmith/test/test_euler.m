% set the Matlab path
% setpath;

% set the baseline parameters
% initModelParam;

function test_euler()
    global MP

    Kguess = 1.25 * MP.KRepSS;

    %== Recover prices ==%
    R      = 1 + netintr(Kguess,1);
    wage   = wagefunc(Kguess, 1);

    %== Compute Policy function for set of prices ==%
    [vSavingsPar, check] = broydn(@eulerres_stst, MP.SavingsParstart, [1e-11,1,1], R, wage);


    % falves notes
    % --------------
    %%% INPUT adjacob
    %       (1) func
    %       (2) X:          point at which to take derivatives
    %       (3) iDeriv :    index of which variables in X are to be differentiated
    %                       other parts of X are assumed constant.
    %       (4) nmax:       tells adjacob how many derivatives to take at a time.  I believe this has to do with
    %                       memory management.

    %== input for eulerres ==%
    X = [vSavingsPar; vSavingsPar; R ; R; wage; wage];

    iDeriv = 2 * length(vSavingsPar) + 1;
    nmax   = 1;

    J1 = adjacob( @equEuler, X, iDeriv,nmax );
end

function resid = equEuler(X)

    global MP

    neps = MP.neps;
    nSavingsPar = MP.nSavingsPar;

    %% Unpack X %%
    % *********************************************************************
    iSavingsPar = 1:neps*nSavingsPar;
    vSavingsPar     = X(iSavingsPar);
    
    iSavingsParNext = iSavingsPar + iSavingsPar(end);
    vSavingsParNext = X(iSavingsParNext);

    iPrices = iSavingsParNext(end)+1: length(X);
    Prices  = X(iPrices);

    R           = Prices(1);
    Rnext       = Prices(2);
    wage        = Prices(3);
    wagenext    = Prices(4);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %== evaluate eulerres ==%
    resid = eulerres(vSavingsPar, vSavingsParNext, R, Rnext, wage, wagenext);
end

%%% Description:
%       Evaluates the euler residuals at stst
function res = eulerres_stst(vSavingsPar, R, wage)

    res = eulerres(vSavingsPar, vSavingsPar, R, R, wage, wage);
end
