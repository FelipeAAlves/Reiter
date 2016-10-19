

%%% Description:
%       Computes the residual of market-clearing condition, using histogram
%       approximation of distribution. The law of motion for histogram is done as in Young (2010);
%
%%% INPUT:
%       (1) K:                      candidate aggregate capital stock
%       (2) (opt) mSavingsPar:      take in policy coefficients, o/w create
%
%%% OUTPUT
%       (1) residual:               residual of market clearing condition
%       (2) (opt) SavingsPar:       solution SavingsParameters for given prices
%       (3) (opt) mHistogramOpt:    mHistogram - histogram for distribution

function [resid, mSavingsParOpt, mHistogramOpt, SSOpt] = ststHistogramResid(K, mSavingsParOpt)

    global MP;

    fprintf(1,' OUTER loop (Prices): \n');
    fprintf(1,'  Kguess = %.4f \n', K);

    %%%  NOTE:  given a GUESS for Kequ, recover prices           %%%
    %%%         and solve for consumption pol by collocation     %%%

    %== Recover prices ==%
    R      = 1 + netintr(K,1);
    wage   = wagefunc(K, 1);

    %=======================================================%
    %%% Compute Policy function for set of prices         %%%
    %=======================================================%
    if nargin == 1

        disp(['   INNER LOOP (Household Policy)'])
        %%%  NOTE:  line 45 determines if Parstart changes or not     %%%
        [vSavingsPar, check] = broydn(@eulerres_stst,MP.SavingsParstart,[1e-11,1,1], R, wage);
        if ( check~=0 )
            warning('broydn not converged');
        end

        mSavingsPar = reshape(vSavingsPar, MP.nSavingsPar, MP.neps);
        MP.SavingsParstart = vSavingsPar;
    else
        mSavingsPar = mSavingsParOpt;
    end

    %=================================================================%
    %%% Compute Stationary Distribution from decision rules         %%%
    %=================================================================%

    %== Compute Transition matrix (iFullMat) ==%
    % Pi1  = forwardmat(1,SavingsPar);
    Pi = sparse( forwardmat(0, mSavingsPar) );

    %== Invariant distribution ==%
    vHistogram = invdistr(Pi);

    %== CHECK ==%
    err_invDist = max( abs(Pi*vHistogram - vHistogram) );
    fprintf('    Invariant distribution Converged to %.4e\n', err_invDist)
    assert(err_invDist<1e-8, 'Invariante distribution is incorrect')
    clear Pi;

    %== Plot Distribution ==%
    %mHistogram = reshape(vHistogram,MP.nHistogram,MP.neps);
    %figure
    %plot(MP.AssetsGridFine, mHistogram(:,1), 'linewidth', 1.5)
    %figure
    %plot(MP.AssetsGridFine, mHistogram(:,2), 'linewidth', 1.5)

    Ksupply = expect_k(vHistogram);

    resid = Ksupply - K;

    fprintf('\n')
    fprintf(1,'  Kguess  = %0.5f;  \n', K);
    fprintf(1,'  Ksupply = %0.5f;  \n', Ksupply);
    fprintf(1,'  resid = %0.5f;  \n', resid);
    disp(['.......................................................................'])


    if (nargout>1)

        %== Set Opt outputs ==%
        mSavingsParOpt = mSavingsPar;
        mHistogramOpt  = reshape( vHistogram, MP.nHistogram, MP.neps);

        if (nargout>3)

            %== states ==%
            SSOpt.vHistogram   = vHistogram;
            SSOpt.K            = Ksupply;

            %== policy ==%
            SSOpt.mSavingsPar = mSavingsPar;

            %== Prices ==%
            SSOpt.R     = R;
            SSOpt.wage  = wage;
        end
    end
end


%%% Description:
%       Evaluates the euler residuals at stst
function res = eulerres_stst(vSavingsPar, R, wage)

    res = eulerres(vSavingsPar, vSavingsPar, R, R, wage, wage);
end
