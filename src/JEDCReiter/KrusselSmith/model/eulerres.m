%
%%% Description
%       Computes the euler residual at collocation points
%
%%% INPUT:
%    Par:               parameter vector for_ current savings function_
%    Parnext:           parameter vector for_ next periods' savings function_
%    R, wage, Transf    prices
%    athis:             grid at which_ to compute Euler residuals; default: knot points of spline
%
function resid = eulerres(Par, Parnext, R, Rnext, wage, wagenext, athis)

    global MP;

    neps = MP.neps;
    nSavingsPar = MP.nSavingsPar;

    %== Reshape ==%
    Par     = reshape(Par    , nSavingsPar,neps);
    Parnext = reshape(Parnext, nSavingsPar,neps);

    %== INITIALIZE resid - allows it to receive an deriv object ==%
    resid = initsize( zeros(nSavingsPar,neps), Par, Parnext, R, Rnext, wage, wagenext);

    %== loop over idio shock ==%
    for ieps=1:neps

        %== Saving policy t ==%
        S = savingspline( Par(:,ieps) );

        if (nargin>6) % athis given
            sthis = interp_savspline(S,athis);
        else
            athis = S.x(1:end-1);  % last point in S is for interpolation
            sthis = S.y(1:end-1);
        end

        %== GET consumption ==%
        C = getConsumption(athis, sthis, R, wage, ieps);

        if ( any(C<0) )  %signal inadmissible value to routine 'broydn.m';
            error('inadmiss')
            %         ieps
            %         I = find(cthis<0)
            %         nthis(I)
            %         xthis(I)
            %         sthis(I)
            resid = 1e100;
            return;
        end

        %== Next period assets ==%
        anext = sthis;
        MUexp = 0;

        %== COMPUTE Exp over idio shock ==%
        for jeps=1:neps
            prob_i2jeps = MP.transpp(ieps,jeps);
            if (prob_i2jeps>0)

                %== Saving policy t+1 ==%
                Sn = savingspline(Parnext(:,jeps));

                %== Saving at values ==%
                snext = interp_savspline(Sn, anext);

                %== Consumption next period ==%
                Cnext = getConsumption(anext, snext, Rnext, wagenext, jeps);

                MUnext = margutil(Cnext);
                MUexp = MUexp + prob_i2jeps * (MUnext);
            end
        end

        %%  NOTE: Euler residual: expressed in relative     %%
        %%        consumption units                         %%
        %%        NOTICE that R is incorporated into MUexp  %%
        resid(:,ieps) = 1 - inv_margutil(MP.beta * Rnext *MUexp)./C;
    end

    %== vectorize resid ==%
    resid = resid(:);
end
