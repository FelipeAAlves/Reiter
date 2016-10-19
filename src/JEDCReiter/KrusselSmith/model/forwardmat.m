%%% Description
%       Function to compute transition matrix
%%% INPUT
%       (1)
%       (2) mSavingsPar :   parameters of the policy
%%% OBS: distribution is over beggining period capital and idio shock
function Pi = forwardmat(iFullMat, mSavingsPar)

    global MP;

    nd = MP.nHistogram;
    IR = []; IC = []; VV = [];

    for (ip=1:MP.neps)  % THIS period's state
        ir{ip} = [];
        ic{ip} = [];
        vv{ip} = [];

        %== Saving policy t ==%
        S = savingspline( mSavingsPar(:,ip) );

        %== Today assets ==%
        Athis = MP.AssetsGridFine;

        %== Interpolate saving ==%
        Aend = interp_savspline(S, Athis);

        %%%  NOTE:  lineartrans is an IMPORTANT fnc         %%%
        %%%      :  allocates to which state I am moving    %%%
        Ti = lineartrans(MP.AssetsGridFine, Aend);
        % Ti is a structure
        %   Ti.iFr:     indicates which position we come from
        %   Ti.iTo:     indicates to which position we go to
        %   Ti.Val:     weight of transition value

        ic{ip} = [ic{ip}; Ti.iFr];  % COLUMN indicates which position we come from
        ir{ip} = [ir{ip}; Ti.iTo];  % ROW indicates to which position we go
        vv{ip} = [vv{ip}; Ti.Val];


        if (iFullMat)
            offsi = (ip-1)*nd;
            for (jp=1:MP.neps)  % NEXT period's state

                pp = MP.transpp(ip,jp);
                if (pp>0)
                    offsj = (jp-1)*nd;
                    IC = [IC; offsi + ic{ip}];  % where come from! take off from offsi
                    IR = [IR; offsj + ir{ip}];  % where to go! take off from offsj
                    VV = [VV; pp*vv{ip}   ];  % intensity
                end
            end
        else
            Pij{ip} = sparse(ir{ip},ic{ip},vv{ip},nd,nd);
        end
    end

    if (iFullMat)
        nn = nd*MP.neps;
        Pi = sparse(IR,IC,VV,nn,nn);
    else

        Pi = op_concat( op_kron(MP.transpp',speye(nd)), blockdiag(Pij) );
        %%%  WARN:  the ordering here matters for sequence of taking saving decision     %%%
        %%%         here we have a different order from what is in Reiter's CODE         %%%

    end
end

%
%%% Description:
%       Construct the transition probability matrix
%
%%% INPUTS:
%   kgrid   : grid for state variable
%   k0      : policy evaluated at the grid
function S = lineartrans(kgrid,k0)

    n = length(kgrid);
    k0 = max(min(k0,kgrid(n)),kgrid(1));

    iPos = lookup(kgrid,k0,0);
    iPos = reshape(iPos,1,length(k0));  % make sure it is column vector
    iPos = min(iPos,n-1);

    pHigh = (k0-kgrid(iPos))./(kgrid(iPos+1)-kgrid(iPos));
    pHigh = reshape(pHigh,1,length(k0));  % make sure it is column vector

    S.iFr = reshape(repmat(1:n,2,1),2*n,1);
    S.iTo = reshape([iPos;iPos+1],2*n,1);
    S.Val = reshape([1-pHigh;pHigh],2*n,1);
end
