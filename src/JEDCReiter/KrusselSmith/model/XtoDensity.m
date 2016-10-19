
%%% Description
%       Converts X to mDensityCoeff
%
%%% INPUTS
%       (1) Xdens
%       (2) Env variable
function [mDensity, mDensityCoeff] =  XtoDensity(mMoments, XDensityCoeff)

    global MP;

    if size(mMoments,2)==1
       mMoments = reshape(mMoments, MP.nMoments, MP.neps);
    end
    %== Initilize aGridMoments,mDensityCoeff,mDensity as deriv object ==%
    aGridMoments = initsize(zeros( MP.nAssetsQuad, MP.nMoments), mMoments);

    mDensityCoeff  = initsize(zeros( MP.nDensityCoeff, MP.neps ), mMoments, XDensityCoeff);
    mDensity       = initsize(zeros( MP.nAssetsQuad, MP.neps ), mMoments, XDensityCoeff);

    %== Allocate XDensityCoeff ==%
    mDensityCoeff(2:end,:) = reshape(XDensityCoeff, MP.nDensityCoeff-1, MP.neps);

    for ieps = 1 : MP.neps

        aGridMoments(:,1)   = MP.AssetsGridQuad - mMoments(1,ieps);

        % Higher order moments (centered)
        for iMoment = 2 : MP.nMoments

            % aMoments on the QuadNodes
            aGridMoments(:,iMoment) = ...
            (MP.AssetsGridQuad - mMoments(1,ieps) ) .^ iMoment - mMoments(iMoment,ieps);

        end

        %== Recover g₀ ==%
        normalization = parametersResidual(mDensityCoeff(2:end,ieps), aGridMoments);

        %== density Parameters ==%
        mDensityCoeff(1,ieps) = 1/normalization;

        %== Density at Nodes ==%
        mDensity(:,ieps)    = 1/normalization * exp( aGridMoments * mDensityCoeff(2:end,ieps) );
    end

end
