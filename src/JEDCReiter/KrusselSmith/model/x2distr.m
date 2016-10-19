
%%% Description
%       Converts x() to a vHistogram
%
%%% INPUTS
%       (1) xdist
%       (2) Env OR SS variable
function D =  x2distr(vHistogramX,vHistogramSS)

    %== Initilize D as deriv1 ==%
    D = initsize(vHistogramSS, vHistogramX);
    D(:) = vHistogramSS;                %stst Histogram Vals

    vHistogramDev = [ - sum(vHistogramX); vHistogramX];

    D = D + vHistogramDev;

    %== Confirm that the Distribution integrates to 1 ==%
    assert(abs(sum(D) - 1)<1e-10,'assert failed in x2distr: distribution does not integrate to one.');
