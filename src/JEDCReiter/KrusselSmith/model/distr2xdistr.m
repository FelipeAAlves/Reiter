%%% Description
%       Converts the vHistogram to vHistogramX
%
%%% INPUT
%       (1) vHistogram  : vHistogram
%       (2) vHistogramSS: 
function vHistogramX =  distr2xdistr(vHistogram, vHistogramSS)

    %== Compute the deviation from Ss ==%
    vHistogramDev = vHistogram - vHistogramSS;

    vHistogramX = vHistogramDev(2:end);
