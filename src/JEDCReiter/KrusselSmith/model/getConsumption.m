
%%% Description:
%       Computes consumption from the budget constraint
% 
function ct = getConsumption(a, savings, R, wage, incomeIndex)

    global MP;
    if ~( size(a) == size(savings) )
        error('Size of assets and Savings must be the same');
    end

    unempInd = (incomeIndex==1);
    ct = R * a + wage * ( MP.mu * unempInd  + (1 - MP.tau) * (1 - unempInd) ) - savings;
end
