function ser = ir_sims(G1, Impact, T, H, scale_shock)

    if (nargin<5 | isempty(scale_shock))
        scale_shock = 0.01;
    end


    nz    = size(Impact,2);
    Sigma = eye(nz);    %irrelevant

    %== loop over shocks ==%
    for i = 1 : nz
        shocks = zeros(nz,T+1);

        if length(scale_shock) == 1
            shocks(i,1) = scale_shock;
        else
            shocks(i,1) = scale_shock(i,i);
        end

        x = simul_sims(G1, Impact, shocks, Sigma, H);

        ser{i} = x(2:end,:);
    end
