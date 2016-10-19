function [ser, shocks] = SimulateSystem(A, B, T, Sigma, exerc, seed)

    if nargin == 4
        exerc = 'simulation';
    end
    
    if nargin == 6
        if verLessThan('matlab', '7.12')
            warning('calcvar:oldmatlab','Setting random seeds in this way is not supported in Matlab releases before R2011a.  Seed will not be set.');
        else
            rng('default')
            rng(seed)
        end

    end

    nvar   = size(B,1);
    nshock = size(B,2);

    switch exerc
        case 'simulation'
            
            % init series             
            ser    = zeros(nvar,T);
            
            % init shocks            
            shocks = randn(T,nshock)'  ;  %the RNG populates the matrix column by column so forming the matrix this way and then transposing it means that we get the same shocks even if we leave off some columns.
            
            state  = zeros(nvar,1);
            for t = 1:T
                state    = A * state + B * Sigma * shocks(:,t);
                ser(:,t) = state;
            end

        case 'irf'
            
            % init series 
            ser    = cell(nshock);
            shocks = cell(nshock);
            
            % iterate on the shock to do IRF
            for ishock =1:nshock
                
                % init series                 
                ser{ishock}    = zeros(nvar,T);
                
                % init shocks                 
                shocks{ishock} = zeros(nshock,T);
                shocks{ishock}(ishock,1) = 1;
                
                state  = zeros(nvar,1);
                for t = 1:T
                    state    = A * state + B * Sigma *shocks{ishock}(:,t);
                    ser{ishock}(:,t) = state;
                end                
                
            end
    end
    
    %== transpose result ==%
    % ser = ser';
