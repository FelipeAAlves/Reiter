%%% Description
%       Finds the invariant distribution
%
%%% INPUTS
%       (1) Pi : transpose of the Transitision matrix
function D = invdistr(Pi)

    % assert(all(abs(sum(Pi)-1)<1e-10));

    opts.disp=0;
    [eigenvector,e] = eigs(Pi,[],1,1+1e-10,opts);
    %[eigenvector,e] = eigs(Pi,[],1,'lm',opts);

    assert(abs(e-1)<1e-10);
    D = eigenvector/sum(eigenvector);

    %== check all elements are positive ==%
    assert(min(D)>-1e-12);
    
    %== Eigenvector METHOD ==%
    % Pi = full(Pi);
    % eta = min(Pi(Pi>0))/(2*size(Pi,1));
    % Pi = (Pi==0)*eta + Pi;
    % Pi = repmat(sum(Pi), size(Pi,1),1) .\ Pi;
    % assert( max(abs(sum(Pi)-1.0)) < 1e-10);
    %
    % [eigenvector2,e] = eigs(Pi,[],1,1+1e-10,opts);
    %
    % assert(abs(e-1)<1e-10);
    % D2 = eigenvector2/sum(eigenvector2);

    %== check all elements are positive ==%
    % assert(min(D2)>-1e-12);

    %== density ==%
    % D = max(D,0);
    % D = D/sum(D);
