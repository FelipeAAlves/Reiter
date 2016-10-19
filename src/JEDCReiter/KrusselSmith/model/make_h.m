function H = make_h(ii,n)
    global MP;
    if (nargin<2)
        n = MP.nHistogram-1;
    end
    H = zeros(n,length(ii));
    for j=ii
        h = expect_k([],j);
        h = reshape(h,MP.nHistogram,MP.neps);
        m = length(h)-1;
        h2 = h(2:end,:) - repmat(h(1,:),m,1);
        H(:, j) = h2(:);
  end
%%%  REVIEW:  check ndstst     %%%
