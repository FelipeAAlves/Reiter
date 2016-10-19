function caliinc(neps)
    global MP;

    MP.neps = neps;


    if (MP.neps==1)
        MP.ypp = 1;
        MP.transpp = 1;

        MP.ydist = 1
    elseif (MP.neps==2)
        %MP.ypp = [1 1];
        MP.ypp = [0.5 1.5];
        MP.transpp = [0.9 0.1; 0.5 0.5];

        [v,L] = eig(MP.transpp );

        MP.ydist = v(:,MP.neps) / sum(v(:,MP.neps));
    elseif (MP.neps==3)
        %MP.ypp = [1 1];
        MP.ypp = [0.2 1 5];
        MP.transpp = [0.97 0.03 0; 0.03 0.94 0.03; 0 0.03 0.97];
        [v,L] = eig(MP.transpp );

        MP.ydist = v(:,MP.neps) / sum(v(:,MP.neps));

    elseif (MP.neps>10)
        n = MP.neps;
        assert(mod(n,2)==1);
        nmed = (n+1)/2;  %index of median



        MP.ypp = exp(linspace(-1.2,1.2,n));
        % make biggest income group very big:
        % MP.ypp(n) = MP.ypp(n)*(MP.ypp(2)/MP.ypp(1))^(fac-1);

        pSwitch = sqrt(0.015)/diff(log(MP.ypp(1:2)))/MP.freq;
        pSwitch2 = pSwitch/2;
        x = [pSwitch;pSwitch2;1.5];
        [x,iflag] = broydn(@cali_z,x,[1e-7,0,1],nmed);
        [resid,Pi] = cali_z(x,nmed);

        % make it global:
        MP.transpp = Pi;
        MP.distpp = invdistr(Pi')';
    else
        error('wrong neps');
    end




function [resid,Pi] = cali_z(x,nmed)
    global MP;

    pSwitch = x(1);
    pSwitch2 = x(2);
    expo = x(3);
    if (expo<1)
    resid = 1e100;
    return;
    end


    % transition probability of survivors
    PiHH = make_pi(pSwitch,pSwitch2,MP.neps);
    Pi4 = PiHH^MP.freq;
    pdeath = 0.02/MP.freq;
    PiDeath = eye(MP.neps)*(1-pdeath);
    PiDeath(:,nmed) = PiDeath(:,nmed) + pdeath;
    Pi = PiDeath*PiHH;
    MP.distpp = invdistr(Pi')';

    frac1 = 0.1476;
    MP.ypp = exp(linspace(-1.2,1.2,MP.neps));
    nbig = ceil(4*MP.neps/5);  %index of median
    MP.ypp(nbig:end) = MP.ypp(nbig:end).^expo;
    ytop = 100*frac1/MP.xY(end);
    ytarget = (1-frac1-MP.distpp(end)*MP.wY(1)*ytop*MP.xY(1)); %/(0.99-MP.distpp(end)*MP.wY(1));
    yactual = dot(MP.distpp(1:end-1),MP.ypp(1:end-1));
    MP.ypp =  MP.ypp * ytarget / yactual;
    MP.ypp(end) = ytop;

    y = log(MP.ypp(:));
    pi = Pi4(nmed,:);
    Ey = pi*y;
    Vary = pi*(y-Ey).^2;
    Stdy = sqrt(Vary);

    Ymat = MP.ypp(:)*MP.xY(:)';
    Pmat = MP.distpp(:)*MP.wY(:)';
    [Yqu,Yshares] = quints(Ymat,Pmat);


    resid = [Vary - 0.015;
    MP.distpp(end)*MP.wY(end)-0.01;
    Yshares(5) - 0.6139];
    %sum(Yshares(7:8)) - (0.1637+0.1476)];
    % Yshares(8) - 0.1476];


function Pi = make_pi(pSwitch,pSwitch2,n)
Pi = eye(n);
for i=1:n
if(i==1)
Pi(i,i+1) = Pi(i,i+1) + pSwitch;
Pi(i,i) = Pi(i,i) - pSwitch;
elseif(i==n)
Pi(i,i-1) = Pi(i,i-1) + pSwitch;
Pi(i,i) = Pi(i,i) - pSwitch;
elseif(i==n-1)
Pi(i,i+1) = Pi(i,i+1) + pSwitch2;
Pi(i,i-1) = Pi(i,i-1) + pSwitch;
Pi(i,i) = Pi(i,i) - pSwitch - pSwitch2;
else
Pi(i,i+1) = Pi(i,i+1) + pSwitch;
Pi(i,i-1) = Pi(i,i-1) + pSwitch;
Pi(i,i) = Pi(i,i) - 2*pSwitch;
end
end


function [resid,w] =  inciid_z(x)
global MP;
if(x(1)>=x(2))
resid = 1e100;
end
if(MP.neps==1)
p = [0.5 0.5];
else
% p = [0.06 0.94];
p = [0.25 0.75];
p = [0.5 0.5];
end
logx = log(x);
Elogx = p*logx;
sigY = sqrt(0.061/MP.freq);
resid = [p*x - 1;0];
if(MP.neps==1)
resid(2) = p*(logx-Elogx).^2 - sigY;
else
resid(2) = x(1) - 0.2;
end
w = p';
