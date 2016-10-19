iEBW = Env.iEquBWS;
iEFW = Env.iEquFW;
iVBW = Env.iVarBWS;
iVFW = Env.iVarDec;

Ss0 = Env.jac0(iEBW,iVBW);
Sd0 = Env.jac0(iEBW,iVFW);
Es0 = Env.jac0(iEFW,iVBW);
Ed0 = Env.jac0(iEFW,iVFW);
Ss1 = Env.jac1(iEBW,iVBW);
Sd1 = Env.jac1(iEBW,iVFW);
Es1 = Env.jac1(iEFW,iVBW);
Ed1 = Env.jac1(iEFW,iVFW);
Se0 = Env.jace(iEBW,:);
Ee0 = Env.jace(iEFW,:);

T = -Ss0\Ss1;
D = -Ss0\Sd0;
Se= -Ss0\Se0;

