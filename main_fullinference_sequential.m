close all; clear all;
% Linh Huynh, May 15, 2024
% Sequential Inference

%% Part 1: Mono-Cultured Data
%% Simulated Data 
% Parameters
b0S     = 1.1/120;
%b0S = 0.5;
rS      = 1/120;
KS      = 1e3;
gammaS  = 0.7;
d0S     = b0S-rS;
alphaS  = 0;
sigmaS  = 0.5;

b0R     = 1.1/120;
rR      = 1/120;
KR      = 1e3;
gammaR = 0.5;
d0R     = b0R-rR;
alphaR  = 0;
sigmaR  = 0.5;

% Functions
bfunc_S = @(S,R) b0S - gammaS*(rS/KS).*S - alphaS*sigmaS*(rS/KS).*R;
dfunc_S = @(S,R) d0S + (1-gammaS)*(rS/KS).*S + alphaS*(1-sigmaS)*(rS/KS).*R;
bfunc_R = @(S,R) b0R - gammaR*(rR/KR).*R - alphaR*sigmaR*(rR/KR).*S;
dfunc_R = @(S,R) d0R + (1-gammaR)*(rR/KR).*R + alphaR*(1-sigmaR)*(rR/KR).*S;

% Simulation
nrun = 50;
tmax = 3*1e3;
dt   = 1/30;
S0 = 1;
R0 = 1;
[Smat,Rmat,t,dt] = datasimulation_Langevin(bfunc_S,dfunc_S,bfunc_R,dfunc_R,nrun,tmax,dt,S0,R0);

% plot
for traj = 1:nrun
    figure(1)
    plot(t,Smat(traj,:))
    hold on
    plot(t,Rmat(traj,:))
    xlabel('Time')
    ylabel('Cell Number')
    hold on
end

%% Inference of birth and death rates from mono-cultured data
Sbinsize = 10;
Rbinsize = 10;
%[binleftindexSRvec] = separatebirthdeathrates_2D(Smat,Rmat,dt,Sbinsize,Rbinsize);
[bSS_estimated,dSS_estimated,bRR_estimated,dRR_estimated,binedgesS,binedgesR,IDS,IDR,dSvarvec,dRvarvec,dSmoment3cen,dRmoment3cen] = separatebirthdeathrates_2D(Smat,Rmat,dt,Sbinsize,Rbinsize);

S_binleft = (binedgesS(IDS))';
R_binleft = (binedgesR(IDR))';

bSS_true = bfunc_S(S_binleft,R_binleft);
bSS_true = bSS_true.*S_binleft;
error_bSS = norm(bSS_true-bSS_estimated)./norm(bSS_true);

figure(2)
plot3(S_binleft,R_binleft,bSS_estimated,'r.')
hold on
plot3(S_binleft,R_binleft,bSS_true,'b.')
xlabel('S-Population Size')
ylabel('R-Population Size')
zlabel('Rate')
legend('inferred','true')

figure(6)
plot(S_binleft,bSS_estimated,'r.')
hold on
plot(S_binleft,bSS_true,'b.')
legend('inferred','true')


%% dSS
dSS_true = dfunc_S(S_binleft,R_binleft);
dSS_true = dSS_true.*S_binleft;
error_dSS = norm(dSS_true-dSS_estimated)./norm(dSS_true);

figure(3)
plot3(S_binleft,R_binleft,dSS_estimated,'r.')
hold on
plot3(S_binleft,R_binleft,dSS_true,'b.')
xlabel('S-Population Size')
ylabel('R-Population Size')
zlabel('Rate')
legend('inferred','true')

figure(7)
plot(S_binleft,dSS_estimated,'r.')
hold on
plot(S_binleft,dSS_true,'b.')
legend('inferred','true')


%% bRR
bRR_true = bfunc_R(S_binleft,R_binleft);
bRR_true = bRR_true.*R_binleft;
error_bRR = norm(bRR_true-bRR_estimated)./norm(bRR_true);

figure(4)
plot3(S_binleft,R_binleft,bRR_estimated,'r.')
hold on
plot3(S_binleft,R_binleft,bRR_true,'b.')
xlabel('S-Population Size')
ylabel('R-Population Size')
zlabel('Rate')
legend('inferred','true')

figure(8)
plot(R_binleft,bRR_estimated,'r.')
hold on
plot(R_binleft,bRR_true,'b.')
legend('inferred','true')


%% dRR
dRR_true = dfunc_R(S_binleft,R_binleft);
dRR_true = dRR_true.*R_binleft;
error_dRR = norm(dRR_true-dRR_estimated)./norm(dRR_true);

figure(5)
plot3(S_binleft,R_binleft,dRR_estimated,'r.')
hold on
plot3(S_binleft,R_binleft,dRR_true,'b.')
xlabel('S-Population Size')
ylabel('R-Population Size')
zlabel('Rate')
legend('inferred','true')

figure(9)
plot(R_binleft,dRR_estimated,'r.')
hold on
plot(R_binleft,dRR_true,'b.')
legend('inferred','true')


%% Inference of r and K from monoculture data (INCORPORATE HEYRIM'S CODE HERE)

 rS_estimated = rS;
 KS_estimated = KS;
 rR_estimated = rR;
 KR_estimated = KR;

%% Inference density dependence parameter and intrinsic rates (rough estimation)
%b0SS_estimated = b0S.*S_binleft;
%b0S_estimated  = b0SS_estimated./S_binleft;
%d0SS_estimated = d0S.*S_binleft;
%d0S_estimated  = d0SS_estimated./S_binleft;
%gammaS_estimated = (KS_estimated).*(b0SS_estimated - bSS_estimated)./(rS_estimated.*S_binleft.*S_binleft);

%b0RR_estimated = b0R.*R_binleft;
%b0R_estimated  = b0RR_estimated./R_binleft;
%d0RR_estimated = d0R.*R_binleft;
%d0R_estimated  = b0RR_estimated./R_binleft;
%gammaR_estimated = (KR_estimated).*(b0RR_estimated - bRR_estimated)./(rR_estimated.*R_binleft.*R_binleft);


%% Inference of intra-species density dependence parameter gammaS and intrinsic rates b0S
fun = @(unknownparavecS) dot(momenteqS(unknownparavecS,rS_estimated,KS_estimated,S_binleft,dSvarvec,dSmoment3cen,dt),momenteqS(unknownparavecS,rS_estimated,KS_estimated,S_binleft,dSvarvec,dSmoment3cen,dt));
b0S_guess = 0.1;
gammaS_guess = 0.1;

x0  = [b0S_guess,gammaS_guess];
A   = [];
b   = [];
Aeq = [];
beq = [];
lb  = [0,0];
ub  = [10,1];
sol_estimatedb0gammaS = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
b0S_estimated = sol_estimatedb0gammaS(1);
gammaS_estimated = sol_estimatedb0gammaS(2);



%% Inference of intra-species density dependence parameter gammaR and intrinsic rates b0R
funR = @(unknownparavecR) dot(momenteqR(unknownparavecR,rR_estimated,KR_estimated,R_binleft,dRvarvec,dRmoment3cen,dt),momenteqR(unknownparavecR,rR_estimated,KR_estimated,R_binleft,dRvarvec,dRmoment3cen,dt));
b0R_guess = 0.1;
gammaR_guess = 0.1;

x0  = [b0R_guess,gammaR_guess];
A   = [];
b   = [];
Aeq = [];
beq = [];
lb  = [0,0];
ub  = [10,1];
sol_estimatedb0gammaR = fmincon(funR,x0,A,b,Aeq,beq,lb,ub);
b0R_estimated = sol_estimatedb0gammaR(1);
gammaR_estimated = sol_estimatedb0gammaR(2);


%% Part 2: Co-cultured Data
%% Simulated Data 
% Parameters
b0S     = 1.1/120;
rS      = 1/120;
KS      = 1e3;
gammaS  = 0.7;
d0S     = b0S-rS;
alphaS  = 0.5;
sigmaS  = 0.5;

b0R     = 1.1/120;
rR      = 1/120;
KR      = 1e3;
gammaR  = 0.5;
d0R     = b0R-rR;
alphaR  = 0.5;
sigmaR  = 0.5;

% Functions
bfunc_S = @(S,R) b0S - gammaS*(rS/KS).*S - alphaS*sigmaS*(rS/KS).*R;
dfunc_S = @(S,R) d0S + (1-gammaS)*(rS/KS).*S + alphaS*(1-sigmaS)*(rS/KS).*R;
bfunc_R = @(S,R) b0R - gammaR*(rR/KR).*R - alphaR*sigmaR*(rR/KR).*S;
dfunc_R = @(S,R) d0R + (1-gammaR)*(rR/KR).*R + alphaR*(1-sigmaR)*(rR/KR).*S;

% Simulation
nrun = 50;
tmax = 3*1e3;
dt   = 1/30;
S0 = 1;
R0 = 1;
[Smat,Rmat,t,dt] = datasimulation_Langevin(bfunc_S,dfunc_S,bfunc_R,dfunc_R,nrun,tmax,dt,S0,R0);

% plot
for traj = 1:nrun
    figure(10)
    plot(t,Smat(traj,:))
    hold on
    plot(t,Rmat(traj,:))
    xlabel('Time')
    ylabel('Cell Number')
    hold on
end

%% Inference of birth and death rates from co-cultured data
Sbinsize = 10;
Rbinsize = 10;
%[binleftindexSRvec] = separatebirthdeathrates_2D(Smat,Rmat,dt,Sbinsize,Rbinsize);
[bSS_estimated,dSS_estimated,bRR_estimated,dRR_estimated,binedgesS,binedgesR,IDS,IDR,dSvarvec,dRvarvec,dSmoment3cen,dRmoment3cen] = separatebirthdeathrates_2D(Smat,Rmat,dt,Sbinsize,Rbinsize);

S_binleft = (binedgesS(IDS))';
R_binleft = (binedgesR(IDR))';

bSS_true = bfunc_S(S_binleft,R_binleft);
bSS_true = bSS_true.*S_binleft;
error_bSS = norm(bSS_true-bSS_estimated)./norm(bSS_true);

figure(11)
plot3(S_binleft,R_binleft,bSS_estimated,'r.')
hold on
plot3(S_binleft,R_binleft,bSS_true,'b.')
xlabel('S-Population Size')
ylabel('R-Population Size')
zlabel('Rate')
legend('inferred','true')

figure(12)
plot(S_binleft,bSS_estimated,'r.')
hold on
plot(S_binleft,bSS_true,'b.')
legend('inferred','true')


%% dSS
dSS_true = dfunc_S(S_binleft,R_binleft);
dSS_true = dSS_true.*S_binleft;
error_dSS = norm(dSS_true-dSS_estimated)./norm(dSS_true);

figure(13)
plot3(S_binleft,R_binleft,dSS_estimated,'r.')
hold on
plot3(S_binleft,R_binleft,dSS_true,'b.')
xlabel('S-Population Size')
ylabel('R-Population Size')
zlabel('Rate')
legend('inferred','true')

figure(14)
plot(S_binleft,dSS_estimated,'r.')
hold on
plot(S_binleft,dSS_true,'b.')
legend('inferred','true')


%% bRR
bRR_true = bfunc_R(S_binleft,R_binleft);
bRR_true = bRR_true.*R_binleft;
error_bRR = norm(bRR_true-bRR_estimated)./norm(bRR_true);

figure(15)
plot3(S_binleft,R_binleft,bRR_estimated,'r.')
hold on
plot3(S_binleft,R_binleft,bRR_true,'b.')
xlabel('S-Population Size')
ylabel('R-Population Size')
zlabel('Rate')
legend('inferred','true')

figure(16)
plot(R_binleft,bRR_estimated,'r.')
hold on
plot(R_binleft,bRR_true,'b.')
legend('inferred','true')


%% dRR
dRR_true = dfunc_R(S_binleft,R_binleft);
dRR_true = dRR_true.*R_binleft;
error_dRR = norm(dRR_true-dRR_estimated)./norm(dRR_true);

figure(17)
plot3(S_binleft,R_binleft,dRR_estimated,'r.')
hold on
plot3(S_binleft,R_binleft,dRR_true,'b.')
xlabel('S-Population Size')
ylabel('R-Population Size')
zlabel('Rate')
legend('inferred','true')

figure(18)
plot(R_binleft,dRR_estimated,'r.')
hold on
plot(R_binleft,dRR_true,'b.')
legend('inferred','true')

%% Inference of inter-species interaction parameter alpha (INCORPORATE HEYRIM'S CODE)

%% Inference of inter-species density dependence parameter sigma
sigmaS_estimated = (b0S_estimated.*S_binleft - gammaS*(rS/KS).*(S_binleft.^2) - bSS_estimated)./(alphaS*(rS/KS).*R_binleft.*S_binleft);
sigmaR_estimated = (b0R_estimated.*R_binleft - gammaR*(rR/KR).*(R_binleft.^2) - bRR_estimated)./(alphaR*(rR/KR).*S_binleft.*R_binleft);

function lhsS = momenteqS(unknownparavecS,rS_estimated,KS_estimated,S_binleft,dSvarvec,dSmoment3cen,dt)
    b0S_unknown    = unknownparavecS(1);
    gammaS_unknown = unknownparavecS(2);
    bS = b0S_unknown - gammaS_unknown*(rS_estimated/KS_estimated).*S_binleft;
    dS = b0S_unknown - rS_estimated +(1-gammaS_unknown)*(rS_estimated/KS_estimated).*S_binleft;
    lhsS = [(bS + dS).*S_binleft*dt - (bS.^2+dS.^2).*S_binleft*(dt^2) - dSvarvec;
            S_binleft*dt.*(1-3*dt*(bS+dS)).*(bS-dS) - dSmoment3cen];
end

function lhsR = momenteqR(unknownparavecR,rR_estimated,KR_estimated,R_binleft,dRvarvec,dRmoment3cen,dt)
    b0R_unknown    = unknownparavecR(1);
    gammaR_unknown = unknownparavecR(2);
    bR = b0R_unknown - gammaR_unknown*(rR_estimated/KR_estimated).*R_binleft;
    dR = b0R_unknown - rR_estimated +(1-gammaR_unknown)*(rR_estimated/KR_estimated).*R_binleft;
    lhsR = [(bR + dR).*R_binleft*dt - (bR.^2+dR.^2).*R_binleft*(dt^2) - dRvarvec;
            R_binleft*dt.*(1-3*dt*(bR+dR)).*(bR-dR) - dRmoment3cen];
end