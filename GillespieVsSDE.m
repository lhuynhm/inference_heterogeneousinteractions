% close all; clear all;
set(0, 'DefaultLineLineWidth', 1);
% Linh Huynh, May 15, 2024
% Heyrim Cho, May 16, 2024
% Sequential Inference

% %% Parameters
% b0S     = 1.1/120;
% rS      = 1/120;
% KS      = 1e3;
% gammaS  = 0.7;
% d0S     = b0S-rS;
% alphaS  = 0.1;
% sigmaS  = 0.5;
% 
% b0R     = 1.1/120;
% rR      = 1/120;
% KR      = 1e3;
% gammaR = 0.5;
% d0R     = b0R-rR;
% alphaR  = 0.1;
% sigmaR  = 0.5;

%% Parameters PLOS compbio 
% b0S     = 2.5;
% d0S     = 0.5;
% rS      = b0S-d0S;
% KS      = 1e3;
% gammaS  = 0.5;
% alphaS  = 1;
% sigmaS  = 0.5;
% 
% b0R     = 2.25;
% d0R     = 0.5;
% rR      = b0R-d0R;
% KR      = 1e3;
% gammaR = 0.5;
% d0R     = b0R-rR;
% alphaR  = -.1;
% sigmaR  = 0.5;

%% Parameters  PC3 
% rS      = 0.293;
% % rS = b0S - d0S; b0S/d0S = 1/18.3 / 0.015; 
% % b0S * 18.3 * 0.015; 
% b0S     = rS/(1 - 18.3*0.015); 
% d0S     = b0S - rS; 
% KS      = 0.843*10^3; 
% alphaS  = 0; 
% gammaS  = 1; 
% sigmaS  = 1; 
% 
% rR      = 0.363; 
% b0R     = rR/(1 - 16.9*0.015); 
% d0R     = b0R - rR; 
% KR      = 2.217*10^3; 
% alphaR  = 0.2; 
% gammaR  = 1; 
% sigmaR  = 1; 


%% Parameters DU145 
rS      = 0.306;
b0S     = rS/(1 - 18.3*0.015);  
d0S     = b0S - rS; 
KS      = 0.724*10^3;
alphaS = -.5;
gammaS = 1;
sigmaS  = 0.5;

rR      = 0.210; 
b0R     = rR/(1 - 16.9*0.015);
d0R     = b0R - rR;
KR      = 1.388*10^3;
alphaR = .25;
gammaR = 1;
sigmaR  = 0.5;

%% Parameters talk example 
% rS      = 0.306;
% b0S     = rS/(1 - 18.3*0.015);  
% d0S     = b0S - rS; 
% KS      = 0.724*10^3;
% alphaS = .5;
% gammaS = 1;
% sigmaS  = 1;
% 
% rR      = 0.210; 
% b0R     = rR/(1 - 16.9*0.015);
% d0R     = b0R - rR;
% KR      = 1.388*10^3;
% alphaR = .5;
% gammaR = 1;
% sigmaR  = 1;

%% Functions
bfunc_S = @(S,R) max(b0S - gammaS*(rS/KS).*S - alphaS*sigmaS*(rS/KS).*R,0);
dfunc_S = @(S,R) max(d0S + (1-gammaS)*(rS/KS).*S + alphaS*(1-sigmaS)*(rS/KS).*R,0);
bfunc_R = @(S,R) max(b0R - gammaR*(rR/KR).*R - alphaR*sigmaR*(rR/KR).*S,0);
dfunc_R = @(S,R) max(d0R + (1-gammaR)*(rR/KR).*R + alphaR*(1-sigmaR)*(rR/KR).*S,0);

%x = 0:.01:1;
%plot(x,d0S + (1-x)*(rS/KS)*723 + alphaS*(1-0)*(rS/KS)*1000);

%disp(b0S - gammaS*(rS/KS)*1000-alphaS*sigmaS*(rS/KS).*100-(d0S + (1-gammaS)*(rS/KS)*997.9334 + alphaS*(1-sigmaS)*(rS/KS)*1000))

%% stochastic simulation parameters 
nrun = 1000;
tmax = 50;
dt   = 1/100;
nplot = true;

clear NGil
%% Gillespie Simulation Comparison 
S0 = KS-1;
R0 = 1;

[Smat,Rmat,t,dt] = datasimulation_Langevin(bfunc_S,dfunc_S,bfunc_R,dfunc_R,nrun,tmax,dt,S0,R0);
[tGil, NGil] = LV_gillespie_multiple(b0S,rS,KS,gammaS,d0S,alphaS,sigmaS,b0R,rR,KR,gammaR,d0R,alphaR,sigmaR,transpose([S0,R0]),tmax,nrun);


% plot
colororder([0 0 1 ; 1 0 0]);
%colororder([0 0 1]); %blue only
%colororder([1 0 0]); %red only

SDEdie = 0;
Gildie = 0;

for traj = 1:nrun
    %figure(1)
    ind = find(tGil(traj,:)); %nonzero indices of tGil

    %plot(tGil(traj,ind),NGil(2*traj-1,ind)) %sensitive cells
    %plot(t,Smat(traj,:))

    if Rmat(traj,end) == 0 
        SDEdie = SDEdie+1;
    end
    if NGil(2*traj,ind(end)) ==0
        Gildie = Gildie + 1;
    end

    plot(t,Rmat(traj,:))
    hold on
    plot(tGil(traj,ind),NGil(2*traj,ind)) %resistant cells
    %plot(t,Rmat(traj,:))
    %legend('S(t)','R(t)')
    xlabel('Time')
    ylabel('Cell Number')
    title('Gillespie (red) vs SDE (blue)')
end
hold off 
disp(['SDE extinct: ',num2str(SDEdie), ' and Gillespie extinct: ', num2str(Gildie)]);
disp(['extinction prob difference (SDE-Gillespie):',num2str((SDEdie-Gildie)/nrun)]);
