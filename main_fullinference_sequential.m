% close all; clear all;
set(0, 'DefaultLineLineWidth', 1);
Tstart = tic; 
% Linh Huynh, May 15, 2024
% Heyrim Cho, May 16, 2024
% Sequential Inference

%% Parameters  PC3 cell line 
rS      = 0.293;
% rS = b0S - d0S; b0S/d0S = 1/18.3 / 0.015; 
% b0S * 18.3 * 0.015; 
b0S     = rS/(1 - 18.3*0.015); 
d0S     = b0S - rS; 
KS      = 0.843*10^3; 
alphaS  = 0.2; 
gammaS  = 0.0; 
sigmaS  = 0.0; 

rR      = 0.363; 
b0R     = rR/(1 - 16.9*0.015); 
d0R     = b0R - rR; 
KR      = 2.217*10^3; 
alphaR  = 0; 
gammaR  = 0.0; 
sigmaR  = 0.0; 


%% Parameters DU145 
% rS      = 0.306;
% b0S     = rS/(1 - 18.3*0.015);  
% d0S     = b0S - rS; 
% KS      = 0.724*10^3;
% alphaS = 0.25;
% gammaS = 0.5;
% sigmaS  = 0.5;
% 
% rR      = 0.210; 
% b0R     = rR/(1 - 16.9*0.015);
% d0R     = b0R - rR;
% KR      = 1.388*10^3;
% alphaR = -0.5;
% gammaR = 0.5;
% sigmaR  = 0.5;


%% Functions
bfunc_S = @(S,R) b0S - gammaS*(rS/KS).*S - alphaS*sigmaS*(rS/KS).*R;
dfunc_S = @(S,R) d0S + (1-gammaS)*(rS/KS).*S + alphaS*(1-sigmaS)*(rS/KS).*R;
bfunc_R = @(S,R) b0R - gammaR*(rR/KR).*R - alphaR*sigmaR*(rR/KR).*S;
dfunc_R = @(S,R) d0R + (1-gammaR)*(rR/KR).*R + alphaR*(1-sigmaR)*(rR/KR).*S;


%% stochastic simulation parameters 
nrun = 100;
tmax = 100;
dt   = 1/10;
global Gtheta
Gtheta = 10; 
binsize = 10; 

nplot = 0; 


%% Part 1: Mono-Cultured Data
%% Simulation data for mono-culture S 
disp( 'Part 1: Mono-Cultured Data S')
S0 = 1;
R0 = 0;
dtsimulate = 0.1; 
[Smat,Rmat,t] = datasimulation_Langevin(bfunc_S,dfunc_S,bfunc_R,dfunc_R,nrun,tmax,dtsimulate,S0,R0);
[Smat,Rmat] = reshape_data_dt(Smat,Rmat,t,dt,tmax); 

% plot
if( nplot )
for traj = 1:nrun
    figure(1)
    plot(t,Smat(traj,:))
    hold on
    plot(t,Rmat(traj,:))
    xlabel('Time')
    ylabel('Cell Number')
end
hold off 
end 

%% Inference of birth and death rates from mono-cultured data S 
Sbinsize = binsize;
[bSS_estimated,dSS_estimated,binedgesS,IDS,dSvarvec,dSmoment3cen] = separatebirthdeathrates_1D(Smat,dt,Sbinsize);
S_binleft = (binedgesS(IDS))';


%% Paramter estimation for S [rS, KS, gammaS, b0S] - minimize SSQ norm 
param_init = [ 1 5000 0.25 1 ];  
ub = [ 10 1e4  1  10 ]; 
lb = [ 0   1   0   0 ]; 

fitdata.b_estimated = bSS_estimated; 
fitdata.d_estimated = dSS_estimated; 
fitdata.xdata = S_binleft; 



opt = optimoptions('fmincon' ,'Display','none');
[param_fit,ss01] = fmincon(@(params)ssq(params, fitdata),param_init,[],[],[],[],lb,ub,[],opt);

% check results 
if( nplot )

figure(12); 

subplot( 1, 4, 1 ); hold on; plot( fitdata.xdata, bSS_estimated, 'r.' );  xlabel('S');  ylabel('Birth rate')
subplot( 1, 4, 2 ); hold on; plot( fitdata.xdata, dSS_estimated, 'r.' );  xlabel('S');  ylabel('Death rate')

fitdata.xdata = [min(S_binleft):10:max(S_binleft)]; 
[bfunc_S_fit, dfunc_S_fit] = mono_model(param_fit, fitdata.xdata);

bSS_true = bfunc_S(fitdata.xdata,zeros(size(fitdata.xdata))); bSS_true = bSS_true.*fitdata.xdata;
dSS_true = dfunc_S(fitdata.xdata,zeros(size(fitdata.xdata))); dSS_true = dSS_true.*fitdata.xdata;

subplot( 1, 4, 1 ); hold on;  plot( fitdata.xdata, bfunc_S_fit, 'b' ); plot( fitdata.xdata, bSS_true, 'k--' ); xlabel('S');  ylabel('Birth rate')
subplot( 1, 4, 2 ); hold on;  plot( fitdata.xdata, dfunc_S_fit, 'b' ); plot( fitdata.xdata, dSS_true, 'k--' ); xlabel('S');  ylabel('Death rate')
end 

disp( 'fitted parameters for [rS, KS, gammaS, b0S]')
disp( param_fit )

param_S_true = [ rS, KS, gammaS, b0S ]; 
disp( 'relative error : ')
err(1:4) = abs(param_fit - param_S_true)./abs(param_S_true); err(3) = abs(param_fit(3) - param_S_true(3)); 
disp( err(1:4) )

% save param_fit for Part 2. 
param_fit_S = param_fit; 

%% Simulation data for mono-culture R 
disp( 'Part 1: Mono-Cultured Data R')
S0 = 0;
R0 = 1;
[Smat,Rmat,t] = datasimulation_Langevin(bfunc_S,dfunc_S,bfunc_R,dfunc_R,nrun,tmax,dtsimulate,S0,R0);
[Smat,Rmat] = reshape_data_dt(Smat,Rmat,t,dt,tmax); 


% plot
if( nplot )
for traj = 1:nrun
    figure(3)
    plot(t,Smat(traj,:))
    hold on
    plot(t,Rmat(traj,:))
    xlabel('Time')
    ylabel('Cell Number')
end
hold off 
end 

%% Inference of birth and death rates from mono-cultured data R 

Rbinsize = binsize;
[bRR_estimated,dRR_estimated,binedgesR,IDR,dRvarvec,dRmoment3cen] = separatebirthdeathrates_1D(Rmat,dt,Rbinsize);
R_binleft = (binedgesR(IDR))';


%% Paramter estimation for R [rR, KR, gammaR, b0R] - minimize SSQ norm 
param_init = [ 1 5000 0.25 1 ];  
ub = [ 10 1e4   1  10 ]; 
lb = [ 0   1    0   0 ]; 

fitdata.b_estimated = bRR_estimated; 
fitdata.d_estimated = dRR_estimated; 
fitdata.xdata = R_binleft; 


opt = optimoptions('fmincon' ,'Display','none');
[param_fit,ss01] = fmincon(@(params)ssq(params, fitdata),param_init,[],[],[],[],lb,ub,[],opt);

% check results 
if( nplot )

figure(12); 

subplot( 1, 4, 3 ); hold on; plot( fitdata.xdata, bRR_estimated, 'r.' );  xlabel('R');  ylabel('Birth rate')
subplot( 1, 4, 4 ); hold on; plot( fitdata.xdata, dRR_estimated, 'r.' );  xlabel('R');  ylabel('Death rate')

fitdata.xdata = [min(R_binleft):10:max(R_binleft)]; 
[bfunc_R_fit, dfunc_R_fit] = mono_model(param_fit, fitdata.xdata);

bRR_true = bfunc_R(zeros(size(fitdata.xdata)),fitdata.xdata); bRR_true = bRR_true.*fitdata.xdata;
dRR_true = dfunc_R(zeros(size(fitdata.xdata)),fitdata.xdata); dRR_true = dRR_true.*fitdata.xdata;

subplot( 1, 4, 3 ); hold on; plot( fitdata.xdata, bfunc_R_fit, 'b' ); plot( fitdata.xdata, bRR_true, 'k--' ); xlabel('R');  ylabel('Birth rate')
subplot( 1, 4, 4 ); hold on; plot( fitdata.xdata, dfunc_R_fit, 'b' ); plot( fitdata.xdata, dRR_true, 'k--' ); xlabel('R');  ylabel('Death rate')
end 

disp( 'fitted parameters for [rR, KR, gammaR, b0R]')
disp( param_fit )

param_R_true = [ rR, KR, gammaR, b0R ]; 
disp( 'relative error : ')
err(5:8) = abs(param_fit - param_R_true)./abs(param_R_true); err(3+4) = abs(param_fit(3) - param_R_true(3)); 
disp( err(5:8) )

% save param_fit for Part 2. 
param_fit_R = param_fit; 


%% Part 2: Mix-Cultured Data
%% Simulation for mix-culture 
disp( 'Part 2: Co-Cultured Data')
S0 = 1;
R0 = 1;
[Smat,Rmat,t] = datasimulation_Langevin(bfunc_S,dfunc_S,bfunc_R,dfunc_R,nrun,tmax,dtsimulate,S0,R0);
[Smat,Rmat] = reshape_data_dt(Smat,Rmat,t,dt,tmax); 

% plot
if( nplot )
figure
for traj = 1:nrun
    plot(t,Smat(traj,:),'b')
    hold on
    plot(t,Rmat(traj,:), 'r')
    xlabel('Time')
    ylabel('Cell Number')
end
hold off 
end 

%% Inference of birth and death rates from mix-cultured data 

Sbinsize = binsize;
Rbinsize = binsize;
[bSS_estimated,dSS_estimated,bRR_estimated,dRR_estimated,binedgesS,binedgesR,IDS,IDR,dSvarvec,dRvarvec,dSmoment3cen,dRmoment3cen] = separatebirthdeathrates_2D(Smat,Rmat,dt,Sbinsize,Rbinsize);

S_binleft = (binedgesS(IDS))';
R_binleft = (binedgesR(IDR))';


%% Paramter estimation for interaction [alphaS, sigmaS, alphaR, sigmaR] - minimize SSQ norm 

param_init = [ 0.1 0.25 0.1 0.25 ];  
ub = [  2  1   2  1 ]; 
lb = [ -2  0  -2  0 ]; 


fitdata.bSS_estimated = bSS_estimated; 
fitdata.dSS_estimated = dSS_estimated; 
fitdata.bRR_estimated = bRR_estimated; 
fitdata.dRR_estimated = dRR_estimated; 
fitdata.S = S_binleft; 
fitdata.R = R_binleft; 

params_fitted = [param_fit_S, param_fit_R]; 

opt = optimoptions('fmincon' ,'Display','none');
[param_fit,ss01] = fmincon(@(params)ssq2(params, fitdata, params_fitted),param_init,[],[],[],[],lb,ub,[],opt);

% check results 
if( nplot )
[fitdata.S,fitdata.R] = ndgrid(  unique( S_binleft ), unique( R_binleft ) ); 

figure; 
[bfunc_S_fit, dfunc_S_fit, bfunc_R_fit, dfunc_R_fit] = mixture_model(param_fit, fitdata.S, fitdata.R, params_fitted); 

bSS_true = bfunc_S(fitdata.S,fitdata.R); bSS_true = bSS_true.*fitdata.S;
dSS_true = dfunc_S(fitdata.S,fitdata.R); dSS_true = dSS_true.*fitdata.S;
bRR_true = bfunc_R(fitdata.S,fitdata.R); bRR_true = bRR_true.*fitdata.R;
dRR_true = dfunc_R(fitdata.S,fitdata.R); dRR_true = dRR_true.*fitdata.R;

subplot( 2, 2, 1 ); hold on; plot3( S_binleft, R_binleft, bSS_estimated, 'r.' ); alpha 0.5; surfo( fitdata.S, fitdata.R, bfunc_S_fit ); colormap gray;  
plot3( fitdata.S(1:10:end,1:10:end), fitdata.R(1:10:end,1:10:end), bSS_true(1:10:end,1:10:end), 'k.' );
subplot( 2, 2, 2 ); hold on; plot3( S_binleft, R_binleft, dSS_estimated, 'r.' ); alpha 0.5; surfo( fitdata.S, fitdata.R, dfunc_S_fit ); colormap gray;  
plot3( fitdata.S(1:10:end,1:10:end), fitdata.R(1:10:end,1:10:end), dSS_true(1:10:end,1:10:end), 'k.' ); 
subplot( 2, 2, 3 ); hold on; plot3( S_binleft, R_binleft, bRR_estimated, 'r.' ); alpha 0.5; surfo( fitdata.S, fitdata.R, bfunc_R_fit ); colormap gray;  
plot3( fitdata.S(1:10:end,1:10:end), fitdata.R(1:10:end,1:10:end), bRR_true(1:10:end,1:10:end), 'k.' ); 
subplot( 2, 2, 4 ); hold on; plot3( S_binleft, R_binleft, dRR_estimated, 'r.' ); alpha 0.5; surfo( fitdata.S, fitdata.R, dfunc_R_fit ); colormap gray;  
plot3( fitdata.S(1:10:end,1:10:end), fitdata.R(1:10:end,1:10:end), dRR_true(1:10:end,1:10:end), 'k.' ); 
end 

ll = {'$b_S(S,R)$', '$d_S(S,R)$', '$b_R(S,R)$', '$d_R(S,R)$'};
for n = 1:4; subplot( 2, 2, n ); 
xlabel('$S$', 'interprete', 'latex');  ylabel('$R$', 'interprete', 'latex'); zlabel(ll{n}, 'interprete', 'latex'); view( [-5, 40] )
box on; grid on; set(gca, 'fontsize', 12); 
end 
disp( 'fitted parameters for [alphaS, sigmaS, alphaR, sigmaR] ')
disp( param_fit )

param_mix_true = [alphaS, sigmaS, alphaR, sigmaR] ; 
disp( 'abs error ')
err(9:12) = abs( param_fit - param_mix_true );  
disp( err(9:12) )


params_fitted = [params_fitted, param_fit]; 


toc(Tstart)

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


function [bfunc_S, dfunc_S] = mono_model(params, S)

b0S     = params(4);
%b0S = 0.5;
rS      = params(1);
KS      = params(2);
gammaS  = params(3);
d0S     = b0S-rS;

bfunc_S = (b0S -   gammaS  *(rS/KS)*S).*S;  
dfunc_S = (d0S + (1-gammaS)*(rS/KS)*S).*S;  
end 


function err = ssq(params, data) 

S = data.xdata; 

[bfunc_S, dfunc_S] = mono_model(params, S); 

err = norm( data.b_estimated(:) - bfunc_S(:) ) + norm( data.d_estimated(:) - dfunc_S(:) ); 

end 


function [bfunc_S, dfunc_S, bfunc_R, dfunc_R] = mixture_model(params, S, R, params_fitted)

b0S     = params_fitted(4);
rS      = params_fitted(1);
KS      = params_fitted(2);
gammaS  = params_fitted(3);
d0S     = b0S-rS;

b0R     = params_fitted(8);
rR      = params_fitted(5);
KR      = params_fitted(6);
gammaR  = params_fitted(7);
d0R     = b0R-rR;

alphaS  = params(1);
sigmaS  = params(2);

alphaR  = params(3);
sigmaR  = params(4);

% Functions
bfunc_S =  (b0S - gammaS*(rS/KS)*S - alphaS*sigmaS*(rS/KS)*R) .*S;
dfunc_S =  (d0S + (1-gammaS)*(rS/KS)*S + alphaS*(1-sigmaS)*(rS/KS)*R) .*S;
bfunc_R =  (b0R - gammaR*(rR/KR)*R - alphaR*sigmaR*(rR/KR)*S) .*R;
dfunc_R =  (d0R + (1-gammaR)*(rR/KR)*R + alphaR*(1-sigmaR)*(rR/KR)*S) .*R;

end 


function err = ssq2(params, data, params_fitted) 

S = data.S; 
R = data.R; 

[bfunc_S, dfunc_S, bfunc_R, dfunc_R] = mixture_model(params, S, R, params_fitted); 

err = norm( data.bSS_estimated(:) - bfunc_S(:), 'fro' ) + norm( data.dSS_estimated(:) - dfunc_S(:), 'fro' ) ... 
    + norm( data.bRR_estimated(:) - bfunc_R(:), 'fro' ) + norm( data.dRR_estimated(:) - dfunc_R(:), 'fro' ); 

end 

function [Smat_,Rmat_] = reshape_data_dt(Smat,Rmat,t,dt,tmax) 


for n = 1:size( Smat, 1 )
    Smat_(n,:) = interp1(t, Smat(n,:), 0:dt:tmax ); 
    Rmat_(n,:) = interp1(t, Rmat(n,:), 0:dt:tmax ); 
end 

end 
