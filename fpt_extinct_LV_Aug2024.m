% Linh Huynh, August 2024

close all; clear all;
bS  = 2/120;
dS  = 1/120;
rS  = bS - dS;
KS  = 1e2; %carrying capacity of S-cells
bR  = 2/120;
dR  = 1/120;
rR  = bR - dR;
KR  = 1e2; %carrying capacity of S-cells
gSR = 0;
mu  = 0;
alphaS = -0.5; %to what extent interactions with R-cells cause death of S-cells
alphaR = -0.5; %to what extent interactions with S-cells cause death of R-cells
% max numbers of cells
maxS = 20; % total number of S-cells
maxR = 20; % total number of R-cells
dim  = (maxS+1)*(maxR+1);

% universal plot
cmax = 2*1e2; 

%% Case 1: gammaS = 0; gammaR = 0; sigmaS = 0; sigmaR = 0
% parameter
gammaS = 0; %density dependence within S-cells
gammaR = 0; %density dependence within R-cells
sigmaS = 0;
sigmaR = 0;

% numbers of cells vectors
numbS     = (0:1:maxS)';
numbS_vec = [];
numbR_vec = [];
for j = 0:maxR
    numbS_vec = [numbS_vec;numbS];
    numbR     = repmat(j,maxS+1,1);
    numbR_vec = [numbR_vec;numbR];
end

% generator matrix
indices = (1:1:dim)';
nS_vec = [];
nR_vec = [];
for k = 2:length(indices)
    i  = indices(k);    
    nS = mod(mod(i,maxS+1)-1,maxS+1); %number of resistant cells corresponding to this index
    nS_vec = [nS_vec;nS];
    %nR = floor(i/(maxS+1)); %number of resistant cells corresponding to this index
    nR = numbR_vec(ceil(i/(maxS+1))*(maxS+1));
    nR_vec = [nR_vec;nR];
    %fprintf('New iteration\n')
    % reaction rates
    r1 = max(dR*nR + (1-gammaR)*(rR/KR)*nR*nR + (1-sigmaR)*alphaR*(rR/KR)*nS*nR,0); %death of R-cells
    r3 = max(dS*nS + (1-gammaS)*(rS/KS)*nS*nS + (1-sigmaS)*alphaS*(rS/KS)*nR*nS,0); %death of S-cells
    r5 = max(bS*nS - gammaS*(rS/KS)*nS*nS - sigmaS*alphaS*(rS/KS)*nR*nS,0); %birth of S-cells
    r6 = max(bR*nR - gammaR*(rR/KR)*nR*nR - sigmaR*alphaR*(rR/KR)*nS*nR,0); %birth of R-cells 
    r4 = -(r1+r3+r5+r6); %stay the same
    rates  = [r1;r3;r4;r5;r6];
    % column indices corresponding to the reactions
    j1 = i-(maxS+1);  %death of resistant cells
    %j2 = i+maxS;      % transition from sensitive to resistant (-1,1)
    j3 = i-1; %death of sensitive cells
    j4 = i;   %stay the same
    j5 = i+1; %birth of sensitive cells
    j6 = i+(maxS+1);
    idcoli = [j1;j3;j4;j5;j6];
    % test interior
    test1 = (numbR_vec(i)-1>=0); % death of resistant cells
    %test2 = (numbS_vec(i)-1>=0)&&(numbR_vec(i)+1<=maxR); % transition from sensitive to resistant
    test3 = (numbS_vec(i)-1>=0); % death of sensitive cells
    test4 = 1; % stay at current state
    test5 = (numbS_vec(i)+1<=maxS); % birth of sensitive cells
    test6 = (numbR_vec(i)+1<=maxR); % birth of resistant cells
    test  = [test1; test3; test4; test5; test6];
    for idtest = 1:length(test)
        if test(idtest)==1
            idcoli(idtest);
            Q(i,idcoli(idtest)) = rates(idtest);
        else
            continue
        end
    end
    Q(i,i) = -(sum(Q(i,:)) - Q(i,i));
end
Q(1,1) = 1;
rhs_vec = -ones(dim,1);
rhs_vec(1) = 0;
tau_vec = Q\rhs_vec;
tau_mat = (reshape(tau_vec,[maxS+1,maxR+1]))';

fg1 = figure(1);
imagesc(tau_mat)
set(gca,'FontSize',15)
set(gca,'Ydir','Normal')
hold on
contour(tau_mat,'k','LineWidth',2)
colorbar
xlabel('Number of Type-S Cells')
ylabel('Number of Type-R Cells')
caxis([0,cmax])
AddLetters2Plots(fg1, {'(A)'},'HShift', 0.55,'FontSize',27)
title(['$\gamma_S$=',num2str(gammaS),'$, \gamma_R$=',num2str(gammaR),'$, \sigma_S$=',num2str(sigmaS),'$, \sigma_R$=',num2str(sigmaR)],'Interpreter','latex')

%% Case 2: gammaS = 0.5; gammaR = 0; sigmaS = 0; sigmaR = 0
% parameter
gammaS = 0.5; %density dependence within S-cells
gammaR = 0; %density dependence within R-cells
sigmaS = 0;
sigmaR = 0;

% numbers of cells vectors
numbS     = (0:1:maxS)';
numbS_vec = [];
numbR_vec = [];
for j = 0:maxR
    numbS_vec = [numbS_vec;numbS];
    numbR     = repmat(j,maxS+1,1);
    numbR_vec = [numbR_vec;numbR];
end

% generator matrix
indices = (1:1:dim)';
nS_vec = [];
nR_vec = [];
for k = 2:length(indices)
    i  = indices(k);    
    nS = mod(mod(i,maxS+1)-1,maxS+1); %number of resistant cells corresponding to this index
    nS_vec = [nS_vec;nS];
    %nR = floor(i/(maxS+1)); %number of resistant cells corresponding to this index
    nR = numbR_vec(ceil(i/(maxS+1))*(maxS+1));
    nR_vec = [nR_vec;nR];
    %fprintf('New iteration\n')
    % reaction rates
    r1 = max(dR*nR + (1-gammaR)*(rR/KR)*nR*nR + (1-sigmaR)*alphaR*(rR/KR)*nS*nR,0); %death of R-cells
    r3 = max(dS*nS + (1-gammaS)*(rS/KS)*nS*nS + (1-sigmaS)*alphaS*(rS/KS)*nR*nS,0); %death of S-cells
    r5 = max(bS*nS - gammaS*(rS/KS)*nS*nS - sigmaS*alphaS*(rS/KS)*nR*nS,0); %birth of S-cells
    r6 = max(bR*nR - gammaR*(rR/KR)*nR*nR - sigmaR*alphaR*(rR/KR)*nS*nR,0); %birth of R-cells 
    r4 = -(r1+r3+r5+r6); %stay the same
    rates  = [r1;r3;r4;r5;r6];
    % column indices corresponding to the reactions
    j1 = i-(maxS+1);  %death of resistant cells
    %j2 = i+maxS;      % transition from sensitive to resistant (-1,1)
    j3 = i-1; %death of sensitive cells
    j4 = i;   %stay the same
    j5 = i+1; %birth of sensitive cells
    j6 = i+(maxS+1);
    idcoli = [j1;j3;j4;j5;j6];
    % test interior
    test1 = (numbR_vec(i)-1>=0); % death of resistant cells
    %test2 = (numbS_vec(i)-1>=0)&&(numbR_vec(i)+1<=maxR); % transition from sensitive to resistant
    test3 = (numbS_vec(i)-1>=0); % death of sensitive cells
    test4 = 1; % stay at current state
    test5 = (numbS_vec(i)+1<=maxS); % birth of sensitive cells
    test6 = (numbR_vec(i)+1<=maxR); % birth of resistant cells
    test  = [test1; test3; test4; test5; test6];
    for idtest = 1:length(test)
        if test(idtest)==1
            idcoli(idtest);
            Q(i,idcoli(idtest)) = rates(idtest);
        else
            continue
        end
    end
    Q(i,i) = -(sum(Q(i,:)) - Q(i,i));
end
Q(1,1) = 1;
rhs_vec = -ones(dim,1);
rhs_vec(1) = 0;
tau_vec = Q\rhs_vec;
tau_mat = (reshape(tau_vec,[maxS+1,maxR+1]))';

fg2 = figure(2);
imagesc(tau_mat)
set(gca,'FontSize',15)
set(gca,'Ydir','Normal')
hold on
contour(tau_mat,'k','LineWidth',2)
colorbar
xlabel('Number of Type-S Cells')
ylabel('Number of Type-R Cells')
caxis([0,cmax])
AddLetters2Plots(fg2, {'(B)'},'HShift', 0.55,'FontSize',27)
title(['$\gamma_S$=',num2str(gammaS),'$, \gamma_R$=',num2str(gammaR),'$, \sigma_S$=',num2str(sigmaS),'$, \sigma_R$=',num2str(sigmaR)],'Interpreter','latex')

%% Case 3: gammaS = 0.5; gammaR = 0; sigmaS = 0.5; sigmaR = 0
% parameter
gammaS = 0; %density dependence within S-cells
gammaR = 0; %density dependence within R-cells
sigmaS = 0.5;
sigmaR = 0;

% numbers of cells vectors
numbS     = (0:1:maxS)';
numbS_vec = [];
numbR_vec = [];
for j = 0:maxR
    numbS_vec = [numbS_vec;numbS];
    numbR     = repmat(j,maxS+1,1);
    numbR_vec = [numbR_vec;numbR];
end

% generator matrix
indices = (1:1:dim)';
nS_vec = [];
nR_vec = [];
for k = 2:length(indices)
    i  = indices(k);    
    nS = mod(mod(i,maxS+1)-1,maxS+1); %number of resistant cells corresponding to this index
    nS_vec = [nS_vec;nS];
    %nR = floor(i/(maxS+1)); %number of resistant cells corresponding to this index
    nR = numbR_vec(ceil(i/(maxS+1))*(maxS+1));
    nR_vec = [nR_vec;nR];
    %fprintf('New iteration\n')
    % reaction rates
    r1 = max(dR*nR + (1-gammaR)*(rR/KR)*nR*nR + (1-sigmaR)*alphaR*(rR/KR)*nS*nR,0); %death of R-cells
    r3 = max(dS*nS + (1-gammaS)*(rS/KS)*nS*nS + (1-sigmaS)*alphaS*(rS/KS)*nR*nS,0); %death of S-cells
    r5 = max(bS*nS - gammaS*(rS/KS)*nS*nS - sigmaS*alphaS*(rS/KS)*nR*nS,0); %birth of S-cells
    r6 = max(bR*nR - gammaR*(rR/KR)*nR*nR - sigmaR*alphaR*(rR/KR)*nS*nR,0); %birth of R-cells 
    r4 = -(r1+r3+r5+r6); %stay the same
    rates  = [r1;r3;r4;r5;r6];
    % column indices corresponding to the reactions
    j1 = i-(maxS+1);  %death of resistant cells
    %j2 = i+maxS;      % transition from sensitive to resistant (-1,1)
    j3 = i-1; %death of sensitive cells
    j4 = i;   %stay the same
    j5 = i+1; %birth of sensitive cells
    j6 = i+(maxS+1);
    idcoli = [j1;j3;j4;j5;j6];
    % test interior
    test1 = (numbR_vec(i)-1>=0); % death of resistant cells
    %test2 = (numbS_vec(i)-1>=0)&&(numbR_vec(i)+1<=maxR); % transition from sensitive to resistant
    test3 = (numbS_vec(i)-1>=0); % death of sensitive cells
    test4 = 1; % stay at current state
    test5 = (numbS_vec(i)+1<=maxS); % birth of sensitive cells
    test6 = (numbR_vec(i)+1<=maxR); % birth of resistant cells
    test  = [test1; test3; test4; test5; test6];
    for idtest = 1:length(test)
        if test(idtest)==1
            idcoli(idtest);
            Q(i,idcoli(idtest)) = rates(idtest);
        else
            continue
        end
    end
    Q(i,i) = -(sum(Q(i,:)) - Q(i,i));
end
Q(1,1) = 1;
rhs_vec = -ones(dim,1);
rhs_vec(1) = 0;
tau_vec = Q\rhs_vec;
tau_mat = (reshape(tau_vec,[maxS+1,maxR+1]))';

fg3 = figure(3);
imagesc(tau_mat)
set(gca,'FontSize',15)
set(gca,'Ydir','Normal')
hold on
contour(tau_mat,'k','LineWidth',2)
colorbar
xlabel('Number of Type-S Cells')
ylabel('Number of Type-R Cells')
caxis([0,cmax])
AddLetters2Plots(fg3, {'(C)'},'HShift', 0.55,'FontSize',27)
title(['$\gamma_S$=',num2str(gammaS),'$, \gamma_R$=',num2str(gammaR),'$, \sigma_S$=',num2str(sigmaS),'$, \sigma_R$=',num2str(sigmaR)],'Interpreter','latex')

%% Case 4: gammaS = 1; gammaR = 1; sigmaS = 1; sigmaR = 1
% parameter
gammaS = 1; %density dependence within S-cells
gammaR = 1; %density dependence within R-cells
sigmaS = 1;
sigmaR = 1;

% numbers of cells vectors
numbS     = (0:1:maxS)';
numbS_vec = [];
numbR_vec = [];
for j = 0:maxR
    numbS_vec = [numbS_vec;numbS];
    numbR     = repmat(j,maxS+1,1);
    numbR_vec = [numbR_vec;numbR];
end

% generator matrix
indices = (1:1:dim)';
nS_vec = [];
nR_vec = [];
for k = 2:length(indices)
    i  = indices(k);    
    nS = mod(mod(i,maxS+1)-1,maxS+1); %number of resistant cells corresponding to this index
    nS_vec = [nS_vec;nS];
    %nR = floor(i/(maxS+1)); %number of resistant cells corresponding to this index
    nR = numbR_vec(ceil(i/(maxS+1))*(maxS+1));
    nR_vec = [nR_vec;nR];
    %fprintf('New iteration\n')
    % reaction rates
    r1 = max(dR*nR + (1-gammaR)*(rR/KR)*nR*nR + (1-sigmaR)*alphaR*(rR/KR)*nS*nR,0); %death of R-cells
    r3 = max(dS*nS + (1-gammaS)*(rS/KS)*nS*nS + (1-sigmaS)*alphaS*(rS/KS)*nR*nS,0); %death of S-cells
    r5 = max(bS*nS - gammaS*(rS/KS)*nS*nS - sigmaS*alphaS*(rS/KS)*nR*nS,0); %birth of S-cells
    r6 = max(bR*nR - gammaR*(rR/KR)*nR*nR - sigmaR*alphaR*(rR/KR)*nS*nR,0); %birth of R-cells 
    r4 = -(r1+r3+r5+r6); %stay the same
    rates  = [r1;r3;r4;r5;r6];
    % column indices corresponding to the reactions
    j1 = i-(maxS+1);  %death of resistant cells
    %j2 = i+maxS;      % transition from sensitive to resistant (-1,1)
    j3 = i-1; %death of sensitive cells
    j4 = i;   %stay the same
    j5 = i+1; %birth of sensitive cells
    j6 = i+(maxS+1);
    idcoli = [j1;j3;j4;j5;j6];
    % test interior
    test1 = (numbR_vec(i)-1>=0); % death of resistant cells
    %test2 = (numbS_vec(i)-1>=0)&&(numbR_vec(i)+1<=maxR); % transition from sensitive to resistant
    test3 = (numbS_vec(i)-1>=0); % death of sensitive cells
    test4 = 1; % stay at current state
    test5 = (numbS_vec(i)+1<=maxS); % birth of sensitive cells
    test6 = (numbR_vec(i)+1<=maxR); % birth of resistant cells
    test  = [test1; test3; test4; test5; test6];
    for idtest = 1:length(test)
        if test(idtest)==1
            idcoli(idtest);
            Q(i,idcoli(idtest)) = rates(idtest);
        else
            continue
        end
    end
    Q(i,i) = -(sum(Q(i,:)) - Q(i,i));
end
Q(1,1) = 1;
rhs_vec = -ones(dim,1);
rhs_vec(1) = 0;
tau_vec = Q\rhs_vec;
tau_mat = (reshape(tau_vec,[maxS+1,maxR+1]))';

fg4 = figure(4);
imagesc(tau_mat)
set(gca,'FontSize',15)
set(gca,'Ydir','Normal')
hold on
contour(tau_mat,'k','LineWidth',2)
colorbar
xlabel('Number of Type-S Cells')
ylabel('Number of Type-R Cells')
caxis([0,cmax])
AddLetters2Plots(fg4, {'(D)'},'HShift', 0.55,'FontSize',27)
title(['$\gamma_S$=',num2str(gammaS),'$, \gamma_R$=',num2str(gammaR),'$, \sigma_S$=',num2str(sigmaS),'$, \sigma_R$=',num2str(sigmaR)],'Interpreter','latex')

