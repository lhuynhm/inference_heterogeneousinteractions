function [bSS_estimated,dSS_estimated,binedgesS,IDS,dSvarvec,dSmoment3cen] = separatebirthdeathrates_1D(Smat,dt,Sbinsize)
%function [binleftindexSRvec] = separatebirthdeathrates_1D(Smat,dt,Sbinsize)
% By Linh Huynh, December 22, 2023
% updated by Heyrim Cho, May 16, 2024 
% Smat, Rmat: each row is a time series trajectory.
% Row i of Smat corresponds to row i of Rmat, as we simulate a pair of
% (S,R) trajectories at a time
% One pair is a realization of the LV stochastic process

    [nrS,ncS] = size(Smat);
    dSmat = Smat(:,2:end) - Smat(:,1:end-1);
    dSvec = reshape(dSmat,[nrS*(ncS-1),1]);
    maxS  = max(Smat,[],"all");
    Sbinpoints = 0:Sbinsize:maxS;
    
    % the syntax 'discretize' gives the indices of left points
    binleftindexSmat = [];
    for r = 1:nrS
        Straj = Smat(r,:);
        [binleftindexSmat(r,:),binedgesS] = discretize(Straj(1:end-1),Sbinpoints);
    end
    binleftindexSvec = reshape(binleftindexSmat,[nrS*(ncS-1),1]);
    % binleftindexSRvec = [];


    % for r = 1:nrS*(ncS-1)
    %     idS = binleftindexSvec(r,1);
    %     idR = binleftindexRvec(r,1);
    %     idSR = str2double(strcat(num2str(idS),num2str(idR)));
    %     binleftindexSRvec(r,1) = idSR;
    % end

    [G,IDS] = findgroups(binleftindexSvec);
    %indices = [IDS,IDR];
    %indices output gives the grid points of the estimated values for birth and death;
%any grid point which is not included in the list of indices did not have
%any data points in it.
%

    dSmeanvec = splitapply(@mean,dSvec,G);
    dSvarvec  = splitapply(@var,dSvec,G);

    thirdcentral = @(x)moment(x,3);
    dSmoment3cen  = splitapply(thirdcentral,dSvec,G);


    % Count number of samples in each group G
    for n = min(G):max(G) 
        numG(n) = length( find( G==n ) ); 
    end 
    % exclude too small groups to improve accuracy 
    global Gtheta 
    ind = find( numG > Gtheta );
    IDS = IDS(ind); 

    dSmeanvec = dSmeanvec(ind); 
    dSvarvec  = dSvarvec(ind); 
    dSmoment3cen = dSmoment3cen(ind); 
    

    bSS_estimated = (dSvarvec+dSmeanvec)./(2*dt);
    dSS_estimated = (dSvarvec-dSmeanvec)./(2*dt);  

    % shift to mid bin 
    binedgesS = binedgesS + (binedgesS(2)-binedgesS(1))/2; 

end



