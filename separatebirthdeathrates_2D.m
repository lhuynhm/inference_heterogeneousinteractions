function [bSS_estimated,dSS_estimated,bRR_estimated,dRR_estimated,binedgesS,binedgesR,IDS,IDR,dSvarvec,dRvarvec,dSmoment3cen,dRmoment3cen] = separatebirthdeathrates_2D(Smat,Rmat,dt,Sbinsize,Rbinsize)
% By Linh Huynh, December 22, 2023
% Smat, Rmat: each row is a time series trajectory.
% Row i of Smat corresponds to row i of Rmat, as we simulate a pair of
% (S,R) trajectories at a time
% One pair is a realization of the LV stochastic process

    [nrS,ncS] = size(Smat);
    [nrR,ncR] = size(Rmat);
    dSmat = Smat(:,2:end) - Smat(:,1:end-1);
    dRmat = Rmat(:,2:end) - Rmat(:,1:end-1);
    dSvec = reshape(dSmat,[nrS*(ncS-1),1]);
    dRvec = reshape(dRmat,[nrR*(ncR-1),1]);
    maxS  = max(Smat,[],"all");
    maxR  = max(Rmat,[],"all");
    Sbinpoints = 0:Sbinsize:maxS;
    Rbinpoints = 0:Rbinsize:maxR;
    
    if nrS ~= nrR
        print('number of S trajectories is not equal to number of R trajectories')
        return
    end
    % the syntax 'discretize' gives the indices of left points
    binleftindexSmat = [];
    binleftindexRmat = [];
    for r = 1:nrS
        Straj = Smat(r,:);
        Rtraj = Rmat(r,:);
        [binleftindexSmat(r,:),binedgesS] = discretize(Straj(1:end-1),Sbinpoints);
        [binleftindexRmat(r,:),binedgesR] = discretize(Rtraj(1:end-1),Rbinpoints);
    end
    binleftindexSvec = reshape(binleftindexSmat,[nrS*(ncS-1),1]);
    binleftindexRvec = reshape(binleftindexRmat,[nrR*(ncR-1),1]);

    [G,IDS,IDR] = findgroups(binleftindexSvec,binleftindexRvec);
    
    dSmeanvec = splitapply(@mean,dSvec,G);
    dSvarvec  = splitapply(@var,dSvec,G);

    dRmeanvec = splitapply(@mean,dRvec,G);
    dRvarvec  = splitapply(@var,dRvec,G);

    thirdcentral = @(x)moment(x,3);
    dSmoment3cen  = splitapply(thirdcentral,dSvec,G);
    dRmoment3cen  = splitapply(thirdcentral,dRvec,G);


    % Count number of samples in each group G
    for n = min(G):max(G) 
        numG(n) = length( find( G==n ) ); 
    end 
    % exclude too small groups to improve accuracy 
    global Gtheta  
    ind = find( numG > Gtheta );
    IDS = IDS(ind); 
    IDR = IDR(ind); 
    dSmeanvec = dSmeanvec(ind); 
    dSvarvec  = dSvarvec(ind); 
    dRmeanvec = dRmeanvec(ind); 
    dRvarvec  = dRvarvec(ind); 
    dSmoment3cen = dSmoment3cen(ind); 
    dRmoment3cen = dRmoment3cen(ind); 
    

    bSS_estimated = (dSvarvec+dSmeanvec)./(2*dt);
    dSS_estimated = (dSvarvec-dSmeanvec)./(2*dt);  

    bRR_estimated = (dRvarvec+dRmeanvec)./(2*dt);
    dRR_estimated = (dRvarvec-dRmeanvec)./(2*dt);  

    % shift to mid bin 
    binedgesS = binedgesS + (binedgesS(2)-binedgesS(1))/2; 
    binedgesR = binedgesR + (binedgesR(2)-binedgesR(1))/2; 


end



