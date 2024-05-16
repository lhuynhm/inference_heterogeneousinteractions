% Linh Huynh, April 2022
function [Smat,Rmat,t,dt] = datasimulation_Langevine(bfunc_S,dfunc_S,bfunc_R,dfunc_R,nrun,tmax,dt,S0,R0)
    tic   
    %tmax = 3000; 
    %dt   = 1/30; %logistic growth
    t    = 0:dt:tmax;

    Smat=nan(nrun,length(t)); % State vector preallocation
    dWS=randn(size(Smat))*sqrt(dt); % Wiener process increments

    Rmat=nan(nrun,length(t)); % State vector preallocation
    dWR=randn(size(Rmat))*sqrt(dt); % Wiener process increments
    
    
    %X(:,1)=poissrnd(3); % logistic growth
    Smat(:,1) = S0;
    Rmat(:,1) = R0;
    
    for j=2:length(t)
        Smat(:,j)=Smat(:,j-1)+...
            dt*Smat(:,j-1).*(bfunc_S(Smat(:,j-1),Rmat(:,j-1))-dfunc_S(Smat(:,j-1),Rmat(:,j-1)))+...
            sqrt(Smat(:,j-1).*(bfunc_S(Smat(:,j-1),Rmat(:,j-1))+dfunc_S(Smat(:,j-1),Rmat(:,j-1)))).*dWS(:,j);
        Smat(:,j) = max(Smat(:,j),0);

        Rmat(:,j)=Rmat(:,j-1)+...
            dt*Rmat(:,j-1).*(bfunc_R(Smat(:,j-1),Rmat(:,j-1))-dfunc_R(Smat(:,j-1),Rmat(:,j-1)))+...
            sqrt(Rmat(:,j-1).*(bfunc_R(Smat(:,j-1),Rmat(:,j-1))+dfunc_R(Smat(:,j-1),Rmat(:,j-1)))).*dWR(:,j);
        Rmat(:,j) = max(Rmat(:,j),0);
    end

    % Ought to protect against X<0 but probably won't happen for these params.
    runtime=toc;
    disp('Run time')
    disp(runtime)

    %% Store
    %save datasimulation_Langevine.mat
end
