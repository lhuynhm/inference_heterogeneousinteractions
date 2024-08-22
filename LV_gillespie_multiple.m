function [tmat,N] = LV_gillespie_multiple(b0S,rS,KS,gammaS,d0S,alphaS,sigmaS,b0R,rR,KR,gammaR,d0R,alphaR,sigmaR,N0,t_f,nruns)
%Output variables: 
%each row of tmat is a separate run; columns show the times of each event 
%odd rows, 2*j-1, of N are the counts of sensitive cells at each of the time points
%in the corresponding column of tmat(j,:)
%even rows, 2*j, of N are the counts of resistant cells at each of the time
%points in the corresponding column of tmat(j,:)
%odd row 2*j-1 and even row 2*j were run in the same simulation (sensitive
%and resistant cells, respectively)

% Define stoichiometric matrix for the GV system:
v_1b = [1; 0];
v_1d = [-1; 0];
v_2b = [0; 1];
v_2d = [0; -1];
V = [v_1b v_1d v_2b v_2d];

% Parameter vector
c(1)  = b0S;
c(2)  = rS;
c(3)  = KS;
c(4)  = gammaS;
c(5)  = d0S;
c(6)  = alphaS;
c(7)  = sigmaS;
c(8)  = b0R;
c(9)  = rR;
c(10) = KR;
c(11) = gammaR;
c(12) = d0R;
c(13) = alphaR;
c(14) = sigmaR;

% Initialize simulations
tmat(1,1) = 0;
N = NaN(nruns*2,1);

for j = 1:nruns
    N(2*j-1,1) = N0(1,1);
    N(2*j,1) = N0(2,1);
    tmat(j,1) = 0;
    i=1;
    while tmat(j,i) < t_f
        i = i+1;
        % Calculate reaction hazards:
        % birth of S-cells
        a(1) = c(1)*N(2*j-1,i-1) - c(4)*(c(2)/c(3))*N(2*j-1,i-1)*N(2*j-1,i-1) - c(7)*c(6)*(c(2)/c(3))*N(2*j-1,i-1)*N(2*j,i-1);
        % death of S-cells
        a(2) = c(5)*N(2*j-1,i-1) + (1-c(4))*(c(2)/c(3))*N(2*j-1,i-1)*N(2*j-1,i-1) + (1 - c(7))*c(6)*(c(2)/c(3))*N(2*j-1,i-1)*N(2*j,i-1);
        % birth of R cells
        a(3) = c(8)*N(2*j,i-1) - c(11)*(c(9)/c(10))*N(2*j,i-1)*N(2*j,i-1) - c(14)*c(13)*(c(9)/c(10))*N(2*j,i-1)*N(2*j-1,i-1);
        % death of R cells
        a(4) = c(12)*N(2*j,i-1) + (1-c(11))*(c(9)/c(10))*N(2*j,i-1)*N(2*j,i-1) + (1-c(14))*c(13)*(c(9)/c(10))*N(2*j,i-1)*N(2*j-1,i-1);
        
        % Draw time of next reaction from exponential distribution:
        a(a<0)=0; 
        asum = sum(a);
        tmat(j,i) = tmat(j,i-1)+log(1/rand)/asum;
        % Convert reaction hazards into reaction probabilities and draw
        % reaction to occur:
        n = rand; k = 1; Phi = a(1)/asum;
        while Phi < n
            k = k+1;
            Phi = Phi+a(k)/asum;
        end
        % Update state vector given drawn reaction:
        N(2*j-1,i) = N(2*j-1,i-1)+V(1,k);
        N(2*j,i) = N(2*j,i-1)+V(2,k);
    end
end
