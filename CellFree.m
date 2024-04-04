clc;
close all;
clear;

K = 20; % Number of users
Map = 20:20:200; % Number of APs
rho_p_cf = 1;
tau_cf = K;
%Beta = ones(M,K); %largescale path loss
rho_d_cf = 1; % Allocating all power to downlink signals

Rd_k = zeros(K,length(Map));
for mm = 1:length(Map)
M = Map(mm);
Beta = ones(M,K); %largescale path loss
%% Pilot generation

% Generate a random Hermitian matrix
A = randn(K) + 1i * randn(K); % Complex random matrix
A = (A + conj(A')) / 2; % Ensuring Hermitian property

% Compute eigenvectors and eigenvalues
[pilots, ~] = eig(A);

check = norm(pilots(:,5))^2
check2 = abs(pilots(:,1)'*pilots(:,2))^2  % tends to 0 as n increases

% Conclusion - My pilots are orthogonal
% we have now K different Kx1 size orthogonal pilot sequences
% tau cf = K

%% Generate channel coefficients MxK
G = sqrt(1/2) * (randn(M,K) + 1i * randn(M,K)); % Complex random matrix

% large scale pathloss Beta = 1 in each case.

Wp = sqrt(1/2) * (randn(K,M) + 1i * randn(K,M)); %Noise while receiving pilot vectors at AP

Yp = zeros(K,M); % Equation 2
for m = 1:M
    temp = zeros(K,1);
    for k=1:K
        temp = temp + G(m,k)*pilots(:,k);
    end
    Yp(:,m) = sqrt(rho_p_cf*tau_cf)*temp + Wp(:,m);
end

Yhatp = zeros(M,K); % Equation 3
for m = 1:M  
    for k=1:K
        Yhatp(m,k) = pilots(:,k)'*Yp(:,m);
    end    
end

C = zeros(M,K);
for m = 1:M  
    for k=1:K
        temp2 = 0;
        for kk = 1:K
            temp2 = temp2 + Beta(m,kk)*(abs(pilots(:,k)'*pilots(:,kk)))^2;
        end
        C(m,k) = (sqrt(rho_p_cf*tau_cf)*Beta(m,k))/(tau_cf*rho_p_cf*temp2 + 1);
    end    
end

Gcap = C.*Yhatp; % Equation 4
Gamma = sqrt(rho_p_cf*tau_cf)* Beta.*C;
Eta = zeros(M,K);
for m = 1:M  
    for k=1:K
        Eta(m,k) = 1/(K*Gamma(m,k));
    end    
end


%% Equation 24 begins 
for kth = 1:K
    num_s = 0;
    for msn = 1:M
        num_s = num_s + sqrt(Eta(msn,kth))*Gamma(msn,kth);
    end

    den_s1 = 0;
    for ksd = 1:K
        if(ksd~=kth)
            den_s_i = 0;
            for msd = 1:M
                den_s_i = den_s_i + sqrt(Eta(msd,ksd))*Gamma(msd,ksd)*(Beta(msd,kth)/Beta(msd,ksd));
            end
            den_s1 = den_s1 + (den_s_i^2)*(abs(pilots(:,ksd)'*pilots(:,kth)))^2;
        end
    end

    den_s2=0;
    for ksd2 = 1:K
        for msd2=1:M
            den_s2 = den_s2 + Eta(msd2,ksd2)*Gamma(msd2,ksd2)*Beta(msd2,kth);
        end
    end


    Rd_k(kth,mm) = log2(1+ (rho_d_cf*(num_s^2))/(rho_d_cf*den_s1 + rho_d_cf*den_s2 +1));

end % end of loop for kth user
end %end of loop for M antennas at AP
Cap = sum(Rd_k,1)/K;
plot(Map,Cap);
