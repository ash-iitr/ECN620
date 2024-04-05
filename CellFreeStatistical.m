clc;
close all;
clear;

Kusers = [10 20]; % Number of users
Map = 20:20:200; % Number of APs
rho_p_cf = 1;
rho_d_cf = 1; % Allocating all power to downlink signals
iters =1; % for statistical achievable rate
for ku = 1:length(Kusers)
    K = Kusers(ku);
    tau_cf = K;
Rs_k = zeros(K,length(Map));
for mm = 1:length(Map)
M = Map(mm);
Beta = ones(M,K); %largescale path loss
%% Pilot generation

pilots = sqrt(1/2)*(eye(K) + 1i * eye(K));


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

Gamma = sqrt(rho_p_cf*tau_cf)* Beta.*C;
Eta = zeros(M,K);
for m = 1:M  
    for k=1:K
        Eta(m,k) = 1/(K*Gamma(m,k));
    end    
end

%% Equation 24 
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


    Rs_k(kth,mm) = log2(1+ (rho_d_cf*(num_s^2))/(rho_d_cf*den_s1 + rho_d_cf*den_s2 +1));
end

end %end of loop for M antennas at AP
Cap = sum(Rs_k,1)/K;
Cap2(ku,:) = sum(Cap,1)/iters;
end

hold on
for p=1:length(Kusers)
    plot(Map,Cap2(p,:),'-.','LineWidth',2,'DisplayName',sprintf('Statistical - K = %d', Kusers(p)));
    legend('-DynamicLegend');
    legend('show');
end
hold off
xlabel('Number of APs');
ylabel('Achievable Rate Per User (bits/s/Hz)');
