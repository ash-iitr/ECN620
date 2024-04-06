clc;
close all;
clear;

Kusers = [10 20]; % Number of users
Map = 20:20:200; % Number of APs
rho_p_cf = 1;
%tau_cf = K;
%Beta = ones(M,K); %largescale path loss
rho_d_cf = 1; % Allocating all power to downlink signals
iters =1; % for statistical achievable rate
for ku = 1:length(Kusers)
    K = Kusers(ku);
    tau_cf = K;
for it = 1:iters %for statistical achievable rate
Rs_k = zeros(K,length(Map));
for mm = 1:length(Map)
M = Map(mm);
Beta = ones(M,K); %largescale path loss
%% Pilot generation

% % Generate a random Hermitian matrix
% A = randn(K) + 1i * randn(K); % Complex random matrix
% A = (A + conj(A')) / 2; % Ensuring Hermitian property
% 
% % Compute eigenvectors and eigenvalues
% [pilots, ~] = eig(A);

pilots = sqrt(1/2)*(eye(K) + 1i * eye(K));

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


%% Equation 24 - pilots and C needed 
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
    % %% Equation 26 begins
    % num_e = 0;
    % for men = 1:M
    %     num_e = num_e + sqrt(Eta(men,kth))*G(men,kth)*Gcap(men,kth);
    % end
    % 
    % den_e1 = 0;
    % for ked = 1:K
    %     if(ked~=kth)
    %         den_e_i = 0;
    %         for med = 1:M
    %             den_e_i = den_e_i + sqrt(Eta(men,ked))*G(men,kth)*Gcap(men,ked);
    %         end
    %         den_e1 = den_e1 + (abs(den_e_i))^2;
    %     end
    % end
    % 
    % Re_k(kth,mm) = log2(1+ ((rho_d_cf*abs(num_e)^2)/(rho_d_cf*den_e1 + 1)));
    itee = 10;
    for ite = 1:itee % for equation 26
        %% Generate channel coefficients MxK
        G = sqrt(1/2) * (randn(M,K) + 1i * randn(M,K)); % Complex random matrix

        % large scale pathloss Beta = 1 in each case.

        Wp = sqrt(1/2) * (randn(K,M) + 1i * randn(K,M)); %Noise while receiving pilot vectors at AP

        Yp = zeros(K,M); % Equation 2
        for m1 = 1:M
            temp = zeros(K,1);
            for k1=1:K
                temp = temp + G(m1,k1)*pilots(:,k1);
            end
            Yp(:,m1) = sqrt(rho_p_cf*tau_cf)*temp + Wp(:,m1);
        end

        Yhatp = zeros(M,K); % Equation 3
        for m2 = 1:M
            for k2=1:K
                Yhatp(m2,k2) = pilots(:,k2)'*Yp(:,m2);
            end
        end

        C = zeros(M,K);
        for m3 = 1:M
            for k3=1:K

                temp2 = 0;
                for kk = 1:K
                    temp2 = temp2 + (abs(pilots(:,k3)'*pilots(:,kk)))^2;
                end
                C(m3,k3) = (sqrt(rho_p_cf*tau_cf))/(tau_cf*rho_p_cf*temp2 + 1);

            end
        end

        %Gcap = C.*Yhatp; % Equation 4
        Gcap = zeros(M,K);
        for m4 = 1:M
            for k4=1:K
                Gcap(m4,k4) = C(m4,k4)*Yhatp(m4,k4);
            end
        end


        Gamma = sqrt(rho_p_cf*tau_cf)*C;
        Eta = zeros(M,K);
        for m5 = 1:M
            for k5=1:K
                Eta(m5,k5) = 1/(K*Gamma(m5,k5));
            end
        end
             %% Equation 26 begins
             for kte=1:K % Calculate for each user
                 num_e = 0;
                 for men = 1:M
                     num_e = num_e + sqrt(Eta(men,kte))*G(men,kte)*conj(Gcap(men,kte));
                 end

                 den_e1 = 0;
                 for ked = 1:K
                     if(ked~=kte)
                         den_e_i = 0;
                         for med = 1:M
                             den_e_i = den_e_i + sqrt(Eta(men,ked))*G(men,kte)*conj(Gcap(men,ked));
                         end
                         den_e1 = den_e1 + (abs(den_e_i))^2;
                     end
                 end

                 Re_k(kte) = log2(1+ ((rho_d_cf*((abs(num_e))^2))/(rho_d_cf*den_e1 + 1)));

             end


    end
TempExactCap = sum(Re_k,2)/itee;
TempExactCap2(mm,ku) = sum(TempExactCap,1)/K;

%end % end of loop for kth user
end %end of loop for M antennas at AP
Cap(it,:) = sum(Rs_k,1)/K;
%Capexact(it,:) = sum(Re_k,1)/K;
end
Cap2(ku,:) = sum(Cap,1)/iters; % Using statistical
%Capexact2(ku,:) = sum(Capexact,1)/iters; % Using exact realization
end

plot(Map,Cap2(1,:),'-.');
hold on;
plot(Map,Cap2(2,:),'-.');
hold on;
% plot(Map,Capexact2(1,:));
% hold on;
% plot(Map,Capexact2(2,:));
legend('Statistical - K=10','Statistical - K=20', 'Exact - K=10', 'Exact - K=20');
xlabel('Number of APs');
ylabel('Achievable rate per user');



function [G,Yp,Yhatp,C,Gcap,Gamma,Eta] = generateChannels(M,K)
G = sqrt(1/2) * (randn(M,K) + 1i * randn(M,K)); % Complex random matrix

        % large scale pathloss Beta = 1 in each case.

        Wp = sqrt(1/2) * (randn(K,M) + 1i * randn(K,M)); %Noise while receiving pilot vectors at AP

        Yp = zeros(K,M); % Equation 2
        for m1 = 1:M
            temp = zeros(K,1);
            for k1=1:K
                temp = temp + G(m1,k1)*pilots(:,k1);
            end
            Yp(:,m1) = sqrt(rho_p_cf*tau_cf)*temp + Wp(:,m1);
        end

        Yhatp = zeros(M,K); % Equation 3
        for m2 = 1:M
            for k2=1:K
                Yhatp(m2,k2) = pilots(:,k2)'*Yp(:,m2);
            end
        end

        C = zeros(M,K);
        for m3 = 1:M
            for k3=1:K

                temp2 = 0;
                for kk = 1:K
                    temp2 = temp2 + (abs(pilots(:,k3)'*pilots(:,kk)))^2;
                end
                C(m3,k3) = (sqrt(rho_p_cf*tau_cf))/(tau_cf*rho_p_cf*temp2 + 1);

            end
        end

        %Gcap = C.*Yhatp; % Equation 4
        Gcap = zeros(M,K);
        for m4 = 1:M
            for k4=1:K
                Gcap(m4,k4) = C(m4,k4)*Yhatp(m4,k4);
            end
        end


        Gamma = sqrt(rho_p_cf*tau_cf)*C;
        Eta = zeros(M,K);
        for m5 = 1:M
            for k5=1:K
                Eta(m5,k5) = 1/(K*Gamma(m5,k5));
            end
        end
end
