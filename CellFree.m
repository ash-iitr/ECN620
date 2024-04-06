clc;
close all;
clear;

Kusers = [10 20]; % Number of users
Map = 20:20:200; % Number of APs
rho_p_cf = 1;
rho_d_cf = 1; % Allocating all power to downlink signals

for ku = 1:length(Kusers)
    K = Kusers(ku);
    tau_cf = K;
    Rs_k = zeros(K,length(Map));
    
    for mm = 1:length(Map)
        M = Map(mm);
        Beta = ones(M,K); %largescale path loss

        [pilots,~,~,~,~,~,Gamma,Eta] = generateChannels(M,K,rho_p_cf,tau_cf);

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
        % %% Equation 26 begins

        itee = 100;
        % for ite = 1:itee % for equation 26
        % 
        %     [pilots,G,Yp,Yhatp,C,Gcap,Gamma,Eta] = generateChannels(M,K,rho_p_cf,tau_cf);
        %     %% Equation 26 begins
        %     for kte=1:K % Calculate for each user
        %         num_e = 0;
        %         for men = 1:M
        %             num_e = num_e + sqrt(Eta(men,kte))*G(men,kte)*conj(Gcap(men,kte));
        %         end
        % 
        %         den_e1 = 0;
        %         for ked = 1:K
        %             if(ked~=kte)
        %                 den_e_i = 0;
        %                 for med = 1:M
        %                     den_e_i = den_e_i + sqrt(Eta(men,ked))*G(men,kte)*conj(Gcap(men,ked));
        %                 end
        %                 den_e1 = den_e1 + (abs(den_e_i))^2;
        %             end
        %         end
        % 
        %         Re_k(kte) = log2(1+ ((rho_d_cf*((abs(num_e))^2))/(rho_d_cf*den_e1 + 1)));
        % 
        %     end
        % end
        Cek = zeros(K,1);
       for kth = 1:K %for kth user, iteration will run, then averaged in the end
           Ce = zeros(itee,1);
           for ite = 1:itee
               [~,G,~,~,~,Gcap,~,Eta] = generateChannels(M,K,rho_p_cf,tau_cf);
               num_e = 0;
               for m1=1:M
               num_e = num_e + sqrt(Eta(m1,kth))*G(m1,kth)*conj(Gcap(m1,kth));
               end

                tempK = 0;
               for kk=1:K
                   if kk~=kth
                       tempM = 0;
                       for mk=1:M
                           tempM = tempM + sqrt(Eta(mk,kk))*G(mk,kth)*conj(Gcap(mk,kk));
                       end
                       tempK = tempK + abs(tempM)^2;
                   else
                       disp('skipped');
                   end
               end

               Ce(ite) = log2(1+ (num_e^2)/(tempK +1));
           end
           Cek(kth) = mean(Ce);
       end
    Ctest(ku,mm) = mean(Cek);
    end %end of loop for M antennas at AP
    CapPerUserS(ku,:) = sum(Rs_k,1)/K; % Using statistical
    

end

plot(Map,CapPerUserS(1,:),'-.');
hold on;
plot(Map,CapPerUserS(2,:),'-.');
hold on;
plot(Map,Ctest(1,:));
hold on;
plot(Map,Ctest(2,:));
legend('Statistical - K=10','Statistical - K=20', 'Exact - K=10', 'Exact - K=20');
xlabel('Number of APs');
ylabel('Achievable rate per user');



function [pilots,G,Yp,Yhatp,C,Gcap,Gamma,Eta] = generateChannels(M,K,rho_p_cf,tau_cf)

%pilot generation
A = sqrt(1/2)*(rand(K,K) + 1j* rand(K,K));
pilots = orth(A);

%channel generation
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
