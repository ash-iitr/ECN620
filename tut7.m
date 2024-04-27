rng(111);
N=10^3; % Number of channel realizations
Nt=1;   % No. of transmit antenna per user
K=16    % No. of users
Nr=64;  % No. of receive antenna                                                              
Nrf= 16; % No. of RF chains at the receiver
B=0.5*10^6;% Bandwidth in Hz

%% OVERALL SYSTEM
% SNR range
EsNo=-20:2.5:20;
Average_SR = zeros(size(EsNo));
% Loop over the SNR values
for iter_1=1:length(EsNo)
Sum_rate=0;
Sum_add = 0;

% Loop over the channel realizations
for iter_2=1:N
% Generate Rayleigh faded channel matrix
H= (randn(Nr, K) + 1j * randn(Nr, K)) / sqrt(2); 
%% Compute SVD of the generated channel matrix
[U, S, V] = svd(H);
%% Sum rate calculation for each user
        Sum_add = zeros(1, K);

Noise_pow = 1/(10^(EsNo(iter_1)/10)); % Compute Noise power
        
for l = 1 : K
   Hl = H(:, l);
   Ul = U(:, l);
            
   Sig_pow = abs(Ul'*Hl)^2; % Compute signal power
 
   % Interference is calculated by considering all columns of matrix H except lth column as it will provide desired signal
   Inter = sum(abs(Ul' * H(:, [1:l-1, l+1:end])).^2);
   

% Compute the sum rate for each user
%%Sum_add(iter_2) = B * log2(1 + Sig_pow / (Inter + Noise_pow));
Sum_add(l) = B * log2(1 + Sig_pow / (Inter + Noise_pow));
end
        
% Compute the sum rate
Sum_rate = Sum_rate + sum(Sum_add);

end
%% Compute the average sum rate per user
Average_SR(iter_1) = Sum_rate/(N*K);
end
figure
plot(EsNo,Average_SR,'r') ;
xlabel('Es/N0');
ylabel('Sum Rate (bits per sec)');
legend('SVD based digital beamcombiner');
title('Average Sum Rate');
