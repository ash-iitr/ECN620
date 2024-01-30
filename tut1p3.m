%% Plot the ergodic capacity of single-input multiple-output (SIMO) channel

NT = 1; % Number of transmit antennas
NR = 5; % Number of receive antennas
N = 1000; % Number of iterations
SNR_dB = -10:2:40; % SNR in dB scale
SNR_linear = 10.^(SNR_dB/10);
% Create a variable to store capacity values
C = zeros(length(SNR_dB),1000);
% Calculate capacity for each SNR and each iteration
for j =1:N
h(:,j) = 1/sqrt(2) .* (randn(NR,1) + 1i*randn(NR,1));
end

for i=1:length(SNR_dB)
    for j=1:N
    C(i,j) = log2(1 + ((norm(h(:,j))^2))* SNR_linear(i) );
    end
end

CE = sum(C,2)/N;
% Plot Capacity vs SNR in dB

plot(SNR_dB,CE); hold on;
xlabel('SNRdB');
ylabel('Ergodic capacity');
