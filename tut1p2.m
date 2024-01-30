%% Plot the capacity versus SNR for a fixed single-input multiple-output (SIMO) channel. Assume Nr=5 receive antennas. 

NR = 5;
SNRdB = -10:2:40;
SNR_linear = 10.^(SNRdB/10);
% Create a variable to store the values of Capacity
C = zeros(length(SNRdB),1);
% Compute capacity for each SNR
h = 1/sqrt(2) .* (randn(NR,1) + 1i*randn(NR,1));

for i=1:length(SNRdB)
    C(i) = log2(1 + ((norm(h)^2))* SNR_linear(i) );
end

% Plot Capacity vs SNR in dB
plot(SNRdB,C); hold on;
xlabel('SNRdB');
ylabel('Capacity');
