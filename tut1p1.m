%%Plot the capacity versus SNR for AWGN single-input single-output (SISO) channel. C= log2(1+ SNR)

SNRdB = -10:2:40; %Range of SNR in dB scale

% Linear Conversion
SNR_linear = 10.^(SNRdB/10);
% Compute Capacity for SISO channel
for i=1:length(SNR_linear)
    C(i) = log2(1 + SNR_linear(i));
end
% Plot capacity vs SNR
plot(SNRdB,C); hold on;
xlabel('SNRdB');
ylabel('Capacity');
