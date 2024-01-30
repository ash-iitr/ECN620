%% Plot the ergodic capacity of the multiple-input multiple-output (MIMO) channel. 

NT = 5; % Number of transmit antennas
NR = 5; % Number of receive antennas
N = 1000; % Number of iterations
SNR_dB = -10:2:40; % SNR in dB scale
SNR_linear = 10.^(SNR_dB/10);
% Create a variable to store Capacity
C = zeros(length(SNR_linear),N);
% Calculate capacity for each iteration and each SNR

% h = zeros(NR,NT,N);
for j =1:N
h(:,:,j) = 1/sqrt(2) .* (randn(NR,NT) + 1i*randn(NR,NT));
end

for i=1:length(SNR_dB)
    for j=1:N
    C(i,j) = log2(det(eye(NR) + (SNR_linear(i)/NT).*(h(:,:,j)*h(:,:,j)')) );
    end
end


CE = sum(real(C),2)./N;
%Plot capacity vs SNR
plot(SNR_dB,CE); hold on;
xlabel('SNRdB');
ylabel('Capacity');
