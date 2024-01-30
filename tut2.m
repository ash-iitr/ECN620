clc;
clear;
close all;

%% parameters
N = 10^4; % Number of symbols/frames
N_sc = 64; % Number of Sub-carriers
N_fft = 64; % FFT size (Assuming standard specifications)
CP = 6; % Cyclic Prefix length
Ch_length = 3; % Channel length [H(0) H(1) H(2)]
N_BitsPerSym = 64; % Number of bits per OFDM symbol (Same as the number of Sub-carriers for BPSK)
snr_dB = -5:5:55; % snr values in dB
snr_linear = 10.^(snr_dB/10);
ber = zeros(length(snr_dB),1);
% loop through snrs
for k=1:length(snr_dB)
%% OFDM Transmitter
% generate input bits and modulate using BPSK
bits = randi([0 1],[N*N_sc,1]);
bpsk_mod = 2*bits -1;
% Taking IFFT, normalize the power
bpsk_sp = reshape(bpsk_mod,N_fft,[]);
bpsk_ifft = ifft(bpsk_sp);
% Add CP
bpsk_ifft_cp = zeros(N_fft+CP,N);
bpsk_ifft_cp(1:N_fft,:) = bpsk_ifft(1:N_fft,:);
bpsk_ifft_cp(N_fft+1:end,:) = bpsk_ifft(1:CP,:);
% generate a multipath channel 
h = 1/sqrt(2)*(randn(Ch_length,1) + 1i*randn(Ch_length,1));
% Frequency response of the channel (To use at receiver)

% Convolution of Each symbol with the Random channel

% Concatenation of the symbols
bpsk_ifft_cp_ps = reshape(bpsk_ifft_cp,[],1);
% Gaussian noise with mean=0 and var=1
noise = 1/sqrt(2)*1/sqrt(snr_linear(k)).*(randn(length(bpsk_ifft_cp_ps),1));
% Adding noise with the input
rx_noisy_op = bpsk_ifft_cp_ps + noise;



%% OFDM Receiver
% Formatting the received vector into symbols
rx_bpsk_ifft_cp_sp = reshape(rx_noisy_op,N_fft+CP,[]);
% Remove CP
rx_bpsk_ifft_sp = rx_bpsk_ifft_cp_sp(1:N_fft,:);
% Converting into frequency domain
rx_bpsk_fft_sp = fft(rx_bpsk_ifft_sp);
% Single Tap Equalizer

% Extracting the required data
rx_bpsk_fft_ps = reshape(rx_bpsk_fft_sp,[],1);
% Demodulation
demod_bits = rx_bpsk_fft_ps>0;
demod_bits_double = zeros(length(demod_bits),1);
% Convertion into bits
for l=1:length(demod_bits)
    if(demod_bits(l))
        demod_bits_double(l)=1;
    end

end
% Count the number of errors
ber(k) = sum(demod_bits_double~=bits)/(N*N_sc);
theory_ber(k) = 0.5*(1-sqrt(((CP/N_sc)*snr_linear(k))/((CP/N_sc)*snr_linear(k)+2)));

end
%% Plotting
% Simulated BER

% Theoretical BER
%theory_ber = 0.5*(1+sqrt(snr_linear./snr_linear+2));
% Plot theoretical BER and simulated BER vs SNR
plot(snr_dB,ber,'-o');hold on;
plot(snr_dB,theory_ber,'-*');
xlabel('snrdB');
ylabel('BER');
legend('Simulation', 'Analytical');
title('BPSK OFDM BER vs SNR');
grid on;
