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
theory_ber = zeros(length(snr_dB),1);
% loop through snrs
for k=1:length(snr_dB)
    %% OFDM Transmitter
    % generate input bits and modulate using BPSK
    bits = randi([0 1],[N*N_sc,1]);
    bpsk_mod = 2*bits -1;
    % Taking IFFT, normalize the power
    bpsk_sp = reshape(bpsk_mod,N_fft,[]);

    for i=1:N
    
    bpsk_ifft = (ifft(bpsk_sp(:,i)));
    % Add CP
    bpsk_ifft_cp = zeros(N_sc+CP,1);
    bpsk_ifft_cp(CP+1:end,1) = bpsk_ifft(1:N_sc,1);
    bpsk_ifft_cp(1:CP,1) = bpsk_ifft(N_sc-CP+1:end,1);
    % % Parallel to serial conversion
    bpsk_ifft_cp_ps = reshape(bpsk_ifft_cp,1,[]);
    % generate a multipath channel 
    
    h(1,1:Ch_length) = (Ch_length/sqrt(2))*(randn(1,Ch_length) + 1i*randn(1,Ch_length));
    
    % % Convolution of Each symbol with the Random channel
    
     y = cconv(h,bpsk_ifft_cp_ps,N_sc+CP);
    
    % % Gaussian noise with mean=0 and var=1 
    
      noise = (1/sqrt(snr_linear(k))).* ( randn(1,N_sc+CP) + 1j* randn(1,N_sc+CP) );
    
    % % Adding noise with the input
    
     y_rx = y + noise;


    %% OFDM Receiver
    % Formatting the received vector into symbols                    
    % % Remove CP
    
    rx_bpsk_ifft_sp = y_rx(CP+1:end);
    
    % % Converting into frequency domain
    rx_bpsk_fft_sp = fft(rx_bpsk_ifft_sp);
    % % Single Tap Equalizer
    Channel_fft = fft(h,N_sc);
    eq_fft_bpsk_rx = rx_bpsk_fft_sp./Channel_fft;
    % % Extracting the required data
    extracted_data(:,i) = eq_fft_bpsk_rx';
    
    end

rx_bpsk_fft_ps = reshape(extracted_data,[],1);
% % Demodulation
demod_bits = rx_bpsk_fft_ps>0;
% % % Convertion into bits 
% % % Count the number of errors
ber(k) = sum(demod_bits~=bits)/(N*N_sc);
%theory_ber(k) = 0.5*(1-sqrt(((Ch_length/N_sc)*snr_linear(k))/((Ch_length/N_sc)*snr_linear(k)+1)));
theory_ber(k) = 0.5*(1-sqrt((snr_linear(k))/((snr_linear(k))+2)));
end
%% Plotting
% Simulated BER

% Theoretical BER
%theory_ber = 0.5*(1+sqrt(snr_linear./snr_linear+2));
% Plot theoretical BER and simulated BER vs SNR

semilogy(snr_dB,ber,'-o');hold on;
semilogy(snr_dB,theory_ber,'-*');
xlabel('snrdB');
ylabel('BER');
legend('Simulation', 'Analytical');
title('BPSK OFDM BER vs SNR');
grid on;
