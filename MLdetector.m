% Test cases didn't pass
% MIMO Maximum likelihood (ML) detector
clc;
clear;
close all;

rng(9);

Nt = 2; % Tx antennas 
Nr = 2;% Rx antennas
num_bits = 2* 10^5; % length of bits 
num_realizations = 3; % Number of realizations

% Define SNR range for plotting
SNR_dB = 0:2:12;
SNR_linear = 10.^(SNR_dB/10);
BER_avg = zeros(length(SNR_dB), 1);

for realization = 1:num_realizations
    % Generate bit stream of dimension: num_bitsx1
    data = randi([0 1],num_bits,1);

    % Map data to BPSK symbols (1 --> +1, 0 --> -1)
    symbols = 2*data -1;

    % Create 2x2 MIMO rayleigh fading channel matrix
    H = sqrt(1/2)*(randn(Nr,Nt)+ 1j*randn(Nr,Nt));

    % Reshape data for MIMO transmission from the Tx antennas
    data_mtx = reshape(symbols, [], Nt);

    
    BER = zeros(length(SNR_dB), 1);

    % Loop through various SNR values to generate BER performance
    for ii = 1:length(SNR_dB)
        noise = ( 1/sqrt(SNR_linear(ii)/2))*(randn(size(data_mtx)) + 1j*randn(size(data_mtx))); % generate AWGN noise of dimension same as data_mtx
        
        received_symbols = data_mtx*H + noise; % write the received symbol equation
        
        detected_data_mtx = zeros(size(received_symbols));

        for jj = 1:size(received_symbols, 1)
            
            detected_data_mtx(jj,:) = ml_detector(received_symbols(jj,:),H);    % call the ML detector function to detect from received_symbols block
            
        end

         detected_bits = (detected_data_mtx+1)./2; % perform BPSK demodulation on detected_data_mtx
         detected_bits_reshaped = reshape(detected_bits,[],1);
         
         errors = sum(detected_bits_reshaped~=data); % compute the total bits found in error
         BER(ii) = errors/num_bits;% compute BER
     end
    % 
     BER_avg = BER_avg + BER;
end

% Average BER over realizations
BER_avg = BER_avg / num_realizations;

% Plot BER vs SNR curve
figure(1);
semilogy(SNR_dB, BER_avg, '-^', 'LineWidth', 3);
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('Average BER vs SNR for BPSK in 2x2 MIMO (ML Detection)');
grid on;

% Function for Maximum Likelihood Detection
function [detected_symbol] = ml_detector(y, H)

    % Creat a lookup table for 2x2 MIMO config considering BPSK modulated symbols
    possible_symbols = [1,1;1,-1;-1,1;-1,-1];
    distances = zeros(size(possible_symbols,1), 1); % store the norm distances

    for i=1:length(distances) % Loop through all combinations and calculate distance
        
    distances(i) = norm(y - H*possible_symbols(i,:)')^2; % compute distances. for eg:- || y - Hx ||^2, here x can be taken from the look up table
    
    end

    % Find index of minimum distance combination
  [~, min_index] = min(distances);
  detected_symbol = possible_symbols(min_index,:)'; % symbol from the look-up table corresponding to the min_index
  
end
