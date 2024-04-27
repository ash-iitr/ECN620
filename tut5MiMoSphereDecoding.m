rng(18);

Nt = 2; % Tx antennas 
Nr = 2; % Rx antennas
num_bits = 2e5; % length of bits 
num_realizations = 1; % Number of realizations

% Define SNR range for plotting
SNR_dB = 0:3:12;

BER_ML = zeros(length(SNR_dB), 1);
BER_SD = zeros(length(SNR_dB), 1);

% Loop through various SNR values to generate BER performance
for ii = 1:length(SNR_dB) 
    ii
    % Generate bit stream
    data = randi([0 1], num_bits, 1);

    % Map data to BPSK symbols (+1, -1)
    symbols = 2 * data - 1;

    
    % Reshape data for MIMO transmission from the Tx antennas
    data_mtx = reshape(symbols, Nt, []);
    errors_ML=0;
    errors_SD=0;
   for realization = 1:num_realizations
    % Create 2x2 MIMO channel matrix
    H = sqrt(1/2) * (randn(Nr,Nt) + 1j*randn(Nr,Nt));
    noise = sqrt(1/2) * (randn(Nr, length(data)/Nt) + 1j*randn(Nr, length(data)/Nt));
    received_symbols =  H*data_mtx + sqrt(Nt)*10^(-SNR_dB(ii)/20) * noise; % y = Hx + n
   
      % Perform QR decomposition on H
      [Q, R] = qr(H);
      
      % Apply QR decomposition to received symbols
      z = Q'*received_symbols;
        
        
       %% Optimal ML Detector 
        detected_ML_data_mtx = zeros(size(received_symbols));

        for jj = 1:size(received_symbols, 2)
            detected_ML_data_mtx(:,jj) = ml_detector(z(:,jj), R); % ML detector
        end

        detected_ML_bits = detected_ML_data_mtx(:) > 0; % BPSK demodulation
        errors_ML = errors_ML+ sum(abs(data - detected_ML_bits)); 
       % BER_ML(ii) = errors_ML / num_bits;
        
      
        %% SD Detector
       
      
      detected_SD_data_mtx = zeros(size(received_symbols));
      for jj = 1:size(received_symbols,2)
        % detected_SD_data_mtx(:,jj) = sphere_decoder(z(:,jj), R); % Sphere decoder
        radius =inf;
        level=Nt;
        cons_syms=[1 -1];
        detected_SD_data_mtx(:,jj) = SD(z(:,jj), R, cons_syms, radius); % An example recursive function with the arguments z(:,jj), R, cons_syms, radius, etc.
       end
      
      detected_SD_bits = detected_SD_data_mtx(:) > 0; % BPSK demodulation
      errors_SD = errors_SD+ sum(abs(data - detected_SD_bits)); 
     % BER_SD(ii) = errors_SD / num_bits;

        
        
   end
    BER_ML(ii) = errors_ML / (num_bits*num_realizations);
     BER_SD(ii) = errors_SD / (num_bits*num_realizations);
    % BER_avg_ML = BER_avg_ML + BER_ML;
    % BER_avg_SD = BER_avg_SD + BER_SD;
end

% Average BER over realizations
% BER_avg_ML = BER_avg_ML / num_realizations;
% BER_avg_SD = BER_avg_SD / num_realizations;
% Plot BER vs SNR curve
semilogy(SNR_dB, BER_ML, '-^', 'LineWidth', 3);
hold on
semilogy(SNR_dB, BER_SD, '-o', 'LineWidth', 3);
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('Average BER vs SNR for BPSK in 2x2 MIMO (ML and SD Detection)');
legend('ML', 'SD');
grid on;

% Function for Maximum Likelihood Detection
function [detected_symbol] = ml_detector(z, R)
    % BPSK lookup table for 2x2 MIMO config
    possible_symbols = [-1 1 -1 1;-1 -1 1 1];

    distances = zeros(size(possible_symbols,1), 1); % store the norm distances

    % Loop through all combinations and calculate distance
    for i = 1:size(possible_symbols, 2)
        received_symbol = possible_symbols(:,i);
        distances(i) = norm(z - R* received_symbol)^2;
    end

    % Find index of minimum distance combination
    [~, min_index] = min(distances);
    detected_symbol = possible_symbols(:,min_index);
end

  
    % write your sphere decoder function here
%function [out_SD, d_CSE]= SD(z(:,jj), R, cons_syms, radius)


 % your code here
function [out_SD, d_CSE]= SD(z, R, cons_syms, radius)
    
    Nt = size(R, 2);
    out_SD = zeros(Nt, 1);
    d_CSE = inf;
    
    function sp_decode(depth, branch, branch_Sym)
        if depth > Nt    
        distance = norm(z - R * branch_Sym)^2;
          if distance < d_CSE
           d_CSE = distance;
           out_SD = branch_Sym;
          end
       else
            
         for symbol = cons_syms
           nBranch_Sym = branch_Sym;
           nBranch_Sym(depth) = symbol;
                
           nBranch = branch + abs(z(depth) - R(depth, depth) * symbol)^2;
                
           if nBranch < d_CSE + radius
            sp_decode(depth + 1, nBranch, nBranch_Sym);
           end
        end
     end
   end

sp_decode(1, 0, zeros(Nt, 1));
end


