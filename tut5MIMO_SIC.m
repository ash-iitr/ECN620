%%% Simulate the symbol error rate performance of MIMO system with SIC detector %%%
rng(100);
N=2*10^3;                                                                    % No. of input Bits
Nt=4;                                                                      % No. of Transmit antenna
Nr=4;                                                                      % No. of Receive antenna
M=4;                                                                       % Modulation Order 
%% OVERALL SYSTEM
EsNo=0:3:30; % SNR in dB
for iter_1=1:length(EsNo)
    sym_err_count=0;
    for iter_2=1:ceil(N/Nt)
        %% Transmit Data
        data=randi([0 M-1],Nt,1);
        Tx_data = qammod(data, M, 'Gray', 'UnitAveragePower', true);
        %% Received Signal
        % Model the MIMO Channel Matrix
        H = 1/sqrt(2)*(randn(Nr,Nt) + sqrt(-1)*randn(Nr,Nt));%Fill the entries of channel matrix using the Rayleigh Fading
        % Model the AWGN noise vecto
        Noise = 1/sqrt(2)*(randn(Nr,1) + sqrt(-1)*randn(Nr,1));%Fill the enties of noise vector assuming complex Gaussian with mean 0 and variance 1
        % Express the MIMO received signal vector ensuring the SNR is maintainted as per EsNo loop
        Rx_data = (H*Tx_data)+sqrt(Nt)*10^(-EsNo(iter_1)/(20))*Noise;  
        %% Receiver
        % Equivalent MIMO System Model
        % Perform QR decomposition of the channle matrix
        [Q,R] = qr(H);
        % fill the code (Hint:- use qr function)
        % Find the equivalent received signal vector
        Z= Q'*Rx_data;% fill the code
        X = zeros(Nt,1);
        %% SIC detector
         for iter=Nt:-1:1
           % Find the symbol estimate starting from bottom most antenna
          temp = qamdemod(Z(iter)/R(iter,iter),M,'Gray','UnitAveragePower', true);
           temp2 = qammod(temp, M, 'Gray', 'UnitAveragePower', true);
           % Map it to the nearest constellation point and save it as the SIC solution for the selected antenna index
           X(iter,1) =  temp2;
           % Use the decoded symbol to remove its effect from the received signal vector
            Z = Z - R(:,iter)*X(iter,1);
           % iterate the process          
         end   
        % Demodulate the SIC solution   
        demod_out_sic=qamdemod(X, M, 'Gray', 'UnitAveragePower', true);   % fill the code
        % Count the total number of symbols in error
        sym_err_count= sym_err_count + symerr(data,demod_out_sic); % fill the code (Hint:- use symerr function)
    end
    SER_SIC(iter_1)=sym_err_count/(ceil(N/Nt)*Nt);
end
figure(1)
semilogy(EsNo,SER_SIC,'b -s'); 
legend('SER SIC');
xlabel('Es/No(dB)');
ylabel('Symbol Error Rate');
title('SER performance of MIMO SIC detector');
grid on

