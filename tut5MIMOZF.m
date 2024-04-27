%%%%------ Simulating SER performance of MIMO ZF Detector and compare it with its real equivalent counter part ------%%%%

rng(100);
N=10^4;                                                                    % No. of input Symbols
Nt=4;                                                                      % No. of Transmit antenna
Nr=4;                                                                      % No. of Receive antenna
M=4;                                                                       % Modulation Order 

%% OVERALL SYSTEM
EsNo=0:3:30; %SNR in dB
SER_ZF=zeros(1,length(EsNo));
SER_ZF_real=zeros(1,length(EsNo));
for iter_1=1:length(EsNo)
    sym_err_count=0;
    sym_err_count_real=0;
    for iter_2=1:ceil(N/Nt)
        %% DATA
        data=randi([0 M-1],Nt,1);
        Tx_data = qammod(data, M, 'Gray', 'UnitAveragePower', true);
        %% Transmiter
        Noise = 1/sqrt(2)*(randn(Nr,1) + sqrt(-1)*randn(Nr,1));
        H = 1/sqrt(2)*(randn(Nr,Nt) + sqrt(-1)*randn(Nr,Nt));
        %% Received Signal
        % Model the MIMO system model ensuring the SNR is same as the runnuing loop of EsNo
        Rx_data = (H*Tx_data)+sqrt(Nt)*10^(-EsNo(iter_1)/(20))*Noise;%
        %% Complex ZF Detector
        % Find the ZF solution
        ZF_sol = ((H'*H)\H'*Rx_data);  % fill the code
        % Demodulate the ZF solition
        demod_out_ZF=  qamdemod(ZF_sol, M, 'Gray', 'UnitAveragePower', true);  % fill the code
        %% Count the total number of symbols in error
        sym_err_count=  sym_err_count + symerr(data,demod_out_ZF);  % fill the code (Hint:- use symerr function)
        %% Real Equivalent ZF Detector 
        % Convert channel matrix to its real equivalent form
        H_real=  [real(H), -imag(H); imag(H), real(H)]; % fill the code
        % Convert received signal to its real equivalent form
        Rx_data_real=  [real(Rx_data); imag(Rx_data)];  % fill the code
        % Find the ZF solution
        ZF_sol_real = ((H_real'*H_real)\H_real'*Rx_data_real);   % fill the code
        % Convert it to complex form
        ZF_sol_eqv=  ZF_sol_real(1:Nt) +1j*ZF_sol_real(Nt+1:2*Nt)  ;  % fill the code
        % Demodulate the ZF solution
        demod_out_ZF_real= qamdemod(ZF_sol_eqv, M, 'Gray', 'UnitAveragePower', true) ;  % fill the code
        %% Count the number of symbols in error
       sym_err_count_real= sym_err_count_real + symerr(data,demod_out_ZF_real);  % fill the code (Hint:- use symerr function)     
    end
    SER_ZF(iter_1)=sym_err_count/(ceil(N/Nt)*Nt);
    SER_ZF_real(iter_1)=sym_err_count_real/(ceil(N/Nt)*Nt);
end

figure(1)
semilogy(EsNo,SER_ZF,'r -o', EsNo,SER_ZF_real,'b');  
legend('SER ZF','SER ZF Real');
xlabel('Es/No(dB)');
ylabel('Symbol Error Rate');
title('SER performance of MIMO ZF detector');
grid on
