clear;
rng('shuffle')

numChan = 10:10:3000;
for l = 1: length(numChan)
%%% Defining Simulation parameters %%%
nUser = 2;          %%%% Number of users %%%%
nBS = numChan(l);          %%%% Number of BS antennas %%%%
nIter = 1000;        %%%% Number of channel iterations %%%%

%%%% Variable initializing %%%%
Chan_Users = zeros(nBS,nUser);                         %% Variable to store channel coefficient
Chan_hardening_Ratio = zeros(nUser,1);                 %% Variable to check the Chan_Hardening condition
Fav_Prop_Ratio = zeros(nUser,1);                       %% Variable to check Favourable Propagation condition

%%% Calculate channel norm and cross correlation %%%
y_iter=0;
g_iter = 0;
for i=1:nIter
    Chan_Users_1 = 1/sqrt(2)*(randn(nBS,nUser) + 1j*randn(nBS,nUser));    %% Create Rayleigh Channel matrix of size BS antennas X nUser of mean 0 and variance 1
    % 100x2 channels created, 100x1 for each user
    y_iter = y_iter + Chan_Users_1(:,1)'*Chan_Users_1(:,1);
    g_iter = g_iter + Chan_Users_1(:,2)'*Chan_Users_1(:,2);
        
end
y_mean = y_iter/nIter;
g_mean = g_iter/nIter;

Chan_Users = 1/sqrt(2)*(randn(nBS,nUser) + 1j*randn(nBS,nUser));    %% Create Rayleigh Channel matrix of size BS antennas X nUser of mean 0 and variance 1
    % 100x2 channels created, 100x1 for each user
    y_h = Chan_Users(:,1)'*Chan_Users(:,1);
    rath(l) = y_h/y_mean ;  % Channel hardening ratio for user 1
    ratg(l) = (Chan_Users(:,2)'*Chan_Users(:,2))/g_mean; % Channel hardening ratio for user 2

    

    %favourable propagation check
    fpc1(l)  = real((Chan_Users(:,1)'*Chan_Users(:,2))/sqrt(y_mean*g_mean));
    fpc2(l)  = real((Chan_Users(:,2)'*Chan_Users(:,1))/sqrt(y_mean*g_mean));
end
figure;
plot(numChan,rath,numChan,ratg); hold on;

plot(numChan,fpc1);hold on;
plot(numChan,fpc2);
xlabel('Number of antennas');
ylabel('ratio');
legend('channel herdening ratio of user 1','channel herdening ratio of user 2','Favourable propagation ratio 1','Favourable propagation ratio 2');


% disp('Channel Hardening Condition check')
% disp(Chan_hardening_Ratio)
% disp('Favourable Propagation Check')
% disp(Fav_Prop_Ratio)
