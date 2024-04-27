rng(111);% Set RNG state for repeatability
Nclmax= 5; % The maximum number of scattering cluster 
Nclmin= 1; % The minimum number of scattering cluster
Nraymax= 4; % The maximum number of rays in each cluster
Nraymin= 1; % The minimum number of rays in each cluster
Nt=  16; % Number of Tx antenna 
Nr= 64 ;  % Number of Rx antenna
sigma= 10*pi/180; % Angular spread
% Compute the number of scattering cluster as per the given distribution
Ncl= randi([Nclmin Nclmax],1,1);
%% Compute the Number of rays in each cluster  as per the given distribution
Nray= randi([Nraymin Nraymax],1,1);
%% Compute the normalization factor
gamma= sqrt((Nr*Nt)/(Ncl*Nray));
H =0;
%% Loop over number of clusters
for c = 1: Ncl
    %% Loop over number of rays in each clsuter
for k = 1:Nray
    %% Generate path gain associated with the k-th ray in i-th clusters as per the given distribution
Beta= abs(randn(1) + 1j * randn(1)) / sqrt(2);
% Compute the average elevation angle as per the given distribution
mue=  2 * pi * rand(1);
% Compute the average azimuth angle as per the given distribution
mua= -pi/2 + pi * rand(1);
%% Compute the elevation angles for the user and BS as per the given distribution
thetaUE = lprnd(1, 1, mue, sigma);
thetaBs =  lprnd(1, 1, mue, sigma);
%% Compute the azimuth angles for the user and BS as per the given distribution
PhiUE   = lprnd(1, 1, mua, sigma);
PhiBs = lprnd(1, 1, mua, sigma);
%% Generate array response vector for the user 

aUE = sqrt(1/Nt) * exp(1j * pi * (0:Nt-1) * (c* sin(thetaUE) * sin(PhiUE) + k*cos(thetaUE))).';

%% Generate array response vector for the base station
 
aBS =sqrt(1/Nr) * exp(1j * pi * (0:Nr-1) * (c* sin(thetaBs) * sin(PhiBs) + k*cos(thetaBs))).';


H= H +   gamma*Beta * aBS * aUE';
%% Generate the channel and accumulate
end
end
%% Display the channel matrtrix
%disp("H: ", H);

%% Laplacian distribution function

function [x,y]  = lprnd(m, n, mu, sigma)
%Check inputs
if nargin < 2
    error('At least two inputs are required');
end
if nargin == 2
    mu = 0; sigma = 1;
end
if nargin == 3
    sigma = 1;
end
% Generate Laplacian noise
u = rand(m, n)-0.5;
b = sigma / sqrt(2);
x = mu - b * sign(u).* log(1- 2* abs(u));
y = mu - b * sign(u).* log(1- 2* abs(u));
end
