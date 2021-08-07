%% script_compare_fcn_AVAR_dFavar_AND_fcn_AVAR_favar.m
% This script compares computation time of 'fcn_AVAR_dFavar' with 'fcn_AVAR_favar'
%
% This script was written on 2021_07_15 by Satya Prasad
% Questions or comments? szm888@psu.edu
%

%% Prepare workspace
clear all %#ok<CLALL>
close all
clc

%% Add path
addpath('.\functions')

%% Intialization
rng('default') % set random seeds

%% Define constants and parameters
random_walk_coefficient = 0.025; % [unit/sqrt(s)]
power_spectral_density  = 0.05; % psd of white noise
mean_white_noise        = 0; % mean of white noise
sampling_frequency      = 20; % [Hz]

elapsed_time_matrix = NaN(2,17);
for j = 2:18
horizon_size         = (2^j)+1;
number_of_time_steps = horizon_size+2^10;

%% Noise generation: Random Walk added to White Noise
time_vector = (1/sampling_frequency)*(0:(number_of_time_steps-1));

random_walk  = fcn_AVAR_generateRandomWalk(random_walk_coefficient,....
               sampling_frequency,number_of_time_steps); % generate random walk
white_noise  = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
               sampling_frequency,mean_white_noise,number_of_time_steps); % generate white noise
noise_signal = random_walk + white_noise; % noise signal

%% Calculation of possible correlation intervals
p = floor(log2(horizon_size));
list_of_correlation_intervals = 2.^(1:(p-1))'; % A vector of correlation intervals

%% Dynamic Allan variance: D-FAVAR
tic
for i = horizon_size:number_of_time_steps
    if i == horizon_size
        clear fcn_AVAR_dFavar
        data = noise_signal(1:horizon_size);
    else
        data = noise_signal(i-horizon_size:i);
    end
    dFavar = fcn_AVAR_dFavar(data,list_of_correlation_intervals); % calculate allan variance
end
dFavar_time = toc;
elapsed_time_matrix(1,j-1) = dFavar_time/(number_of_time_steps-horizon_size+1);

%% Dynamic Allan variance: FAVAR
tic
for i = horizon_size:number_of_time_steps
    data = noise_signal(i-horizon_size+1:i);
    favar = fcn_AVAR_favar(data,list_of_correlation_intervals); % calculate allan variance
end
favar_time = toc;
elapsed_time_matrix(2,j-1) = favar_time/(number_of_time_steps-horizon_size+1);

end

%% Plots
%% Ratio of wall time
figure(12345)
plot(2.^(2:18)+1,elapsed_time_matrix(2,:)./elapsed_time_matrix(1,:),'k--','LineWidth',1.2)
hold on
yline(1,'k','LineWidth',1.2);
grid on
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('Horizon size, $N$','Interpreter','Latex','FontSize',13)
ylabel('Ratio of wall time','Interpreter','Latex','FontSize',13)
title('D-FAVAR vs FAVAR','Interpreter','Latex','FontSize',13)
% Get handle to current axes.
ax = gca;
% Set x and y font sizes.
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;
