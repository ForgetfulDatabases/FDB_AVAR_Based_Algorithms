%% script_plot_fcn_AVAR_dFavar_AND_fcn_AVAR_favar.m
%  This script plots DAVAR estimated using 'fcn_AVAR_dFavar' and
%  'fcn_AVAR_favar'
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

horizon_size         = (2^17)+1;
number_of_time_steps = (2^(17+1))+1;
dFavar_matrix = NaN(16,number_of_time_steps-horizon_size+1);
favar_matrix  = NaN(16,number_of_time_steps-horizon_size+1);

%% Noise generation: Random Walk added to White Noise
time_vector = (1/sampling_frequency)*(0:(number_of_time_steps-1));

random_walk1  = fcn_AVAR_generateRandomWalk(2*random_walk_coefficient,....
               sampling_frequency,horizon_size); % generate random walk
white_noise1  = fcn_AVAR_generateWhiteNoise(10*power_spectral_density,...
               sampling_frequency,mean_white_noise,horizon_size); % generate white noise
random_walk2  = fcn_AVAR_generateRandomWalk(random_walk_coefficient,....
               sampling_frequency,number_of_time_steps-horizon_size); % generate random walk
white_noise2  = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
               sampling_frequency,mean_white_noise,number_of_time_steps-horizon_size); % generate white noise
noise_signal = [random_walk1; random_walk2] + [white_noise1; white_noise2]; % noise signal

%% Calculation of possible correlation intervals
p = floor(log2(horizon_size));
list_of_correlation_intervals = 2.^(1:(p-1))'; % A vector of correlation intervals

%% Dynamic Allan variance: D-FAVAR
for i = horizon_size:number_of_time_steps
    if i == horizon_size
        clear fcn_AVAR_dFavar
        data = noise_signal(1:horizon_size);
    else
        data = noise_signal(i-horizon_size:i);
    end
    dFavar = fcn_AVAR_dFavar(data,list_of_correlation_intervals); % calculate allan variance
    dFavar_matrix(:,i-horizon_size+1) = dFavar;
end

%% Dynamic Allan variance: FAVAR
for i = horizon_size:number_of_time_steps
    data  = noise_signal(i-horizon_size+1:i);
    favar = fcn_AVAR_favar(data,list_of_correlation_intervals); % calculate allan variance
    favar_matrix(:,i-horizon_size+1) = favar;
end

%% Plots
figure(12345)
clf
plot(list_of_correlation_intervals, favar_matrix(:,1), 'ko', 'Markersize', 13)
hold on
plot(list_of_correlation_intervals, dFavar_matrix(:,1), 'k*', 'Markersize', 9)
plot(list_of_correlation_intervals, favar_matrix(:,end), 'kd', 'Markersize', 13)
plot(list_of_correlation_intervals, dFavar_matrix(:,end), 'k.', 'Markersize', 13)
plot(list_of_correlation_intervals, fcn_AVAR_favar(noise_signal,list_of_correlation_intervals),...
     'k--', 'LineWidth', 1.2)
grid on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Correlation interval, $m$', 'Interpreter', 'Latex', 'FontSize', 13)
ylabel('AVAR', 'Interpreter', 'Latex', 'FontSize', 13)
legend('FAVAR', 'D-FAVAR', 'FAVAR', 'D-FAVAR', 'FAVAR', ...
       'Interpreter', 'Latex', 'FontSize', 13, 'Location', 'best', ...
       'NumColumns', 3)
% Get handle to current axes.
ax = gca;
% Set x and y font sizes.
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;
