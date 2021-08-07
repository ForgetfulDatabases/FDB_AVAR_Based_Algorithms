%% script_compare_fcn_AVAR_dFavarI_AND_fcn_AVAR_favarI.m
% This script compares computation time of 'fcn_AVAR_dFavarI' with 'fcn_AVAR_favarI' 
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
upsampling_factor       = 10;

elapsed_time_matrix = NaN(2,17);
for j = 2:18
ireg_time_vector     = NaN;
horizon_size         = (2^j)+1;
number_of_time_steps = horizon_size+2^10;

%% Noise generation: Random Walk added to white noise
time_vector = (1/sampling_frequency)*(0:(number_of_time_steps-1))';

[ireg_random_walk, ireg_time_vector] = ...
    fcn_AVAR_generateIrregularRandomWalk(random_walk_coefficient, ...
    sampling_frequency, number_of_time_steps, ...
    upsampling_factor, ireg_time_vector); % generate random walk
ireg_white_noise = fcn_AVAR_generateIrregularWhiteNoise(power_spectral_density, ...
    sampling_frequency, mean_white_noise, number_of_time_steps, ...
    upsampling_factor, ireg_time_vector); % generate white noise
noise_signal = ireg_random_walk+ireg_white_noise; % random walk plus white noise (irregularly sampled)

%% Calculation of possible correlation intervals
p = floor(log2(horizon_size));
list_of_correlation_intervals = 2.^(1:(p-1))'; % A vector of correlation interval(s)

%% Data Preprocessing into Intervals
input_data    = zeros(number_of_time_steps,1); % initialize input data with zeros
input_weights = zeros(number_of_time_steps,1); % initialize input weights with zeros

for i = 1:number_of_time_steps
    if i < number_of_time_steps
        data = noise_signal(ireg_time_vector>=time_vector(i) & ireg_time_vector<time_vector(i+1));
    else
        data = noise_signal(ireg_time_vector>=time_vector(i));
    end
    input_weights(i) = numel(data);
    if ~isempty(data)
        input_data(i) = mean(data);
    end
end % END: for loop over the all the data

%% Dynamic Allan variance: D-FAVAR-I
tic
for i = horizon_size:number_of_time_steps
    if i == horizon_size
        clear fcn_AVAR_dFavarI
        data    = input_data(1:horizon_size);
        weights = input_weights(1:horizon_size);
    else
        data    = input_data(i-horizon_size:i);
        weights = input_weights(i-horizon_size:i);
    end
    dFavar      = fcn_AVAR_dFavarI(data,weights,list_of_correlation_intervals); % calculate allan variance
end
dFavar_time = toc;
elapsed_time_matrix(1,j-1) = dFavar_time/(number_of_time_steps-horizon_size+1);

%% Dynamic Allan variance: FAVAR-I
tic
for i = horizon_size:number_of_time_steps
    data    = noise_signal(i-horizon_size+1:i);
    weights = input_weights(i-horizon_size+1:i);
    favar   = fcn_AVAR_favarI(data,weights,list_of_correlation_intervals); % calculate allan variance
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
xlabel('Horizon length, $N$','Interpreter','Latex','FontSize',13)
ylabel('Ratio of wall time','Interpreter','Latex','FontSize',13)
title('D-FAVAR-I vs FAVAR-I','Interpreter','Latex','FontSize',13)
% Get handle to current axes.
ax = gca;
% Set x and y font sizes.
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;
