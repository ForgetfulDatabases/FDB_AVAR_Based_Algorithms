%% script_test_fcn_AVAR_favarI.m
% This script tests the function 'fcn_AVAR_favarI' on random walk+white noise
%
% This script was written on 2021_05_15 by Satya Prasad
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

number_of_time_steps = 16385; % Length of the data
p = floor(log2(number_of_time_steps));
list_of_correlation_intervals = 2.^(1:(p-1))'; % list of correlation intervals

%% Example 1: Random Walk+White Noise
random_walk_coefficient = 0.025; % [unit/sqrt(s)]
power_spectral_density  = 0.0025; % PSD of white noise [unit^2 s]
mean_white_noise        = 0; % mean of white noise
sampling_frequency      = 20; % [Hz]
upsampling_factor       = 10;
ireg_time_vector        = NaN;

[ireg_white_noise, ireg_time_vector] = ...
    fcn_AVAR_generateIrregularWhiteNoise(power_spectral_density, ...
    sampling_frequency, mean_white_noise, number_of_time_steps, ...
    upsampling_factor, ireg_time_vector); % generate white noise
ireg_random_walk = ...
    fcn_AVAR_generateIrregularRandomWalk(random_walk_coefficient, ...
    sampling_frequency, number_of_time_steps, ...
    upsampling_factor, ireg_time_vector); % generate random walk
ireg_noise_signal = ireg_random_walk+ireg_white_noise; % noise signal

% Preprocessing step to arrange irregularly sampled signal into regularly
% sampled weighted data
time_vector = (1/sampling_frequency)*(0:number_of_time_steps-1)';
length_of_input_data = length(ireg_noise_signal); % length of the input data
data    = nan(length_of_input_data,1); % initialize input data
weights = nan(length_of_input_data,1); % initialize input weights
for i = 1:length_of_input_data
    if i < length_of_input_data
        temp_data = ireg_noise_signal(ireg_time_vector>=time_vector(i) & ...
                                     ireg_time_vector<time_vector(i+1));
    else
        temp_data = ireg_noise_signal(ireg_time_vector>=time_vector(i));
    end
    weights(i) = numel(temp_data);
    if isempty(temp_data)
        data(i) = 0;
    else
        data(i) = mean(temp_data);
    end
end % END: for loop over the all the data

clear number_of_time_steps p
clear random_walk_coefficient power_spectral_density mean_white_noise sampling_frequency upsampling_factor
clear ireg_white_noise ireg_random_walk ireg_noise_signal time_vector ireg_time_vector

fcn_AVAR_favarI(data, weights, list_of_correlation_intervals, 12345); % calculate allan variance of irregularly sampled data using FAVAR

whos
