%% script_test_fcn_AVAR_avarI.m
% This script tests the function 'fcn_AVAR_avarI' on random walk+white noise
%
% This script was written on 2021_05_14 by Satya Prasad
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

number_of_time_steps = 16385;
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
    fcn_AVAR_generateIrregularWhiteNoise(power_spectral_density,...
    sampling_frequency,mean_white_noise,number_of_time_steps,...
    upsampling_factor,ireg_time_vector); % generate white noise
ireg_random_walk = ...
    fcn_AVAR_generateIrregularRandomWalk(random_walk_coefficient,...
    sampling_frequency,number_of_time_steps,...
    upsampling_factor,ireg_time_vector); % generate random walk
ireg_noise_signal = ireg_random_walk+ireg_white_noise; % noise signal

list_of_correlation_time = list_of_correlation_intervals/sampling_frequency;
sampling_interval = 1/sampling_frequency;
min_time = 0;
max_time = (number_of_time_steps-1)*sampling_interval;

% clear unused variables
clear number_of_time_steps p list_of_correlation_intervals
clear random_walk_coefficient power_spectral_density mean_white_noise sampling_frequency upsampling_factor
clear ireg_white_noise ireg_random_walk

lol = fcn_AVAR_avarI(ireg_noise_signal,ireg_time_vector,list_of_correlation_time,...
    sampling_interval,min_time,max_time,12345); % calculate allan variance of irregularly sampled data

whos
