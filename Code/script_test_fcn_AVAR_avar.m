%% script_test_fcn_AVAR_avar.m
% This script tests 'fcn_AVAR_avar' on different random walk+white noise
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

random_walk = fcn_AVAR_generateRandomWalk(random_walk_coefficient,....
              sampling_frequency,number_of_time_steps); % generate random walk
white_noise = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
              sampling_frequency,mean_white_noise,number_of_time_steps); % generate white noise
noise_signal = random_walk+white_noise; % noise signal

% clear unused variables
clear number_of_time_steps p
clear random_walk_coefficient power_spectral_density mean_white_noise sampling_frequency
clear random_walk white_noise

fcn_AVAR_avar(noise_signal,list_of_correlation_intervals,12345); % calculate allan variance

whos
