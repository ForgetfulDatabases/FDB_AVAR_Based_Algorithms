%% script_compare_fcn_AVAR_favarI_AND_fcn_AVAR_avarI.m
%  This script compares 'fcn_AVAR_favarI' with 'fcn_AVAR_avarI'
%
% This script was written on 2021_06_30 by Satya Prasad
% Questions or comments? szm888@psu.edu
%

%% Prepare workspace
clear all %#ok<CLALL>
close all
clc

%% Add path
addpath('.\functions')

%% Intialization
rng('default')  % set random seeds

%% Define constants and parameters
random_walk_coefficient = 0.025; % [unit/sqrt(s)]
power_spectral_density  = 0.05; % [unit^2 s]
mean_white_noise        = 0; % mean of white noise
sampling_frequency      = 20; % [Hz]
upsampling_factor       = 25; % No units
sampling_interval       = 1/sampling_frequency; % [s]
% increase iterations to decrease the error
number_of_iterations    = 1; % number of monte-carlo simulations

favar_matrix = NaN(17,17);
avar_matrix  = NaN(17,17);
elapsed_time_matrix = NaN(3,17);
for p = 2:18
number_of_time_steps = (2^p)+1;

%% Noise generation: Random Walk added to White Noise
time_vector = (1/sampling_frequency)*(0:(number_of_time_steps-1))';

ireg_noise_signal = nan(number_of_time_steps,number_of_iterations); % variable to store noise signals
ireg_time_vector  = nan(number_of_time_steps,number_of_iterations); % variable to store time vector
for k = 1:number_of_iterations % BEGIN: for loop for Noise Signal iterations
[ireg_random_walk, ireg_time_vector(:,k)] = ...
    fcn_AVAR_generateIrregularRandomWalk(random_walk_coefficient, ...
    sampling_frequency, number_of_time_steps, upsampling_factor, ...
    NaN); % generate random walk
ireg_white_noise = fcn_AVAR_generateIrregularWhiteNoise(power_spectral_density, ...
    sampling_frequency, mean_white_noise, number_of_time_steps, ...
    upsampling_factor, ireg_time_vector(:,k)); % generate white noise
ireg_noise_signal(:,k) = ireg_random_walk+ireg_white_noise; % random walk plus white noise (irregularly sampled)
end % END: for loop for Noise Signal iterations

%% Calculation of possible correlation intervals
list_of_correlation_intervals = 2.^(1:(p-1))'; % A vector of correlation intervals
list_of_correlation_time      = list_of_correlation_intervals/sampling_frequency; % A vector of correlation times

%% Data Preprocessing into Intervals and calculation of weights and average
input_data    = zeros(number_of_time_steps,number_of_iterations); % initialize input data
input_weights = zeros(number_of_time_steps,number_of_iterations); % initialize input weights
tic
for k = 1:number_of_iterations
for i = 1:number_of_time_steps
    if i < number_of_time_steps
        data = ireg_noise_signal(ireg_time_vector(:,k)>=time_vector(i) & ...
                                 ireg_time_vector(:,k)<time_vector(i+1), k);
    else
        data = ireg_noise_signal(ireg_time_vector(:,k)>=time_vector(i) & ...
                                 ireg_time_vector(:,k)<time_vector(i)+sampling_interval, k);
    end
    input_weights(i,k) = numel(data);
    if ~isempty(data)
        input_data(i,k) = mean(data);
    end
end % END: for loop over the all the data
end
elapsed_time = toc;
preprocessing_time = elapsed_time/number_of_iterations;

%% Irregular Allan variance: FAVAR Algorithm
favar = 0;
tic
for k = 1:number_of_iterations % BEGIN: for loop for FAVAR iterations
favar_temp = fcn_AVAR_favarI(input_data(:,k), input_weights(:,k), ...
    list_of_correlation_intervals); % calculate allan variance of irregularly sampled data using FAVAR
favar = favar+favar_temp;
end % END: for loop for FAVAR iterations
elapsed_time = toc;
favar = favar/number_of_iterations;
favar_time = elapsed_time/number_of_iterations;

%% Irregular Allan variance: Normal Algorithm
min_time = time_vector(1);
max_time = time_vector(end);
avar = 0;
tic
for k = 1:number_of_iterations % BEGIN: for loop for AVAR iterations
avar_temp = fcn_AVAR_avarI(ireg_noise_signal(:,k), ireg_time_vector(:,k), ...
    list_of_correlation_time, sampling_interval, min_time, max_time); % calculate allan variance using normal algorithm
avar = avar+avar_temp;
end % END: for loop for AVAR iterations
elapsed_time = toc;
avar = avar/number_of_iterations;
avar_time = elapsed_time/number_of_iterations;

favar_matrix(1:(p-1),p-1) = favar;
avar_matrix(1:(p-1),p-1)  = avar;
elapsed_time_matrix(1,p-1)= favar_time;
elapsed_time_matrix(2,p-1)= avar_time;
elapsed_time_matrix(3,p-1)= preprocessing_time;
end

%% Plots
%% Computational accuracy
figure(12345)
plot(list_of_correlation_intervals,avar,'ko','Markersize',13)
hold on
plot(list_of_correlation_intervals,favar,'k.','Markersize',13)
grid on
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('Correlation interval, $m$','Interpreter','Latex','FontSize',13)
ylabel('AVAR','Interpreter','Latex','FontSize',13)
title('FAVAR-I vs Standard','Interpreter','Latex','FontSize',13)
legend('Normal','FAVAR-I','Interpreter','Latex','FontSize',13,'Location','best')
% Get handle to current axes.
ax = gca;
% Set x and y font sizes.
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;

%% Ratio of wall time
figure(12346)
plot((2.^(2:18)+1), elapsed_time_matrix(2,:)./elapsed_time_matrix(1,:), ...
    'k--','LineWidth',1.2)
hold on
plot((2.^(2:18)+1), elapsed_time_matrix(2,:)./(elapsed_time_matrix(1,:)+elapsed_time_matrix(3,:)), ...
    'r-.','LineWidth',1.2)
yline(1,'k','LineWidth',1.2);
grid on
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('Data length, $N$','Interpreter','Latex','FontSize',13)
ylabel('Ratio of wall time','Interpreter','Latex','FontSize',13)
title('FAVAR-I vs Standard','Interpreter','Latex','FontSize',13)
legend('w/o pre-processing','w/ pre-processing','Interpreter',...
    'Latex','FontSize',13,'Location','best')
% Get handle to current axes.
ax = gca;
% Set x and y font sizes.
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;
