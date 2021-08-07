%% script_compare_fcn_AVAR_favar_AND_fcn_AVAR_avar.m
%  This script compares 'fcn_AVAR_favar' with 'fcn_AVAR_avar'
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
% increase iterations to decrease the error
number_of_iterations    = 1; % number of monte-carlo simulations

favar_matrix = NaN(17,17);
avar_matrix  = NaN(17,17);
elapsed_time_matrix = NaN(2,17);
for p = 2:18
number_of_time_steps = (2^p)+1;

%% Noise generation: White Noise added to Random Walk
time_vector = (1/sampling_frequency)*(0:(number_of_time_steps-1))';

noise_signal = nan(number_of_time_steps,number_of_iterations); % variable to store noise signals
for k = 1:number_of_iterations % BEGIN: for loop for Noise Signal iterations
random_walk = fcn_AVAR_generateRandomWalk(random_walk_coefficient,....
              sampling_frequency,number_of_time_steps); % generate random walk
white_noise = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
              sampling_frequency,mean_white_noise,number_of_time_steps); % generate white noise
noise_signal(:,k) = random_walk+white_noise; % noise signal
end % END: for loop for Noise Signal iterations

%% Calculation of possible correlation intervals
list_of_correlation_intervals = 2.^(1:(p-1))'; % A vector of correlation intervals
list_of_correlation_time      = list_of_correlation_intervals/sampling_frequency; % A vector of correlation time

%% Allan variance: FAVAR Algorithm
favar = 0;
tic
for k = 1:number_of_iterations % BEGIN: for loop for FAVAR iterations
favar_temp = fcn_AVAR_favar(noise_signal(:,k),list_of_correlation_intervals); % calculate allan variance
favar = favar+favar_temp;
end % END: for loop for FAVAR iterations
elapsed_time = toc;
favar_time = elapsed_time/number_of_iterations;
favar = favar/number_of_iterations;

%% Allan variance: Normal Algorithm
avar = 0;
tic
for k = 1:number_of_iterations % BEGIN: for loop for AVAR iterations
avar_temp = fcn_AVAR_avar(noise_signal(:,k),list_of_correlation_intervals); % calculate allan variance
avar = avar+avar_temp;
end % END: for loop for AVAR iterations
elapsed_time = toc;
avar_time = elapsed_time/number_of_iterations;
avar = avar/number_of_iterations;

favar_matrix(1:(p-1),p-1) = favar;
avar_matrix(1:(p-1),p-1)  = avar;
elapsed_time_matrix(1,p-1)= favar_time;
elapsed_time_matrix(2,p-1)= avar_time;
end

%% Plots
%% Computational accuracy
figure(12345)
subplot(2,1,1)
plot(list_of_correlation_intervals,avar,'ko','Markersize',13)
hold on
plot(list_of_correlation_intervals,favar,'k.','Markersize',13)
grid on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Correlation interval, $m$','Interpreter','Latex','FontSize',13)
ylabel('AVAR','Interpreter','Latex','FontSize',13)
title('(a)','VerticalAlignment','bottom','HorizontalAlignment',...
    'center','Interpreter','Latex','FontSize',13)
legend('Normal','FAVAR','Interpreter','Latex','FontSize',13,...
    'Location','best')
% Get handle to current axes.
ax = gca;
% Set x and y font sizes.
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;

subplot(2,1,2)
plot(list_of_correlation_intervals,abs(avar-favar),'k--','LineWidth',1.2)
grid on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Correlation interval, $m$','Interpreter','Latex','FontSize',13)
ylabel('Absolute error, $|e|$','Interpreter','Latex','FontSize',13)
title('(b)','VerticalAlignment','bottom','HorizontalAlignment','center',...
      'Interpreter','Latex','FontSize',13)
% Get handle to current axes.
ax = gca;
% Set x and y font sizes.
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;

%% Wall time ratio
figure(12346)
plot((2.^(2:18)+1),elapsed_time_matrix(2,:)./elapsed_time_matrix(1,:),...
    'k--','LineWidth',1.2)
hold on
yline(1,'k','LineWidth',1.2);
grid on
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('Data length, $N$','Interpreter','Latex','FontSize',13)
ylabel('Ratio of wall time','Interpreter','Latex','FontSize',13)
title('FAVAR vs Standard','Interpreter','Latex','FontSize',13)
% Get handle to current axes.
ax = gca;
% Set x and y font sizes.
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;
