function [random_walk,white_noise] = fcn_AVAR_generateRandomWalk(...
                                     random_walk_coefficient,...
                                     sampling_frequency,...
                                     number_of_time_steps,varargin)
%% fcn_AVAR_generateRandomWalk
%   This function generates random walk noise characterized by 
%   'random_walk_coefficient'. Random walk in measurement needs White noise 
%   in measurement rate.
%
% FORMAT:
%
%   [random_walk,white_noise] = fcn_AVAR_generateRandomWalk(...
%                               random_walk_coefficient,...
%                               sampling_frequency,number_of_time_steps)
%
% INPUTS:
%
%   random_walk_coefficient: Noise coefficient for random walk [unit/sqrt(s)].
%   sampling_frequency: Sampling frequency of the output [Hz].
%   number_of_time_steps: Desired length of output.
%   varargin: figure number for debugging.
%
% OUTPUTS:
%
%   random_walk: A 'number_of_time_steps x 1' vector of random walk.
%   white_noise: A 'number_of_time_steps x 1' vector of white noise.
%
% This function was written on 2021_05_14 by Satya Prasad
% Questions or comments? szm888@psu.edu
%

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 0; % Flag to perform input checking

st = dbstack; %#ok<*UNRCH>
if flag_do_debug
    fprintf(1, 'STARTING function: %s, in file: %s\n', st(1).name, st(1).file);
end

%% Check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_check_inputs
    % Are there the right number of inputs?
    if 3>nargin || 4<nargin
        error('Incorrect number of input arguments')
    end
    
    % Check input type and domain
    fcn_AVAR_checkInputsToFunctions(random_walk_coefficient,'positive');
    fcn_AVAR_checkInputsToFunctions(sampling_frequency,'positive');
    fcn_AVAR_checkInputsToFunctions(number_of_time_steps,'positive integer');
end

if 4 == nargin
    fig_num = varargin{1};
    flag_do_debug = 1;
elseif 1 == flag_do_debug
    fig = figure;
    fig_num = fig.Number;
end

%% Generate Random Walk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Derived parameters
power_spectral_density = random_walk_coefficient^2;
% variance of white noise
variance_white_noise   = power_spectral_density*sampling_frequency;
sampling_interval      = 1/sampling_frequency; % [seconds]

%% Noise generation: Random Walk
mean_white_noise = 0; % default value
white_noise = normrnd(mean_white_noise,sqrt(variance_white_noise),...
                      number_of_time_steps,1); % white noise
random_walk = cumsum([0; white_noise(1:number_of_time_steps-1)]...
                     *sampling_interval); % random walk

%% Any debugging?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_debug
    time_vector = sampling_interval*(0:(number_of_time_steps-1))';
    
    figure(fig_num)
    subplot(2,1,1)
    plot(time_vector,random_walk)
    grid on
    ylabel('Measurement')
    title('Random Walk')
    
    subplot(2,1,2)
    plot(time_vector,white_noise)
    grid on
    xlabel('Time [s]')
    ylabel('Measurement Rate')
    
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end