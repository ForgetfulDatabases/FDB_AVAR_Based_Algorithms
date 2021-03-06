function [ireg_random_walk,ireg_time_vector,reg_random_walk] = ...
    fcn_AVAR_generateIrregularRandomWalk(random_walk_coefficient,...
    sampling_frequency,number_of_time_steps,upsampling_factor,...
    ireg_time_vector,varargin)
%% fcn_AVAR_generateIrregularRandomWalk
%   This function generates irregularly sampled random walk noise 
%   characterized by random_walk_coefficient. Irregularly sampled data is 
%   generated by randomly down-sampling the regularly-upsampled signal.
%
%   A regularly sampled signal is at a frequency of
%   'upsampling_factor*sampling_frequency'. Later it is downsampled both
%   regularly and irregularly to a frequency of 'sampling_frequency'.
% FORMAT:
%
%   [ireg_random_walk,ireg_time_vector,reg_random_walk] = ...
%       fcn_AVAR_generateIrregularRandomWalk(random_walk_coefficient,...
%       sampling_frequency,number_of_time_steps,upsampling_factor,...
%       ireg_time_vector)
%
% INPUTS:
%
%   random_walk_coefficient: Noise coefficient for random walk [unit/sqrt(s)].
%   sampling_frequency: Target sampling frequency [Hz].
%   number_of_time_steps: Desired length of output.
%   upsampling_factor: Factor to upsampled actual random walk.
%   ireg_time_vector: A 'number_of_time_steps x 1' vector of time 
%   corresponding to 'ireg_random_walk'. It can be NaN.
%   varargin: figure number for debugging.
%
% OUTPUTS:
%
%   ireg_random_walk: A 'number_of_time_steps x 1' vector of irregularly 
%   sampled random walk.
%   ireg_time_vector: A 'number_of_time_steps x 1' vector of time 
%   corresponding to 'ireg_random_walk'.
%   reg_random_walk: A 'number_of_time_steps x 1' vector of regularly 
%   sampled random walk.
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
    if 5>nargin || 6<nargin
        error('Incorrect number of input arguments')
    end
    
    % Check input type and domain
    fcn_AVAR_checkInputsToFunctions(random_walk_coefficient,'positive');
    fcn_AVAR_checkInputsToFunctions(sampling_frequency,'positive');
    fcn_AVAR_checkInputsToFunctions(number_of_time_steps,'positive integer');
    fcn_AVAR_checkInputsToFunctions(upsampling_factor,'positive integer');
    if ~isnan(ireg_time_vector)
        fcn_AVAR_checkInputsToFunctions(ireg_time_vector,'time vector');
    end
end

if 6 == nargin
    fig_num = varargin{1};
    flag_do_debug = 1;
elseif 1 == flag_do_debug
    fig = figure;
    fig_num = fig.Number;
end

%% Generate irregularly sampled Random Walk
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
% actual length of random walk generated
actual_number_of_time_steps = upsampling_factor*number_of_time_steps;
% actual frequency at which random walk is generated
actual_frequency     = upsampling_factor*sampling_frequency;
% variance of white noise
variance_white_noise = power_spectral_density*actual_frequency;
sampling_interval    = 1/actual_frequency; % [seconds]

%% Noise generation: Random Walk
mean_white_noise = 0; % default value
% white noise in measurement rate
measurement_rate = normrnd(mean_white_noise,sqrt(variance_white_noise),...
                           actual_number_of_time_steps,1);
random_walk = cumsum([0; measurement_rate(1:actual_number_of_time_steps-1)]...
                     *sampling_interval); % random walk

if isnan(ireg_time_vector)
    indices_of_irregularly_sampled_data = ...
        sort(randperm(actual_number_of_time_steps-upsampling_factor+1,...
                      number_of_time_steps))';
    % time at which data is sampled
    ireg_time_vector = sampling_interval*...
        (indices_of_irregularly_sampled_data-1);
else
    indices_of_irregularly_sampled_data = ...
        round(actual_frequency*ireg_time_vector)+1;
end
% irregularly sampled random walk
ireg_random_walk = random_walk(indices_of_irregularly_sampled_data);
% regularly sampled random walk
reg_random_walk = random_walk(1:upsampling_factor:actual_number_of_time_steps);

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
    figure(fig_num)
    plot(ireg_time_vector,ireg_random_walk)
    grid on
    xlabel('Time [s]')
    ylabel('Measurement')
    title('Random Walk')
    
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end