function white_noise = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
                       sampling_frequency,mean_white_noise,...
                       number_of_time_steps,varargin)
%% fcn_AVAR_generateWhiteNoise
%   This function generates white noise characterized by
%   'power_spectral_density'.
%
% FORMAT:
%
%   white_noise = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
%                 sampling_frequency,mean_white_noise,...
%                 number_of_time_steps)
%
% INPUTS:
%
%   power_spectral_density: Power spectral density of white noise [unit^2 s].
%   sampling_frequency: Sampling frequency of the output [Hz].
%   mean_white_noise: Mean of white noise.
%   number_of_time_steps: Desired length of output.
%   varargin: figure number for debugging.
%
% OUTPUTS:
%
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
    if 4>nargin || 5<nargin
        error('Incorrect number of input arguments')
    end
    
    % Check input type and domain
    fcn_AVAR_checkInputsToFunctions(power_spectral_density,'positive');
    fcn_AVAR_checkInputsToFunctions(sampling_frequency,'positive');
    fcn_AVAR_checkInputsToFunctions(mean_white_noise,'number');
    fcn_AVAR_checkInputsToFunctions(number_of_time_steps,'positive integer');
end

if 5 == nargin
    fig_num = varargin{1};
    flag_do_debug = 1;
elseif 1 == flag_do_debug
    fig = figure;
    fig_num = fig.Number;
end

%% Generate White Noise
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
% variance of white noise
variance_white_noise = power_spectral_density*sampling_frequency;

%% Noise generation: White Noise
white_noise = normrnd(mean_white_noise,sqrt(variance_white_noise),...
                      number_of_time_steps,1); % white noise

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
    time_vector = (1/sampling_frequency)*(0:(number_of_time_steps-1))';
    
    figure(fig_num)
    plot(time_vector,white_noise)
    grid on
    ylabel('Measurement')
    xlabel('Time [s]')
    title('White Noise')
    
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end