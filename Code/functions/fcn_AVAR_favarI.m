function [allan_variance,total_weights] = fcn_AVAR_favarI(data,weights,...
                                          list_of_correlation_intervals,...
                                          varargin)
%% fcn_AVAR_favarI
%   This function computes allan variance of irregularly sampled data 'data' 
%   for all the correlation intervals in 'list_of_correlation_intervals'.
%   It uses a recursive algorithm, inspired from FFT, over simple averages 
%   along correlation intervals.
%
% FORMAT:
%
%   allan_variance = fcn_AVAR_favarI(data,weights,list_of_correlation_intervals)
%
% INPUTS:
%
%   data: A Nx1 vector of data points. It contains weighted average of data
%   in a sampling interval.
%   weights: A Nx1 vector containing weights for the data.
%   list_of_correlation_intervals: A Mx1 vector containing list of 
%   correlation intervals. Correlation intervals must be in increasing
%   order and also power of 2.
%   varargin: figure number for debugging.
%
% OUTPUTS:
%
%   allan_variance: A Mx1 vector containing allan varaince corresponding to 
%   the correlation intervals.
%   total_weights: A Mx1 vector containing total weights or the
%   denominator in the allan variance calculation.
%
% This function was written on 2021_05_15 by Satya Prasad
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
    
    % Check the input type and domain
    fcn_AVAR_checkInputsToFunctions(data,'favar data');
    fcn_AVAR_checkInputsToFunctions(weights,'favar weights');
    % Check the compatibility between 'data' and 'weights' input
    if numel(data)~=numel(weights)
        error('The %s and %s input must be of same length', inputname(1), inputname(2));
    end
    fcn_AVAR_checkInputsToFunctions(list_of_correlation_intervals,'favar interval');
end

if 4 == nargin
    fig_num = varargin{1};
    flag_do_debug = 1;
elseif 1 == flag_do_debug
    fig = figure;
    fig_num = fig.Number;
end

%% Calculate Allan Variance of irregularly sampled data using FAVAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of data points in the INPUT data
number_of_datapoints            = numel(data);
% number of correlation intervals
number_of_correlation_intervals = numel(list_of_correlation_intervals);

% intialize vector of weights with 'weights'
vector_of_weights = weights(1:(number_of_datapoints-1));
vector_of_means   = data(1:(number_of_datapoints-1)); % intialize vector of means with 'data'
length_of_vector_of_means = numel(vector_of_means); % length of vector_of_means

% initialize variable to store allan variance
allan_variance = nan(number_of_correlation_intervals,1);
% initialize variable to store total weights
total_weights  = nan(number_of_correlation_intervals,1);
for i = 1:number_of_correlation_intervals % loop over the list of correlation_intervals
    correlation_interval = list_of_correlation_intervals(i); % correlation_interval
    
    vector_of_sums    = vector_of_weights(1:(length_of_vector_of_means-correlation_interval/2)).* ...
                        vector_of_means(1:(length_of_vector_of_means-correlation_interval/2))+...
                        vector_of_weights((1+correlation_interval/2):length_of_vector_of_means).* ...
                        vector_of_means((1+correlation_interval/2):length_of_vector_of_means);
    vector_of_weights = vector_of_weights(1:(length_of_vector_of_means-correlation_interval/2)) + ...
                        vector_of_weights((1+correlation_interval/2):length_of_vector_of_means);
    temp_vector_of_weights = vector_of_weights;
    % zeros are replaced with one in the temporary varaible to avoid division by zero
    temp_vector_of_weights(0==vector_of_weights) = 1;
    vector_of_means   = vector_of_sums./temp_vector_of_weights; % Recurring step
    length_of_vector_of_means = numel(vector_of_means); % length of vector_of_means
    
    vector_of_weights_front = vector_of_weights((1+correlation_interval):length_of_vector_of_means);
    vector_of_means_front   = vector_of_means((1+correlation_interval):length_of_vector_of_means);
    vector_of_weights_back  = vector_of_weights(1:(length_of_vector_of_means-correlation_interval));
    vector_of_means_back    = vector_of_means(1:(length_of_vector_of_means-correlation_interval));
    
    product_of_weights  = vector_of_weights_front.*vector_of_weights_back;
    difference_of_means = vector_of_means_front-vector_of_means_back;
    allan_variance_sum  = sum(product_of_weights.*(difference_of_means.^2));
    
    total_weights(i)  = sum(product_of_weights);
    allan_variance(i) = 0.5*allan_variance_sum/total_weights(i); % write Allan Variance to the output
end % END: For loop over correlation_intervals

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
    plot(list_of_correlation_intervals,allan_variance,'Linewidth',1.2)
    grid on
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    ylabel('Allan Variance','Interpreter','Latex','FontSize',13)
    xlabel('Correlation Interval [No Units]','Interpreter','Latex','FontSize',13)
    title('FAVAR Algorithm','Interpreter','Latex','FontSize',13)
    
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end