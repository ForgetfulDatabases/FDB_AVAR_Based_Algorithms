function allan_variance = fcn_AVAR_favar(data,list_of_correlation_intervals,...
                          varargin)
%% fcn_AVAR_favar
%   This function computes allan variance of regularly sampled data 'data'
%   for all the correlation intervals in 'list_of_correlation_intervals'.
%   It uses a recursive algorithm, inspired from FFT, over simple averages 
%   along correlation intervals.
%
% FORMAT:
%
%   allan_variance = fcn_AVAR_favar(data,list_of_correlation_intervals)
%
% INPUTS:
%
%   data: A Nx1 vector of data points. N should be of form 2^p+1 (p >= 2).
%   list_of_correlation_intervals: A Mx1 vector containing list of 
%   correlation intervals. Each interval must be of the form 2^p (p >= 1).
%   varargin: figure number for debugging.
%
% OUTPUTS:
%
%   allan_variance: A Mx1 vector containing allan variance corresponding to 
%   the correlation intervals.
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
    if 2>nargin || 3<nargin
        error('Incorrect number of input arguments')
    end
    
    % Check the input type and domain
    fcn_AVAR_checkInputsToFunctions(data,'favar data');
    fcn_AVAR_checkInputsToFunctions(list_of_correlation_intervals,'favar interval');
end

if 3 == nargin
    fig_num = varargin{1};
    flag_do_debug = 1;
elseif 1 == flag_do_debug
    fig = figure;
    fig_num = fig.Number;
end

%% Calculate Allan Variance using FAVAR
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

vector_of_means = data(1:number_of_datapoints-1); % intialize vector of means with the data
length_of_vector_of_means = numel(vector_of_means); % length of vector_of_means

% initialize variable to store allan variance
allan_variance = nan(number_of_correlation_intervals, 1);
for i = 1:number_of_correlation_intervals  % loop over the list of correlation_intervals
    correlation_interval = list_of_correlation_intervals(i); % correlation_interval
    
    % Recurring step
    vector_of_means = 0.5*(vector_of_means(1:(length_of_vector_of_means-correlation_interval/2))+...
                           vector_of_means((1+correlation_interval/2):length_of_vector_of_means));
    length_of_vector_of_means = numel(vector_of_means); % length of vector_of_means
    
    vector_of_means_front = vector_of_means((1+correlation_interval):length_of_vector_of_means);
    vector_of_means_back  = vector_of_means(1:(length_of_vector_of_means-correlation_interval));
    
    allan_variance_sum = sum((vector_of_means_front-vector_of_means_back).^2);
    
    % write Allan Variance to the output
    allan_variance(i) = 0.5*allan_variance_sum/...
                        (number_of_datapoints-2*correlation_interval);
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