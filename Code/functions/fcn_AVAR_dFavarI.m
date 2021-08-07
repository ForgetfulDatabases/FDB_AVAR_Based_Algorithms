function dynamic_allan_variance = fcn_AVAR_dFavarI(data,weights,...
                                  list_of_correlation_intervals)
% fcn_AVAR_dFavarI
%   This function computes dynamic allan variance of irregularly sampled 
%   data 'data' for all the correlation intervals in 
%   'list_of_correlation_intervals'
%
% FORMAT:
%
%   dynamic_allan_variance = fcn_AVAR_dFavarI(data,weights,...
%                            list_of_correlation_intervals)
%
% INPUTS:
%
%   data: A Nx1 vector of data points. It contains weighted average of data
%   in a sampling interval.
%   weights: A Nx1 vector containing weights for the 'data'.
%   list_of_correlation_intervals: A Mx1 vector containing list of 
%   correlation intervals.
%
% OUTPUTS:
%
%   dynamic_allan_variance: A Mx1 vector containing dynamic allan varaince 
%   corresponding to the correlation intervals.
%
% This function was written on 2021_05_15 by Satya Prasad
% Questions or comments? szm888@psu.edu
%

persistent allan_variance total_weights
flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 0; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
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
    if 3~=nargin
        error('Incorrect number of input arguments')
    end
    
    if ~isempty(allan_variance)
        % Check the input type and domain
        fcn_AVAR_checkInputsToFunctions(data,'dfavar data');
        fcn_AVAR_checkInputsToFunctions(weights,'dfavar weights');
        % Check the compatibility between 'data' and 'weights' input
        if numel(data)~=numel(weights)
            error('The data and weights input must be of same length');
        end
        fcn_AVAR_checkInputsToFunctions(list_of_correlation_intervals,'favar interval');
    end
end

%% Calculate Dynamic Allan Variance of irregularly sampled data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(allan_variance)
    % Initialization Step
    [allan_variance, total_weights] = fcn_AVAR_favarI(data,weights,....
        list_of_correlation_intervals); % Calculate allan variance
else
    % Recurring Step
    number_of_datapoints            = numel(data); % number of data points in the INPUT data
    number_of_correlation_intervals = numel(list_of_correlation_intervals); % number of correlation intervals
    
    vector_of_weights = weights(1:number_of_datapoints-2)+weights(2:number_of_datapoints-1);
    temp_vector_of_weights = vector_of_weights; % store weights in a temporary variable to avoid division by zero
    temp_vector_of_weights(0==vector_of_weights) = 1; % replace zero weights with one to avoid division by zero
    
    vector_of_sums = weights(1:number_of_datapoints-2).*data(1:number_of_datapoints-2)+...
                     weights(2:number_of_datapoints-1).*data(2:number_of_datapoints-1);
    
    vector_of_means = vector_of_sums./temp_vector_of_weights;
    length_of_vector_of_means = numel(vector_of_means);
    
    odd_vector_of_weights  = vector_of_weights(1:2:end);
    odd_vector_of_means    = vector_of_means(1:2:end);
    even_vector_of_weights = vector_of_weights(length_of_vector_of_means:-2:2);
    even_vector_of_means   = vector_of_means(length_of_vector_of_means:-2:2);
    for i = 1:number_of_correlation_intervals
        % loop over the list of correlation_intervals
        correlation_interval = list_of_correlation_intervals(i); % correlation_interval
        
        if 2 == correlation_interval
            add_weight_front = even_vector_of_weights(1);
            add_mean_front   = even_vector_of_means(1);
            
            add_weight_back = even_vector_of_weights(2);
            add_mean_back   = even_vector_of_means(2);
            
            sub_weight_back = odd_vector_of_weights(1);
            sub_mean_back   = odd_vector_of_means(1);
            
            sub_weight_front = odd_vector_of_weights(2);
            sub_mean_front   = odd_vector_of_means(2);
            
        else
            add_sum_front    = add_weight_front*add_mean_front+add_weight_back*add_mean_back;
            add_weight_front = add_weight_front+add_weight_back;
            if 0 ~= add_weight_front
                add_mean_front = add_sum_front/add_weight_front;
            else
                add_mean_front = 0;
            end
            
            add_sum_back    = sum(even_vector_of_weights(1+correlation_interval/2:correlation_interval).* ...
                                  even_vector_of_means(1+correlation_interval/2:correlation_interval));
            add_weight_back = sum(even_vector_of_weights(1+correlation_interval/2:correlation_interval));
            if 0 ~= add_weight_back
                add_mean_back = add_sum_back/add_weight_back;
            else
                add_mean_back = 0;
            end
            
            sub_sum_back    = sub_weight_front*sub_mean_front+sub_weight_back*sub_mean_back;
            sub_weight_back = sub_weight_front+sub_weight_back;
            if 0 ~= sub_weight_back
                sub_mean_back = sub_sum_back/sub_weight_back;
            else
                sub_mean_back = 0;
            end
            
            sub_sum_front    = sum(odd_vector_of_weights(1+correlation_interval/2:correlation_interval).* ...
                                   odd_vector_of_means(1+correlation_interval/2:correlation_interval));
            sub_weight_front = sum(odd_vector_of_weights(1+correlation_interval/2:correlation_interval));
            if 0 ~= sub_weight_front
                sub_mean_front = sub_sum_front/sub_weight_front;
            else
                sub_mean_front = 0;
            end
            
        end
        total_weight = total_weights(i)-...
                       sub_weight_front*sub_weight_back+...
                       add_weight_front*add_weight_back;
        change_in_allan_variance = 0.5*(add_weight_front*add_weight_back*(add_mean_front-add_mean_back)^2 - ...
                                        sub_weight_front*sub_weight_back*(sub_mean_front-sub_mean_back)^2);
        
        allan_variance(i) = (total_weights(i)*allan_variance(i)+...
                             change_in_allan_variance)/total_weight;
        total_weights(i)  = total_weight;
    end
end
dynamic_allan_variance = allan_variance;

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
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end