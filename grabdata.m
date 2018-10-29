function [data, names] = grabdata(startdate, enddate, units, frequency, agg_method)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRABDATA - Creates data used in examples
    %
    % Inputs
    % startdate - string; 'YYYY-MM-DD'
    % enddate - string; 'YYYY-MM-DD'
    % units - string; 'lin' = levels
    % frequency - string; 'q' = quarterly, 'm' = monthly
    % agg_method - 'eop' = end of period
    %
    % Outputs
    % data - data matrix Y
    % names - string vector of variable names
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
    names = {'GDPDEF', 'GDP', 'DFF'};
 
    nvar = size(names, 2);
 
    % GDPDEF
    d = getFredData(names{1}, startdate, enddate, units, frequency);
    % Log levels
    data(:, 1) = log(d.Data(:, 2));
 
    % GDP
    d = getFredData(names{2}, startdate, enddate, units, frequency);
    % Log levels
    data(:, 2) = log(d.Data(:, 2));
 
    % DFF
    d = getFredData(names{3}, startdate, enddate, units, frequency, agg_method);
    % Raw FF
    data(:, 3) = d.Data(1:end, 2);
 
end

