%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to run everything in (Machine) Learning Parameter Regions (2018)
% José Luis Montiel Olea and James Nesbit Email: jmn425@nyu.edu
% 
% Code is broken into three tasks
%
% Task 1 replicates section 3.1 of the paper - Summarizing the identified
% set in set-identified SVARs
%
% Task 1 replicates section 3.2 of the paper - Summarizing a Wald Ellipse
% in point-identified SVARs
% 
% Task 1 replicates section 3.3 of the paper - Highest posterior density
% credible set in SVARs
%
% Additional produces the iso-draw curves for Figure 2 and Appendix B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

clc
clear
 
% Working directory
addpath(genpath(pwd))

% Save pictures
pic_dir = './Pictures';

%% Settings

% Data settings
startdate = '1982-10-01';
enddate = '2007-10-01';
units = 'lin'; % levels
frequency = 'q';
agg_method = 'eop'; % end of period

% Picture settings
fontsize = 12;
height = 2.5;
width = 6;
picture_config = pic_config(fontsize, height, width, pic_dir);

%% Data
data = grabdata(startdate, enddate, units, frequency, agg_method);

names = {'\ln p_t', 'GDP', 'DFF'};

%% VAR settings
level = .95; % significance level for analytical bands
nlag = 4; % number of lags in estimation
constant = 1; %include a constant 
nhorizon = 17;
epsilon = 0.1;
delta = 0.1;
                         
%% Task 1 - Section 3.1 - Summarizing the identified set in set-identified
% SVARs

% Construct a restriction cell By columns {variable sign horizon} Can also
% set horizons to inf Set empty rows with variable 0
t1_restrictionCell(:,:,1) = {1   -1   0;    % GDPDEF negative on impact
                             2   -1   0;    % GDP negative on impact
                             3    1   0};   % DFF positive on impact

% Matrix for number of draws for each task. Each row is a different number
% of draws. First column is number of draws, second is whether this is
% total draws = 0; or draws from inside the set = 1.
t1_drawMatrix = [100 1;
                 500 1;
                 optimaldraws(nhorizon, epsilon, delta) 1]; % this row is the optimal number of draws for delta and epsilon
 
% Instance of config class
t1_VAR_config = SVAR_config(data, names, constant, nlag, nhorizon, t1_restrictionCell, epsilon, delta, t1_drawMatrix);

% Analytical bounds from Gafarov, Meier and Montiel Olea(2018)
% [ana_upper, ana_lower] = analyticalbounds(data, names, t1_restrictionCell, nlag); 

% Bounds using RRWZ algorithm
[t1_IRF, t1_collector] = task1(t1_VAR_config); 

% Plot - response of variable 1 to shock 1
task1_plot(t1_IRF, t1_collector, 1, 1, t1_VAR_config, false, picture_config);

% Plot - Comparsion with analytical bands
%comparison_plot(t1_IRF, ana_lower, ana_upper, t1_collector, 1, 1, t1_VAR_config, false, picture_config);
 
%% Task 2 - Section 3.2 - Summarizing a Wald Ellipse in point-identified
% SVARs

% 68% confidence interval
t2_alpha = 0.68;
% Draw matrix
t2_drawMatrix = [100 0; 
                 2000 0];

t2_VAR_config = SVAR_config(data, names, constant, nlag, nhorizon, t1_restrictionCell, epsilon, delta, t2_drawMatrix);

% Compute task 2
[t2_IRF_og, t2_IRF, t2_indexmatrix, t2_collector] = task2(t2_VAR_config, t2_alpha);

% Plot (response of variable 1 to shock 1)
task2_plot(t2_IRF_og, t2_IRF, t2_indexmatrix, t2_collector, 1, 3, t2_VAR_config, false, picture_config, false);


%% Task 3 - Section 3.3 - Highest posterior density credible set in SVARs
t3_alpha = 0.68;
t3_drawMatrix = [100 0; 
                 2000 0];

t3_VAR_config = SVAR_config(data, names, constant, nlag, nhorizon, t1_restrictionCell, epsilon, delta, t3_drawMatrix);

[t3_IRF, t3_indexmatrix, t3_collector] = task3(t3_VAR_config, t3_alpha);

task3_plot(t3_IRF, t3_indexmatrix, t3_collector, 1, 3, t3_VAR_config, false, picture_config, false);

%% Isodraw curves - Figure 2
iso_config = pic_config(fontsize, 3, 4, pic_dir);

% Different numbers of draws
Mchoices = [500; 1000; 2500; 5000; 10000];
isodraw_plot(25, Mchoices, iso_config, 'isodrawall', false);


 %% Individual isocurve - Figure 5
 
isodraw_plot(nhorizon, 1360, iso_config, 'isodraw_single');

%% Compute where eps = delta

equalepsdelta = fzero(@(delta) epsdelta(nhorizon, delta, 1360),0.1);

display(equalepsdelta);