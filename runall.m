%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to run everything in (Machine) Learning Parameter Regions (2020)
% José Luis Montiel Olea and James Nesbit Email: jmn425@nyu.edu
% 
% Code is broken into three tasks
%
% Task 1 replicates section 3.1 of the paper - Summarizing the identified
% set in set-identified SVARs
%
% Task 3 replicates section 3.2 of the paper - Summarizing a Wald Ellipse
% in point-identified SVARs
% 
% Task 1 replicates section 3.3 of the paper - Highest posterior density
% credible set in SVARs
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
nhorizon = 17; %number of IRF horizons to plot
              
rng(123)
%% Task 1 - Section 3.1 - Summarizing the identified set in set-identified
t1_epsilon = 0.1; 
t1_delta = 0.1;

% Construct a restriction cell By columns {variable sign horizon} Can also
% set horizons to inf Set empty rows with variable 0
t1_restrictionCell(:,:,1) = {1   -1   0;    % GDPDEF negative on impact
                             2   -1   0;    % GDP negative on impact
                             3    1   0};   % DFF positive on impact

% Matrix for number of draws for each task. Each row is a different number
% of draws. First column is number of draws, second is whether this is
% total draws = 0; or draws from inside the set = 1.
t1_drawMatrix = [100 1;
                 optimaldraws(nhorizon, t1_epsilon, t1_delta) 1]; 
 
% Instance of config class
t1_VAR_config = SVAR_config(data, names, constant, nlag, nhorizon, ...
                t1_restrictionCell, t1_epsilon, t1_delta, t1_drawMatrix);

% Analytical bounds from Gafarov, Meier and Montiel Olea(2018)
[ana_upper, ana_lower] = analyticalbounds(data, names, t1_restrictionCell, nlag); 

% Bounds using RRWZ algorithm
[t1_IRF, t1_collector, C_check] = task1(t1_VAR_config); 

% C_check is the 3*3*nhorizon array that stores the C_k IRF coefficients
% defined in footnote 23 in Section 3
% The SVAR function for shock i is injective if there are at least three 
% linearly independent vectors of the ith row of C at each horizon
% We check the first 3 horizons, and we report success in section 3.1

r = zeros(3,1);
for i = 1:3
    r(i) = rank(squeeze(C_check(i,:,1:3)));
end
r % the first three horizons have rank 3, hence SVAR function is injective

% Plot response of deflator to MP shock (Figure 3)
task1_plot(t1_IRF, t1_collector, 1, 1, t1_VAR_config, false, picture_config);

% Plot - Comparsion with analytical bands (Figure 4)
rel_haus = comparison_plot(t1_IRF, ana_lower, ana_upper, t1_collector, 1, 1, ...
    t1_VAR_config, false, picture_config);
% rel_haus is the relative hausdorf distance computed for each number of 
% draws in t1_drawMatrix. The distance in the second element of rel_haus 
% is reported in section 3.1

 
% Empirical CDF's (Figure 5)
f = zeros(size(t1_IRF,4) + 1, nhorizon);
x = zeros(size(t1_IRF,4) + 1, nhorizon);
for i= 1:nhorizon
    % Compute empirical CDF
    [f(:,i), dom] = ecdf(squeeze(t1_IRF(1,1,i,:)));
    % Normalize domain between 0 and 1
    dom = (dom - min(dom)) /(max(dom) - min(dom)) ;
    x(:,i) = dom;
    subplot(6,3,i)
    plot(x(:,i), f(:,i), 'Color','blue')
    hold on;
    plot(x(:,i), x(:,i), '--', 'Color','black')
    title(['$h = ' num2str(i - 1) '$'], 'Interpreter', 'latex', 'fontsize', 8)
end
latex_fig(12, 5, 6);
tightfig()
print(gcf, '-depsc2', fullfile(pic_dir, ['uniform.eps']))


%% Task 2 - Section 3.2 - Summarizing a Wald Ellipse in point-identified
t2_epsilon = 0.1; 
t2_delta = 0.1;
% 68% confidence interval
t2_alpha = 0.68;
% Draw matrix
t2_drawMatrix = [100 0; 
                 2000 0];

t2_VAR_config = SVAR_config(data, names, constant, nlag, nhorizon, ...
    t1_restrictionCell, t2_epsilon, t2_delta, t2_drawMatrix);

% Compute task 2
[t2_IRF_og, t2_IRF, t2_indexmatrix, t2_collector] = task2(t2_VAR_config, ...
    t2_alpha);

% Plot response of variable deflator to MP shock (Figure 6)
% We pass shock 3 (instead of shock 1 in the set-id case), because with
% the cholesky ordering the MP shock is the third shock
task2_plot(t2_IRF_og, t2_IRF, t2_indexmatrix, t2_collector, 1, 3, ...
    t2_VAR_config, false, picture_config, false);


%% Task 3 - Section 3.3 - Highest posterior density credible set in SVARs
t3_epsilon = 0.1; 
t3_delta = 0.1;

t3_alpha = 0.68;
t3_drawMatrix = [100 0; 
                 2000 0];

t3_VAR_config = SVAR_config(data, names, constant, nlag, nhorizon, ...
    t1_restrictionCell, t3_epsilon, t3_delta, t3_drawMatrix);

% Compute task 3
[t3_IRF, t3_indexmatrix, t3_collector] = task3(t3_VAR_config, t3_alpha);

% Plot response of deflator to MP shock (Figure 7)
% We pass shock 3 (instead of shock 1 in the set-id case), because with
% the cholesky ordering the MP shock is the third shock
task3_plot(t3_IRF, t3_indexmatrix, t3_collector, 1, 3, t3_VAR_config, ...
    false, picture_config, false);

%% Isodraw curves - Figure 2
iso_config = pic_config(fontsize, 3, 4, pic_dir);

% Different numbers of draws
Mchoices = [500; 1000; 2500; 5000; 10000];
isodraw_plot(25, Mchoices, iso_config, 'isodrawall', false);


 %% Individual isocurve - Figure 7
 
isodraw_plot(nhorizon, 1360, iso_config, 'isodraw_single');

%% Optimal draws when d = 25, eps = delta = 0.01 (section 2.3)
[up, low] = optimaldraws(25, 0.01, 0.01)

% When d = 17, eps= delta = 0.01
%[up, low] = optimaldraws(17, 0.1, 0.1)
%% Compute where eps = delta (Introduction)

equalepsdelta = fzero(@(delta) epsdelta(nhorizon, delta, 1360),0.1);

display(equalepsdelta);

1 - equalepsdelta



%% Individual isocurve (text) - Figure 11

isodraw_plot(1, 120, iso_config, 'isodraw_text');
equalepsdelta = fzero(@(delta) epsdelta(1, delta, 120),0.1);
equalepsdelta % presented in appendix C
1 - equalepsdelta