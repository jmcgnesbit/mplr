%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to produce slides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Task 1 (slides 2 and 6)
task1_slides_plot(t1_IRF, t1_collector, 1, 1, t1_VAR_config, false, picture_config);

%% Isodraw (slide 7)

isodraw_plot(nhorizon, 1000, iso_config, 'slides_isodraw_example');

%% Isodrawall (slide 21)

isodraw_plot(25, Mchoices, iso_config, 'slides_isodraw_all', false);

%% Task 2
task2_slides_plot(t2_IRF_og, t2_IRF, t2_indexmatrix, t2_collector, 1, 3, t2_VAR_config, false, picture_config);

% IK isodraw
isodraw_plot(nhorizon, 1360, iso_config, 'slides_isodraw_t23');

%% Task 3
task3_slides_plot(t3_IRF, t3_indexmatrix, t3_collector, 1, 3, t3_VAR_config, false, picture_config);


%% Isodraw 1982
isodraw_plot(nhorizon, 1982, iso_config, 'slides_isodraw_1982');
equalepsdelta = fzero(@(delta) epsdelta(nhorizon, delta, 1982),0.06);