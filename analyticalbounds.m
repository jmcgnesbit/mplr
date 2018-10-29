function [ana_upper, ana_lower] = analyticalbounds(data, names, restrictionCell, nlag)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ANALYTICALBOUNDS Quick comparsion to Gafarov, Meier and Montiel
    %   Olea(2018) using their code (requires analytical_bounds folder)
    %
    % Inputs
    % data - Y matrix
    % names - string vector of variable names
    % restrictionCell - cell of restrictions
    % nlag - number of lags
    %
    % Output
    % ana_upper - upper bounds
    % ana_lower - lower bounds
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % A
    %R = t1_restrictionMatrix;
    %              TS  ho                       sign                  cum
    restMat = [1 2 * restrictionCell{1, 3} restrictionCell{1, 2} 0; % GDPDEF
    2 2 * restrictionCell{2, 3} restrictionCell{2, 2} 0; % GDP
    3 2 * restrictionCell{3, 3} restrictionCell{3, 2} 0]; % DFF
 
    UMPbaselineConfiguration = SVARconfiguration;
    UMPbaselineConfiguration.isCumulativeIRF = 'no';
    UMPbaselineConfiguration.SVARlabel = 'UMPbaseline';
    UMPbaselineConfiguration.nLags = nlag;
    UMPbaselineConfiguration.nLagsMax = nlag;
 
    nTimeSeries = size(names, 2);
    UMPlabelsOfTimeSeries = names; % variable name tags
    unitsOfMeasurement = repmat({'% change'}, 1, nTimeSeries);
 
    UMPTSDiscription = [names; unitsOfMeasurement];
    UMPtsInColumns = data;
 
    dataset = multivariateTimeSeries(UMPtsInColumns, UMPTSDiscription);
    varEstimates = estimatedVecAR(UMPbaselineConfiguration, dataset);
 
    shockName = 'SR_MP';
 
    IDscheme = IDassumptions(restMat, shockName);
 
    umpSVAR = SVAR(varEstimates, IDscheme);
 
    % Save Theta and Sigma for later
    Theta = umpSVAR.getTheta;
    Sigma = umpSVAR.getSigma;
 
    IRFidSet = estimatedIdentifiedSet(umpSVAR);
 
    for i = 1 : nTimeSeries
        ana_upper(:, i) = IRFidSet(1, i).Values;
        ana_lower(:, i) = IRFidSet(2, i).Values;
    end
 
end


