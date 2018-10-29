function [IRF_new] = recursivebootstrap(A, uhat, Y, nlag, nhorizon, constant)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % recursivebootstrap Perform recusrive bootstrap for task 2 - Wald
    %    confidence set
    %
    % Inputs
    % A - OLS coefficients
    % uhat - residuals from OLS
    % Y - Y data matrix
    % nlag - number of lags
    % nhorizon - number of horizons
    % constant - number of constants
    %
    % Outputs
    % IRF_new - new IRF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [t, nvar] = size(Y);
 
    Y_new = zeros(t + nlag, nvar);
    X_new = zeros(t, constant + nvar * nlag);
 
    % First nlag periods
    Y_new(1: nlag, :) = Y(1 : nlag, :);
 
    % Resample residuals
    u_new = uhat(floor(rand(t, 1) * t) + 1, :);
 
    % Generate data recursively
    for ii = 1:t
        if constant == 0
            constantvec = [];
        elseif constant == 1
            constantvec = 1;
        elseif constant == 2
            constantvec = [1 ii ^ 2];
        end
        X_new(ii, :) = [constantvec reshape(Y_new([nlag + ii - 1 : - 1 : ii], :)', 1, nvar * nlag)];
        Y_new(ii + nlag, :) = X_new(ii, :) * A + u_new(ii, :);
    end
 
    % Chop off first nlag periods
    Y_new = Y_new(nlag + 1 : end, :);
 
    % Get estimates
    [A_new, Sigma_new, ~] = VARestimate(Y_new, X_new);
 
    cholSigma_new = chol(Sigma_new, 'lower');
    IRF_new = IRF_coeff(A_new, Sigma_new, nlag, nhorizon, constant);
 
end

