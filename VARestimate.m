function [A_OLS, sig_OLS, uhat] = VARestimate(Y, X)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VARestimate Estimate VAR using OLS
    %
    % Inputs
    % Y - Y matrix
    % X - X matrix
    %
    % Outputs
    % A_OLS - OLS estimates
    % Sigma - covariance matrix
    % uhat - residuals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [nobs, nvars] = size(Y);
    
    % Create block diagonal matrix
    A_OLS = (X' * X) \ (X' * Y);
    
    
    % Create residuals
    uhat = Y - X * A_OLS; 
    SSE = uhat' * uhat;
    
    % Sigma
    %sig_OLS = SSE / (nobs - nvars);
    sig_OLS = SSE / (nobs);
end