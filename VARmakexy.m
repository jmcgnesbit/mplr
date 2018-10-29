function [Y, X] = VARmakexy(data, nlag, constant)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VARmakexy Produce matrices for VAR estiamation
    %
    % Inputs
    % data
    % nlag - number of lags
    % constant - whether to include constant. 0 = no trend; 1 = constant; 
    %   2 = linear trend; 3 = linear + quadtratic
    %
    % Outputs
    % Y - Y matrix
    % X - X matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Get dimesion of Y
    [nobs, nvars] = size(data);

    % Y matrix 
    Y = data(nlag + 1:end, :);

    % Lags
    X = zeros(nobs - nlag, nvars * nlag);
    for jj = 1 : nlag
            X(:, nvars * (jj - 1) + 1: nvars * jj) = data(nlag - jj + 1 : nobs - jj , :);
    end

    % constant
    if constant == 1 
        X = [ones(nobs-nlag,1) X];

    % constant, linear trend
    elseif constant == 2 
        trend = 1:nobs - nlag;
        X = [ones(nobs-nlag,1) trend' X];

    % constant, linear trend, quad trend
    elseif constant == 3
        trend = 1:nobs - nlag;
        X = [ones(nobs-nlag,1) trend' trend'.^2 X];
    end
end