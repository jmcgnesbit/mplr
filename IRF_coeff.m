function [IRF_C, LR_C] = IRF_coeff(A, cholsig, nlag, nhorizon, constant)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % irf_coeff Generate IRF's that haven't been multiplied by rotation
    % matrices. This is precomputed to save time
    %
    % Inputs
    % A - matrix of OLS coefficients
    % cholsig - lower triangular cholesky decomp of Sigma
    % nlag - number of lags
    % nhorizon - number of horizons to compute
    % constant - logical for whether a constant has been included
    %
    % Outputs
    % IRF_C - IRF matrices
    % LR_C - Long run IRF matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nvar = size(A, 2);
 
    % Split up A
    A_S = zeros(nvar, nvar, nlag);
    for ii = 1:nlag
        % TRANSPOSE!!
        A_S(:, :, ii) = A(constant + 1 + (ii - 1) * nvar: constant + (ii) * nvar, :)';
    end
 
    % Pre-allocate and initialize
    C = zeros(nvar, nvar, nhorizon);
    C(:, :, 1) = eye(nvar);
 
    % Recursive computation
    for k = 2:nhorizon
        for m = 1:k - 1
            if m > nlag
                A_mat = zeros(nvar, nvar);
            else
                A_mat = A_S(:, :, m);
            end
            C(:, :, k) = C(:, :, k) + C(:, :, k - m) * A_mat;
        end
    end
 
    % Long run
    A_sum = zeros(nvar, nvar);
    for k = 1:nlag
        A_sum = A_sum + A_S(:, :, k);
    end
    LR_C = inv(eye(nvar) - A_sum);
 
    % Times by cholesky
    for jj = 1 : nhorizon
        IRF_C(:, :, jj) = C(:, :, jj) * cholsig;
    end
    IRF_LR = LR_C * cholsig;
 
end


