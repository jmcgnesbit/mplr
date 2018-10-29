function orth = OrthNorm(n)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Orthnorm Generate an orthognal using QR decomposition such that R has
    %   positive diagonal
    % 
    % Inputs 
    % n - size of matrix
    %
    % Outputs 
    % orth
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    X = randn(n,n);

    % QR decomposition of X
    [Q, R] = qr(X);

    % Check precision
    if sum(sum(Q*Q'))> n + 1.0e-14
        error('Q*transpose(Q) is not equal to identity')
    end

    % Normalization: if diag(R)_i is negative, reverse sign of i'th column of Q diagonal 
    for ii = 1:n
        if R(ii, ii) < 0
            Q(:, ii) = -Q(:, ii);
        end
    end

    % Random orthonormal matrix such that Q*Q'=I
    orth = Q; 

end