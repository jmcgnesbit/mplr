function zero = epsdelta(d, epsilon_equals_delta, M)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % epsdelta Function for fzero. Solve for iso-draw curve intersection
    %   with 45 degree line (i.e. fix M and d, solve for epsilon = delta
    %   such that upper bound of theorem 3 holds)
    %
    % Inputs
    % dim - dimension of parameter region
    % delta 
    % M 
    %
    % Outputs
    % zero - want to be equal to zero
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    draws = optimaldraws(d, epsilon_equals_delta, epsilon_equals_delta);
    
    zero = draws - M;
end