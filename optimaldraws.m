function [upperbound, lowerbound] = optimaldraws(d, epsilon, delta)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optimaldraws Compute the optimal number of draws from the inside using
    %   Theorem 1.
    %
    % Inputs
    % d - dimension of parameter regin
    % epsilon
    % delta
    %
    % Outputs 
    % upperbound 
    % lowerbound
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    lowerbound_1 = (1- epsilon)* log(1/delta) /(epsilon);
    lowerbound_2 = (3/16) * (d / epsilon);
    lowerbound = ceil(max(lowerbound_1, lowerbound_2));
    
    upperbound_1 = ((2*d) / epsilon) * log((2*d)/delta);
    upperbound_2 = (exp(1)/epsilon) * (2*d + log(1/delta)); 
    upperbound = ceil(min(upperbound_1, upperbound_2));
end