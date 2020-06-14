function epsilon = isodraw_calc(d, delta, M)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % isodraw_calc Function to compute iso-draw curves. Upper bound in
    %   theorem 3 rearranged for epsilon
    %
    % Inputs
    % dim - dimension of parameter region
    % delta 
    % M 
    %
    % Outputs 
    % epsilon
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    old_epsilon = ((2*d) / M) * log((2*d)/delta);
    new_epsilon = (exp(1)/M) * (2*d + log(1/delta));
    epsilon = min(old_epsilon, new_epsilon);
end