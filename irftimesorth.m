function IRF = irftimesorth(IRF_coef, orth)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % irftimesorth Multiply IRF matrices with orthogonal matrices across
    % horizons
    %
    % Arguments
    % IRF_coef - IRF_coef - IRF's without being multiplied by orthogonal matrix. Generated
    % by IRF_coeff.m
    % orth - orthongal matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IRF = zeros(size(IRF_coef));
    for jj = 1 : size(IRF_coef, 3)
        IRF(:, :, jj) = IRF_coef(:, :, jj) * orth;
    end
 
end
