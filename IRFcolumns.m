% Take IRF function and return columns of shocks

function IRF_col = IRFcolumns(IRF)

[nvar, nshock, nhorizon] = size(IRF);


IRF_col = zeros(nvar * nhorizon, nshock);

for jj = 1:nshock
        IRF_col(:, jj) = reshape(IRF(:, jj, :), nvar * nhorizon, 1);
end

end
