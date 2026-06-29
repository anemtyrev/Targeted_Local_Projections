function Z = constructLaggedMatrix(Y, p)
    % constructLaggedMatrix constructs a lagged matrix Z from Y.
    %
    % Inputs:
    %   Y - Observables [T x K] (T: observations, K: variables)
    %   p - Number of lags
    %
    % Outputs:
    %   Z - Lagged matrix [T-p+1 x K*p]

    [T, K] = size(Y);  % T: time points, K: variables
    Z = zeros(T - p + 1, K * p);  % Initialize Z

    % Construct each row of Z
    for t = p:T
        laggedRow = [];
        for lag = 0:(p-1)
            laggedRow = [laggedRow, Y(t - lag, :)];  % Append lags
        end
        Z(t - p + 1, :) = laggedRow;
    end
end
