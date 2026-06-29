function output = irf_VAR_function(data, str, P_VARq)
%==========================================================================
% irf_VAR_function - Estimates Impulse Response Functions (IRFs) from a 
% Vector Autoregression (VAR) model and computes confidence intervals.
%
% Inputs:
%   data     : A Txn matrix containing time series data, where T is the 
%              number of observations and n is the number of variables.
%   str      : A structure containing settings and parameters for the IRF:
%              .k                   : Number of variables in the VAR
%              .H_min, .H_max      : Minimum and maximum horizon for IRFs
%              .BootN              : Number of bootstrap replications
%              .alpha              : Significance level (e.g., 10%)
%              .which_irf_var      : Index of the variable of interest for IRF
%              .random_seed_boot   : Seed for bootstrap reproducibility
%   P_VARq   : Lag order (q) of the VAR model
%
% Outputs:
%   output   : A structure with the following fields:
%              .irf_all : A matrix of all impulse responses (flattened)
%                         Dimensions: (k^2 x (H_max - H_min + 1))
%              .irf     : Selected IRF for variable(s) of interest
%                         Dimensions: (length(which_irf_var) x horizons)
%              .CI      : Confidence intervals of the IRFs from bootstrapping
%
% Example usage:
%   result = irf_VAR_function(myData, settings, 4);
%==========================================================================

    % Estimate the VAR model of lag P_VARq
    VAR_est  = estim_VAR(data, P_VARq, 0);

    % Compute IRFs using estimated VAR coefficients
    IRF_est  = irf_function(VAR_est.A_hat_f, VAR_est.B_0_inv2, str);

    % Reshape IRF results for easier manipulation:
    % The result is a (k^2 x horizon) matrix where each column is a vectorized kxk IRF matrix
    irf_estim = reshape(permute(IRF_est, [2 1 3]), str.k^2, str.H_max - str.H_min + 1);

    % Compute confidence intervals using bootstrap
    CI = conf_int_VAR(VAR_est, str.BootN, str, str.alpha, []);

    % Store results in the output structure
    output.irf_all = irf_estim;                              % All IRFs (vectorized)
    output.irf     = irf_estim(str.which_irf_var, :);        % Selected IRFs
    output.CI      = CI;                                     % Confidence intervals
    output.fullVAR = VAR_est;

end