function v_lambda = compute_TLP_weight(beta_lp, beta_var, variance_matrix, h)
%COMPUTE_TLP_WEIGHT LP-share weight from the TLP empirical risk proof.
%
% The proof gives v = 1 - (Sigma_LP - Sigma_COV) / A, where
% A = T*(beta_LP - beta_VAR)^2 and the Sigmas are asymptotic variances.
% The bootstrap covariance matrices in this code are finite-sample
% variances, so the equivalent expression is:
%   v = 1 - (Var_LP - Cov_LP_VAR) / (beta_LP - beta_VAR)^2.

    if nargin < 4
        h = 1;
    end

    bias_sq = (beta_lp - beta_var)^2;
    var_lp = variance_matrix(1, 1, h);
    cov_lp_var = variance_matrix(1, 2, h);
    variance_gain = var_lp - cov_lp_var;

    if ~isfinite(bias_sq) || ~isfinite(variance_gain)
        v_lambda = NaN;
        return;
    end

    if bias_sq <= eps
        v_lambda = double(variance_gain <= 0);
    else
        v_lambda = 1 - variance_gain / bias_sq;
    end

    v_lambda = max(0, min(1, v_lambda));
end
