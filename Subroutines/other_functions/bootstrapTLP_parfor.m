function output = bootstrapTLP_parfor(input, bootN, str, use_iid_bootstrap)
%==========================================================================
% Function: bootstrapTLP_parfor
% Purpose:
% Single-level bootstrap for TLP estimation with quantile-based confidence
% intervals. This is a simpler alternative to double_bootstrapTLP_parfor
% that does not require nested bootstrap for standard error estimation.
%
% Process:
% 1. Run bootstrap loop for LP and VAR estimates (B iterations, e.g., 2000)
% 2. Calculate variance-covariance matrix from bootstrap distribution
% 3. Calculate v_lambda using the variance-covariance matrix
% 4. Calculate TLP estimates using v_lambda (treating it as fixed)
% 5. For LP: compute analytical standard errors using sigma^2*(X'X)^{-1}
% 6. For LP: compute t-scores using three centering methods
% 7. For TLP/VAR: use quantile-based confidence intervals directly
%
% Key Differences from double_bootstrapTLP_parfor:
% - No second-level bootstrap for SE estimation
% - LP uses analytical SEs for t-score computation
% - TLP and VAR use direct quantile CIs (no t-score inversion)
% - Typically uses more bootstrap iterations (e.g., 2000 vs 200)
%
%==========================================================================

if nargin < 4
    use_iid_bootstrap = 0;
end

% Unpack input
Y      = input.Y;         
X      = input.X;         
A_hat  = input.A_hat;         

[T, k] = size(Y);
p = length(A_hat) / k;
block_length_T = ceil(T^(1/3));

alpha = str.alpha * 0.01;

bootstrap_seed_vector_B1 = [];
if isfield(str, 'bootstrap_seed_vector_B1') && ~isempty(str.bootstrap_seed_vector_B1)
    bootstrap_seed_vector_B1 = str.bootstrap_seed_vector_B1(:);
    if numel(bootstrap_seed_vector_B1) < bootN
        error('bootstrapTLP_parfor:bootstrapSeedVectorTooShort', ...
            'str.bootstrap_seed_vector_B1 must have at least bootN entries.');
    end
    bootstrap_seed_vector_B1 = bootstrap_seed_vector_B1(1:bootN);
end

% Step 1: Original estimates
IRF_original = irf_function(input.A_hat_f, input.B_0_inv2, str);

data_reshaped_lp2 = PrepareLP_VAR(input, str.H_max, str.H_min, str.P_LP, 0, str.which_irf_y, str.which_irf_x);
H_horizons = str.H_max + 1 - str.H_min;

% Compute original LP estimates and analytical standard errors
Beta_LP_original = zeros(H_horizons, 1);
Beta_LP_SE_analytical = zeros(H_horizons, 1);  % Analytical SE: sqrt(sigma^2 * (X'X)^{-1})

for h = 1:H_horizons
    X_lp = data_reshaped_lp2.X{h};
    Y_lp = data_reshaped_lp2.Y{h};
    n_obs = size(X_lp, 1);
    
    X1 = X_lp(:, 1);
    X_mat = [ones(n_obs, 1) X_lp(:, 2:end)];
    
    % Frisch-Waugh-Lovell projection
    Xpr = X1 - X_mat * (X_mat \ X1);
    Ypr = Y_lp - X_mat * (X_mat \ Y_lp);
    
    Beta_LP_original(h) = Xpr \ Ypr;
    
    % Compute analytical standard error: sqrt(sigma^2 / (X'X))
    % After FWL, this simplifies to sqrt(sigma^2 / (Xpr'Xpr))
    residuals = Ypr - Xpr * Beta_LP_original(h);
    sigma2 = (residuals' * residuals) / (n_obs - size(X_mat, 2) - 1);  % df adjustment
    XprXpr_inv = 1 / (Xpr' * Xpr);
    Beta_LP_SE_analytical(h) = sqrt(sigma2 * XprXpr_inv);
end

Beta_VAR_original = reshape(IRF_original(str.which_irf_y, str.which_irf_x, :), H_horizons, 1);
U = Y - X * A_hat;

% Preallocate arrays for bootstrap
Beta_LP_boot = zeros(H_horizons, bootN);
IRF_boot = zeros(k, k, H_horizons, bootN);

Beta_LP_boot(:, 1) = Beta_LP_original;
IRF_boot(:, :, :, 1) = IRF_original;

%==========================================================================
% SINGLE PARFOR LOOP: Bootstrap LP and VAR estimates
%==========================================================================

parfor b = 1:bootN
    if ~isempty(bootstrap_seed_vector_B1)
        rng(bootstrap_seed_vector_B1(b), 'twister');
    end

    % Local copies
    Beta_LP_star = zeros(H_horizons, 1);
    IRF_star = zeros(k, k, H_horizons);

    if b == 1
        Y_b = Y; 
        X_b = X; 
        est_boot_b = input;
    else
        [Y_b, X_b] = generate_bootstrap_sample_single(Y, X, A_hat, U, T, p, k, str, use_iid_bootstrap, block_length_T);
        select = isfinite(X_b);
        Y_b = Y_b(select(:, end), :);
        X_b = X_b(select(:, end), :);
        est_boot_b = estim_VAR(Y_b, str.P_VAR, X_b);
    end

    IRF_star = irf_function(est_boot_b.A_hat_f, est_boot_b.B_0_inv2, str);

    % LP estimation
    data_boot_local = struct();
    data_boot_local.y_t = Y_b;
    data_reshaped_lp = PrepareLP_VAR(data_boot_local, str.H_max, str.H_min, str.P_LP, 0, str.which_irf_y, str.which_irf_x);
    
    for h = 1:H_horizons
        Xh = data_reshaped_lp.X{h};
        Yh = data_reshaped_lp.Y{h};
        X1 = Xh(:, 1);
        X_mat = [ones(size(Xh, 1), 1) Xh(:, 2:end)];
        Xpr = X1 - X_mat * (X_mat \ X1);
        Ypr = Yh - X_mat * (X_mat \ Yh);
        Beta_LP_star(h) = Xpr \ Ypr;
    end

    % Store results
    Beta_LP_boot(:, b) = Beta_LP_star;
    IRF_boot(:, :, :, b) = IRF_star;
end

%==========================================================================
% Calculate variance-covariance matrix from bootstrap distribution
%==========================================================================

% Extract VAR bootstrap estimates
Beta_VAR_boot = zeros(H_horizons, bootN);
for b = 1:bootN
    Beta_VAR_boot(:, b) = reshape(IRF_boot(str.which_irf_y, str.which_irf_x, :, b), H_horizons, 1);
end

% Compute variance-covariance matrix at each horizon
cov_matrix = nan(2, 2, H_horizons);
cor_matrix = nan(2, 2, H_horizons);

for h = 1:H_horizons
    vector = [Beta_LP_boot(h, :); Beta_VAR_boot(h, :)]';
    valid_rows = ~any(isnan(vector), 2);
    if sum(valid_rows) > 10
        cov_matrix(:, :, h) = cov(vector(valid_rows, :));
        cor_matrix(:, :, h) = corr(vector(valid_rows, :));
    end
end

Variance_matrix = cov_matrix;

%==========================================================================
% Calculate v_lambda and TLP estimates
%==========================================================================

v_lambda_original = zeros(H_horizons, 1);
Beta_TLP_original = zeros(H_horizons, 1);

for h = 1:H_horizons
    if ~isnan(Beta_LP_original(h)) && ~isnan(Beta_VAR_original(h))
        v_lambda_original(h) = compute_TLP_weight(Beta_LP_original(h), Beta_VAR_original(h), Variance_matrix, h);

        % Calculate original TLP estimate
        Beta_TLP_original(h) = v_lambda_original(h) * Beta_LP_original(h) + ...
                               (1 - v_lambda_original(h)) * Beta_VAR_original(h);
    else
        v_lambda_original(h) = NaN;
        Beta_TLP_original(h) = NaN;
    end
end

%==========================================================================
% Calculate TLP variance using closed-form formula
%==========================================================================

var_TLP_original = zeros(H_horizons, 1);
Beta_TLP_SE_original = zeros(H_horizons, 1);

for h = 1:H_horizons
    v_lam = v_lambda_original(h);
    var_lp = Variance_matrix(1, 1, h);
    var_var = Variance_matrix(2, 2, h);
    cov_lp_var = Variance_matrix(1, 2, h);
    
    % Var(TLP) = v^2 * Var(LP) + (1-v)^2 * Var(VAR) + 2*v*(1-v)*Cov(LP,VAR)
    var_TLP_original(h) = v_lam^2 * var_lp + (1 - v_lam)^2 * var_var + ...
                          2 * v_lam * (1 - v_lam) * cov_lp_var;
    Beta_TLP_SE_original(h) = sqrt(max(0, var_TLP_original(h)));
end

%==========================================================================
% Calculate TLP bootstrap estimates (treating v_lambda as fixed)
%==========================================================================

Beta_TLP_boot = zeros(H_horizons, bootN);

for b = 1:bootN
    for h = 1:H_horizons
        if ~isnan(Beta_LP_boot(h, b)) && ~isnan(Beta_VAR_boot(h, b)) && ~isnan(v_lambda_original(h))
            Beta_TLP_boot(h, b) = v_lambda_original(h) * Beta_LP_boot(h, b) + ...
                                  (1 - v_lambda_original(h)) * Beta_VAR_boot(h, b);
        else
            Beta_TLP_boot(h, b) = NaN;
        end
    end
end

%==========================================================================
% Calculate t-scores for LP using analytical SEs (three centering methods)
%==========================================================================

% Compute bootstrap means for centering
Beta_LP_boot_mean = nanmean(Beta_LP_boot, 2);

% Preallocate t-scores for LP
t_scores_LP = struct();
t_scores_LP.method1 = zeros(H_horizons, bootN);  % Centered at bootstrap mean
t_scores_LP.method2 = zeros(H_horizons, bootN);  % Centered at VAR original
t_scores_LP.method3 = zeros(H_horizons, bootN);  % Centered at LP original

for b = 1:bootN
    for h = 1:H_horizons
        if Beta_LP_SE_analytical(h) > 1e-10
            % Method 1: Center at bootstrap mean (Cavaliere et al. 2023)
            t_scores_LP.method1(h, b) = (Beta_LP_boot(h, b) - Beta_LP_boot_mean(h)) / Beta_LP_SE_analytical(h);
            
            % Method 2: Center at VAR original (Montiel Olea-Plagborg-Møller 2021)
            t_scores_LP.method2(h, b) = (Beta_LP_boot(h, b) - Beta_VAR_original(h)) / Beta_LP_SE_analytical(h);
            
            % Method 3: Center at LP original (standard)
            t_scores_LP.method3(h, b) = (Beta_LP_boot(h, b) - Beta_LP_original(h)) / Beta_LP_SE_analytical(h);
        else
            t_scores_LP.method1(h, b) = NaN;
            t_scores_LP.method2(h, b) = NaN;
            t_scores_LP.method3(h, b) = NaN;
        end
    end
end

%==========================================================================
% Compute t-score quantiles for LP
%==========================================================================

t_quantiles_LP = struct();
methods = {'method1', 'method2', 'method3'};

for m = 1:length(methods)
    method = methods{m};
    t_quantiles_LP.(method) = zeros(H_horizons, 2);
    
    for h = 1:H_horizons
        valid_t = t_scores_LP.(method)(h, ~isnan(t_scores_LP.(method)(h, :)));
        if length(valid_t) > 10
            t_quantiles_LP.(method)(h, 1) = quantile(valid_t, alpha/2);
            t_quantiles_LP.(method)(h, 2) = quantile(valid_t, 1 - alpha/2);
        else
            t_quantiles_LP.(method)(h, 1) = norminv(alpha/2);
            t_quantiles_LP.(method)(h, 2) = norminv(1 - alpha/2);
        end
    end
end

%==========================================================================
% Compute direct quantile confidence intervals for TLP and VAR
%==========================================================================

Quantile_bands = struct();

% LP quantiles (for comparison, though we'll use t-score based CIs)
Quantile_bands.LP_lower = zeros(H_horizons, 1);
Quantile_bands.LP_upper = zeros(H_horizons, 1);

% TLP quantiles
Quantile_bands.TLP_lower = zeros(H_horizons, 1);
Quantile_bands.TLP_upper = zeros(H_horizons, 1);

% VAR quantiles
Quantile_bands.VAR_lower = zeros(H_horizons, 1);
Quantile_bands.VAR_upper = zeros(H_horizons, 1);

for h = 1:H_horizons
    % LP
    valid_lp = Beta_LP_boot(h, ~isnan(Beta_LP_boot(h, :)));
    if length(valid_lp) > 10
        Quantile_bands.LP_lower(h) = quantile(valid_lp, alpha/2);
        Quantile_bands.LP_upper(h) = quantile(valid_lp, 1 - alpha/2);
    end
    
    % TLP
    valid_tlp = Beta_TLP_boot(h, ~isnan(Beta_TLP_boot(h, :)));
    if length(valid_tlp) > 10
        Quantile_bands.TLP_lower(h) = quantile(valid_tlp, alpha/2);
        Quantile_bands.TLP_upper(h) = quantile(valid_tlp, 1 - alpha/2);
    end
    
    % VAR
    valid_var = Beta_VAR_boot(h, ~isnan(Beta_VAR_boot(h, :)));
    if length(valid_var) > 10
        Quantile_bands.VAR_lower(h) = quantile(valid_var, alpha/2);
        Quantile_bands.VAR_upper(h) = quantile(valid_var, 1 - alpha/2);
    end
end

%==========================================================================
% Compute bias-corrected bands for LP using t-scores (method3 as default)
%==========================================================================

bias_corrected_bands = struct();

% LP: use inverted t-score method
bias_corrected_bands.LP_lower = Beta_LP_original - t_quantiles_LP.method3(:, 2) .* Beta_LP_SE_analytical;
bias_corrected_bands.LP_upper = Beta_LP_original - t_quantiles_LP.method3(:, 1) .* Beta_LP_SE_analytical;

% TLP and VAR: use direct quantiles (no t-score inversion)
bias_corrected_bands.TLP_lower = Quantile_bands.TLP_lower;
bias_corrected_bands.TLP_upper = Quantile_bands.TLP_upper;
bias_corrected_bands.VAR_lower = Quantile_bands.VAR_lower;
bias_corrected_bands.VAR_upper = Quantile_bands.VAR_upper;

%==========================================================================
% Output
%==========================================================================

output.Beta_LP_original = Beta_LP_original;
output.Beta_TLP_original = Beta_TLP_original;
output.Beta_VAR_original = Beta_VAR_original;
output.v_lambda_original = v_lambda_original;
output.Variance_matrix_original = Variance_matrix;

% LP analytical standard errors
output.Beta_LP_SE_analytical = Beta_LP_SE_analytical;
output.Beta_TLP_SE_original = Beta_TLP_SE_original;

% Compute VAR standard error from bootstrap distribution
Beta_VAR_SE_from_boot = zeros(H_horizons, 1);
for h = 1:H_horizons
    valid_var = Beta_VAR_boot(h, ~isnan(Beta_VAR_boot(h, :)));
    if length(valid_var) > 10
        Beta_VAR_SE_from_boot(h) = std(valid_var);
    else
        Beta_VAR_SE_from_boot(h) = sqrt(Variance_matrix(2, 2, h));
    end
end
output.Beta_VAR_SE_from_boot = Beta_VAR_SE_from_boot;

% t-scores for LP (three centering methods)
output.t_scores_LP = t_scores_LP;
output.t_quantiles_LP = t_quantiles_LP;

% Covariance matrices
output.cov_matrix = cov_matrix;
output.cor_matrix = cor_matrix;
output.cov_matrix_avg = cov_matrix;  % Same as cov_matrix for single bootstrap

% Confidence bands
output.bias_corrected_bands = bias_corrected_bands;
output.Quantile_bands = Quantile_bands;

% Bootstrap distributions
output.Beta_LP_boot = Beta_LP_boot;
output.Beta_TLP_boot = Beta_TLP_boot;
output.Beta_VAR_boot = Beta_VAR_boot;

% For compatibility with double bootstrap output structure
% Use t-scores for LP, create placeholder structures for TLP/VAR
output.z_scores_full_LP = t_scores_LP;
output.z_scores_full_TLP = struct('method1', [], 'method2', [], 'method3', []);
output.z_scores_full_VAR = struct('method1', [], 'method2', [], 'method3', []);

output.z_quantiles_LP = t_quantiles_LP;
output.z_quantiles_TLP = struct('method1', [], 'method2', [], 'method3', []);
output.z_quantiles_VAR = struct('method1', [], 'method2', [], 'method3', []);

% SE estimates
output.Beta_LP_SE_boot = repmat(Beta_LP_SE_analytical, 1, bootN);  % Analytical SE repeated
output.Beta_TLP_SE_boot = repmat(Beta_TLP_SE_original, 1, bootN);  % Closed-form SE repeated
output.Beta_VAR_SE_boot = repmat(Beta_VAR_SE_from_boot, 1, bootN);  % Bootstrap SE repeated

end

%==========================================================================
% Helper function: Generate bootstrap sample
%==========================================================================
function [Y_boot, X_boot] = generate_bootstrap_sample_single(Y, X, A_hat, U, T, p, k, str, use_iid_bootstrap, block_length_T)
    if use_iid_bootstrap == 1
        block_length_T = 1;
        number_of_blocks = ceil(T / block_length_T);
        indx_boot2 = randi(number_of_blocks, number_of_blocks, 1);
    else
        number_of_blocks = ceil(T / block_length_T);
        temp = randi(number_of_blocks - 2, number_of_blocks, 1) * block_length_T + 1;
        temp2 = temp + (0:block_length_T-1);
        indx_boot2 = reshape(temp2', block_length_T * number_of_blocks, 1);
        indx_boot2(T+1:end, :) = [];
    end    
    
    U_boot = U(indx_boot2, :);
    U_boot_demean = U_boot - mean(U_boot);

    y_t = [];
    y_t(1:p, :) = U_boot_demean(1:p, :);
    for t = (p+1):T
        aaa = [];
        for i = 1:p
            aaa(i, :) = (A_hat(k*i - k + 1:k*i, :)' * y_t(t-i, :)')';
        end
        y_t(t, :) = sum(aaa, 1) + U_boot_demean(t, :);
    end    

    Y_boot = y_t(:, :); 
    X_boot = lagmatrix(Y_boot, 1:str.P_VAR);
end
