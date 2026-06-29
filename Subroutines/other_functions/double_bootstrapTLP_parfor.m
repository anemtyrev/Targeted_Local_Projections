function output = double_bootstrapTLP_parfor(input, bootN, str, use_iid_bootstrap, compute_SLP, use_bootstrap_bias, use_second_level_tlp_variance)
%==========================================================================
% Function: double_bootstrapTLP_parfor
% Purpose:
%   Parallel version that computes full z-score distributions for
%   bias-corrected confidence bands using the Quantile-t method
%
% INPUTS:
%   input            - VAR estimation output
%   bootN            - Number of bootstrap replications
%   str              - Settings structure
%   use_iid_bootstrap - 1 for iid, 0 for block bootstrap
%   compute_SLP      - FLAG: 1 to compute SLP, 0 to skip (default: 0)
%   use_bootstrap_bias - FLAG: 1 to compute TLP weights using bootstrap-layer
%                        bias estimates (Beta_LP_b1 - Beta_VAR_b1), 
%                        0 to use original sample bias (default: 0)
%   use_second_level_tlp_variance - FLAG: 1 to estimate each first-level
%                        TLP SE directly from second-level TLP draws whose
%                        weights are recomputed from second-level LP/VAR
%                        estimates. The weight variance term is still
%                        Variance_matrix_mean, matching first-level weights.
%                        0 keeps the weighted LP/VAR/COV formula (default: 1).
%
% NOTES ON use_bootstrap_bias:
%   When use_bootstrap_bias = 0 (default):
%     - v_lambda is computed once using (Beta_LP_original - Beta_VAR_original)^2
%     - Same weights applied to all bootstrap replications
%   When use_bootstrap_bias = 1:
%     - For each b1, v_lambda_b1 is computed using (Beta_LP_b1 - Beta_VAR_b1)^2
%     - Variance terms still come from Variance_matrix_mean (averaged across all b1)
%     - Each bootstrap replication has its own weight
%
%==========================================================================

if nargin < 4
    use_iid_bootstrap = 0;
end
if nargin < 5
    compute_SLP = 0;
end
if nargin < 6
    use_bootstrap_bias = 0;
end
if nargin < 7
    if isfield(str, 'use_second_level_tlp_variance')
        use_second_level_tlp_variance = str.use_second_level_tlp_variance;
    else
        use_second_level_tlp_variance = 1;
    end
end

Y = input.Y;
X = input.X;
A_hat = input.A_hat;
[T,k] = size(Y);
p = length(A_hat)/k;
block_length_T = ceil(T^(1/3));
alpha = str.alpha*0.01;
B2 = 100;

H_horizons = str.H_max + 1 - str.H_min;

bootstrap_seed_vector_B1 = [];
if isfield(str, 'bootstrap_seed_vector_B1') && ~isempty(str.bootstrap_seed_vector_B1)
    bootstrap_seed_vector_B1 = str.bootstrap_seed_vector_B1(:);
    if numel(bootstrap_seed_vector_B1) < bootN
        error('double_bootstrapTLP_parfor:bootstrapSeedVectorTooShort', ...
            'str.bootstrap_seed_vector_B1 must have at least bootN entries.');
    end
    bootstrap_seed_vector_B1 = bootstrap_seed_vector_B1(1:bootN);
end

%==========================================================================
% Step 0: Estimate lambda_opt for SLP on ORIGINAL data (only if compute_SLP)
%==========================================================================
data_original = struct();
data_original.y_t = input.y_t;

if compute_SLP
    data_reshaped_slp = PrepareLP_VAR(data_original, str.H_max, str.H_min, str.P_LP, 1, str.which_irf_y, str.which_irf_x);
    X_slp_orig = data_reshaped_slp.X;
    Y_slp_orig = data_reshaped_slp.Y;
    idx_slp = data_reshaped_slp.idx;
    Beta_OLS_slp = X_slp_orig \ Y_slp_orig;

    C = Create_C(H_horizons, 1, 0);

    Struct = [];
    Struct.C = C;
    Struct.h = str.H_max;
    Struct.H_min = str.H_min;
    Struct.T = T;
    Struct.X = X_slp_orig;
    Struct.Y = Y_slp_orig;
    Struct.Z = X_slp_orig;
    Struct.idx = idx_slp;
    Struct.theta = Beta_OLS_slp;
    Struct.S = inv(X_slp_orig' * X_slp_orig);
    Struct.nlag = H_horizons;
    Struct.AR = data_original;

    HAC = var_cov_quadratic(Struct, 0, str.alpha);
    Struct.Sigma = HAC.VAR;

    opts = optimset('Display', 'off', 'TolCon', 1E-10, 'TolFun', 1E-10, 'TolX', 1E-10);
    lambda_start = 1000;
    x_L = 0;

    quadratic_method = 2;
    if isfield(str, 'quadratic_method')
        quadratic_method = str.quadratic_method;
    end

    objfun = @(lambda) Select_lambda(lambda, Struct, quadratic_method);
    lambda_opt = fmincon(objfun, lambda_start, [], [], [], [], x_L, [], [], opts);

    XtX_orig = X_slp_orig' * X_slp_orig;
    Beta_SLP_original = (XtX_orig + lambda_opt * C) \ (X_slp_orig' * Y_slp_orig);
else
    lambda_opt = NaN;
    C = [];
    Beta_SLP_original = NaN(H_horizons, 1);
end

%==========================================================================
% Step 1: Original estimates (LP, VAR)
%==========================================================================
IRF_original = irf_function(input.A_hat_f, input.B_0_inv2, str);
data_reshaped_lp2 = PrepareLP_VAR(data_original, str.H_max, str.H_min, str.P_LP, 0, str.which_irf_y, str.which_irf_x);

Beta_LP_original = zeros(H_horizons, 1);
for h = 1:H_horizons
    X_lp = data_reshaped_lp2.X{h};
    Y_lp = data_reshaped_lp2.Y{h};
    X1 = X_lp(:,1);
    X_mat = [ones(size(X_lp,1),1) X_lp(:,2:end)];
    Xpr = X1 - X_mat * (X_mat\X1);
    Ypr = Y_lp - X_mat * (X_mat\Y_lp);
    Beta_LP_original(h) = Xpr\Ypr;
end

Beta_VAR_original = reshape(IRF_original(str.which_irf_y, str.which_irf_x, :), H_horizons, 1);
U = Y - X*A_hat;

% Preallocate
Beta_LP_boot = zeros(H_horizons, bootN);
Beta_SLP_boot = zeros(H_horizons, bootN);
IRF_boot = zeros(k, k, H_horizons, bootN);
Beta_LP_SE_boot = zeros(H_horizons, bootN);
Beta_VAR_SE_boot = zeros(H_horizons, bootN);
Beta_SLP_SE_boot = zeros(H_horizons, bootN);
cov_matrix_boot = nan(2, 2, H_horizons, bootN);
cor_matrix_boot = nan(2, 2, H_horizons, bootN);
Beta_LP_second_level_boot = nan(H_horizons, B2, bootN);
Beta_VAR_second_level_boot = nan(H_horizons, B2, bootN);

Beta_LP_boot(:,1) = Beta_LP_original;
Beta_SLP_boot(:,1) = Beta_SLP_original;
IRF_boot(:,:,:,1) = IRF_original;

%==========================================================================
% MAIN BOOTSTRAP LOOP
%==========================================================================
parfor b = 1:bootN
    if ~isempty(bootstrap_seed_vector_B1)
        rng(bootstrap_seed_vector_B1(b), 'twister');
    end

    Beta_LP_star = zeros(H_horizons,1);
    Beta_SLP_star = zeros(H_horizons,1);
    IRF_star = zeros(k,k,H_horizons);
    S_matrix_b = [];
    
    if b == 1
        Y_b = Y;
        X_b = X;
        est_boot_b = input;
        data_boot = data_original;
    else
        [Y_b, X_b] = generate_bootstrap_sample(Y, X, A_hat, U, T, p, k, str, [], use_iid_bootstrap, block_length_T);
        select = isfinite(X_b);
        Y_b = Y_b(select(:,end),:);
        X_b = X_b(select(:,end),:);
        est_boot_b = estim_VAR(Y_b, str.P_VAR, X_b);
        data_boot = struct();
        data_boot.y_t = Y_b;
    end
    
    IRF_star = irf_function(est_boot_b.A_hat_f, est_boot_b.B_0_inv2, str);
    
    % LP estimation (horizon-by-horizon)
    data_reshaped_lp = PrepareLP_VAR(data_boot, str.H_max, str.H_min, str.P_LP, 0, str.which_irf_y, str.which_irf_x);
    for h = 1:H_horizons
        Xh = data_reshaped_lp.X{h};
        Yh = data_reshaped_lp.Y{h};
        X1 = Xh(:,1);
        X_mat = [ones(size(Xh,1),1) Xh(:,2:end)];
        Xpr = X1 - X_mat * (X_mat\X1);
        Ypr = Yh - X_mat * (X_mat\Yh);
        Beta_LP_star(h) = Xpr\Ypr;
    end
    
    % SLP estimation with FIXED lambda_opt (only if compute_SLP)
    if compute_SLP
        data_reshaped_slp_b = PrepareLP_VAR(data_boot, str.H_max, str.H_min, str.P_LP, 1, str.which_irf_y, str.which_irf_x);
        X_slp_b = data_reshaped_slp_b.X;
        Y_slp_b = data_reshaped_slp_b.Y;
        XtX_b = X_slp_b' * X_slp_b;
        Beta_SLP_star = (XtX_b + lambda_opt * C) \ (X_slp_b' * Y_slp_b);
        S_matrix_b = (XtX_b + lambda_opt * C) \ XtX_b;
    else
        Beta_SLP_star = NaN(H_horizons, 1);
    end
    
    % Second-level bootstrap for SE estimation
    U_star = Y_b - X_b*est_boot_b.A_hat;
    Beta_LP_second_level = zeros(H_horizons, B2);
    Beta_VAR_second_level = zeros(H_horizons, B2);
    
    for b2 = 1:B2
        [Y_b2, X_b2] = generate_bootstrap_sample(Y_b, X_b, est_boot_b.A_hat, U_star, size(Y_b,1), p, k, str, [], use_iid_bootstrap, block_length_T);
        select2 = isfinite(X_b2);
        Y_b2 = Y_b2(select2(:,end),:);
        X_b2 = X_b2(select2(:,end),:);
        
        if size(Y_b2,1) > size(X_b2,2)
            est_boot2 = estim_VAR(Y_b2, str.P_VAR, X_b2);
            IRF_temp = irf_function(est_boot2.A_hat_f, est_boot2.B_0_inv2, str);
            Beta_VAR_second_level(:,b2) = reshape(IRF_temp(str.which_irf_y, str.which_irf_x,:), H_horizons, 1);
            
            data_boot2_local = struct();
            data_boot2_local.y_t = Y_b2;
            
            % LP second level (horizon-by-horizon)
            data_reshaped_lp2_b = PrepareLP_VAR(data_boot2_local, str.H_max, str.H_min, str.P_LP, 0, str.which_irf_y, str.which_irf_x);
            for h = 1:H_horizons
                Xh2 = data_reshaped_lp2_b.X{h};
                Yh2 = data_reshaped_lp2_b.Y{h};
                if size(Xh2,1) > size(Xh2,2)
                    X1_2 = Xh2(:,1);
                    X_mat2 = [ones(size(Xh2,1),1) Xh2(:,2:end)];
                    Xpr2 = X1_2 - X_mat2*(X_mat2\X1_2);
                    Ypr2 = Yh2 - X_mat2*(X_mat2\Yh2);
                    Beta_LP_second_level(h,b2) = Xpr2\Ypr2;
                else
                    Beta_LP_second_level(h,b2) = NaN;
                end
            end
        else
            Beta_VAR_second_level(:,b2) = NaN;
            Beta_LP_second_level(:,b2) = NaN;
        end
    end
    
    % Compute SEs for LP and VAR
    SE_LP = zeros(H_horizons, 1);
    SE_VAR = zeros(H_horizons, 1);
    
    for h = 1:H_horizons
        valid = ~isnan(Beta_LP_second_level(h,:));
        if sum(valid) > 10
            SE_LP(h) = std(Beta_LP_second_level(h,valid));
        else
            SE_LP(h) = NaN;
        end
        
        valid_var = ~isnan(Beta_VAR_second_level(h,:));
        if sum(valid_var) > 10
            SE_VAR(h) = std(Beta_VAR_second_level(h,valid_var));
        else
            SE_VAR(h) = NaN;
        end
        
        vec = [Beta_LP_second_level(h,:); Beta_VAR_second_level(h,:)]';
        valid_rows = ~any(isnan(vec),2);
        if sum(valid_rows) > 10
            cov_matrix_boot(:,:,h,b) = cov(vec(valid_rows,:));
            cor_matrix_boot(:,:,h,b) = corr(vec(valid_rows,:));
        else
            cov_matrix_boot(:,:,h,b) = NaN(2);
            cor_matrix_boot(:,:,h,b) = NaN(2);
        end
    end
    
    % SLP SE via transformation: Var(SLP) = S * Var(LP) * S'
    % Keep the cross-horizon LP covariance from the second-level draws.
    if compute_SLP
        valid_lp_draws = ~any(isnan(Beta_LP_second_level), 1);
        if sum(valid_lp_draws) > 10
            Cov_LP = cov(Beta_LP_second_level(:, valid_lp_draws)');
            Cov_SLP = S_matrix_b * Cov_LP * S_matrix_b';
            SE_SLP = sqrt(max(0, diag(Cov_SLP)));
        else
            SE_SLP = NaN(H_horizons, 1);
        end
    else
        SE_SLP = NaN(H_horizons, 1);
    end
    
    Beta_LP_boot(:,b) = Beta_LP_star;
    Beta_SLP_boot(:,b) = Beta_SLP_star;
    IRF_boot(:,:,:,b) = IRF_star;
    Beta_LP_SE_boot(:,b) = SE_LP;
    Beta_VAR_SE_boot(:,b) = SE_VAR;
    Beta_SLP_SE_boot(:,b) = SE_SLP;
    Beta_LP_second_level_boot(:,:,b) = Beta_LP_second_level;
    Beta_VAR_second_level_boot(:,:,b) = Beta_VAR_second_level;
end

%==========================================================================
% Calculate mean variance matrix
%==========================================================================
Variance_matrix_mean = nanmean(cov_matrix_boot, 4);

%==========================================================================
% Extract VAR bootstrap estimates
%==========================================================================
Beta_VAR_boot = zeros(H_horizons, bootN);
for b = 1:bootN
    Beta_VAR_boot(:,b) = reshape(IRF_boot(str.which_irf_y, str.which_irf_x,:,b), H_horizons, 1);
end

%==========================================================================
% Calculate v_lambda and TLP estimates
% TWO APPROACHES based on use_bootstrap_bias flag
%==========================================================================

if use_bootstrap_bias == 0
    %======================================================================
    % ORIGINAL APPROACH: v_lambda computed once from original sample
    %======================================================================
    v_lambda_original = zeros(H_horizons, 1);
    Beta_TLP_original = zeros(H_horizons, 1);

    for h = 1:H_horizons
        if ~isnan(Beta_LP_original(h)) && ~isnan(Beta_VAR_original(h))
            v_lambda_original(h) = compute_TLP_weight(Beta_LP_original(h), Beta_VAR_original(h), Variance_matrix_mean, h);
            Beta_TLP_original(h) = v_lambda_original(h) * Beta_LP_original(h) + (1-v_lambda_original(h)) * Beta_VAR_original(h);
        else
            v_lambda_original(h) = NaN;
            Beta_TLP_original(h) = NaN;
        end
    end

    % TLP bootstrap estimates with FIXED weights
    Beta_TLP_boot = zeros(H_horizons, bootN);
    Beta_TLP_SE_boot = zeros(H_horizons, bootN);
    Beta_TLP_boot(:,1) = Beta_TLP_original;
    
    % Store v_lambda for each bootstrap (all same in this case)
    v_lambda_boot = repmat(v_lambda_original, 1, bootN);

    for b = 1:bootN
        for h = 1:H_horizons
            if ~isnan(Beta_LP_boot(h,b)) && ~isnan(Beta_VAR_boot(h,b)) && ~isnan(v_lambda_original(h))
                Beta_TLP_boot(h,b) = v_lambda_original(h)*Beta_LP_boot(h,b) + (1-v_lambda_original(h))*Beta_VAR_boot(h,b);
            else
                Beta_TLP_boot(h,b) = NaN;
            end
        end
        
        for h = 1:H_horizons
            if ~isnan(Beta_LP_SE_boot(h,b)) && ~isnan(Beta_VAR_SE_boot(h,b)) && ~isnan(v_lambda_original(h))
                v_lam = v_lambda_original(h);
                var_lp = Beta_LP_SE_boot(h,b)^2;
                var_var = Beta_VAR_SE_boot(h,b)^2;
                if ~isnan(cov_matrix_boot(1,2,h,b))
                    cov_lp_var = cov_matrix_boot(1,2,h,b);
                else
                    cov_lp_var = Variance_matrix_mean(1,2,h);
                end
                var_tlp = v_lam^2*var_lp + (1-v_lam)^2*var_var + 2*v_lam*(1-v_lam)*cov_lp_var;
                Beta_TLP_SE_boot(h,b) = sqrt(max(0, var_tlp));
            else
                Beta_TLP_SE_boot(h,b) = NaN;
            end
        end
    end

else
    %======================================================================
    % NEW APPROACH: v_lambda computed for each bootstrap replication
    % using (Beta_LP_b1 - Beta_VAR_b1)^2 for bias term
    % but Variance_matrix_mean for variance terms
    %======================================================================
    
    % First compute v_lambda for original sample (for output and original TLP)
    v_lambda_original = zeros(H_horizons, 1);
    Beta_TLP_original = zeros(H_horizons, 1);

    for h = 1:H_horizons
        if ~isnan(Beta_LP_original(h)) && ~isnan(Beta_VAR_original(h))
            v_lambda_original(h) = compute_TLP_weight(Beta_LP_original(h), Beta_VAR_original(h), Variance_matrix_mean, h);
            Beta_TLP_original(h) = v_lambda_original(h) * Beta_LP_original(h) + (1-v_lambda_original(h)) * Beta_VAR_original(h);
        else
            v_lambda_original(h) = NaN;
            Beta_TLP_original(h) = NaN;
        end
    end
    
    % Now compute bootstrap-specific weights and TLP estimates
    Beta_TLP_boot = zeros(H_horizons, bootN);
    Beta_TLP_SE_boot = zeros(H_horizons, bootN);
    v_lambda_boot = zeros(H_horizons, bootN);
    
    Beta_TLP_boot(:,1) = Beta_TLP_original;
    v_lambda_boot(:,1) = v_lambda_original;

    for b = 1:bootN
        % Compute v_lambda for this bootstrap replication
        v_lambda_b = zeros(H_horizons, 1);
        
        for h = 1:H_horizons
            if ~isnan(Beta_LP_boot(h,b)) && ~isnan(Beta_VAR_boot(h,b))
                v_lambda_b(h) = compute_TLP_weight(Beta_LP_boot(h,b), Beta_VAR_boot(h,b), Variance_matrix_mean, h);
                
                % Compute TLP for this bootstrap using its own weight
                Beta_TLP_boot(h,b) = v_lambda_b(h)*Beta_LP_boot(h,b) + (1-v_lambda_b(h))*Beta_VAR_boot(h,b);
            else
                v_lambda_b(h) = NaN;
                Beta_TLP_boot(h,b) = NaN;
            end
        end
        
        v_lambda_boot(:,b) = v_lambda_b;
        
        % Compute TLP SE using bootstrap-specific weight
        for h = 1:H_horizons
            if ~isnan(Beta_LP_SE_boot(h,b)) && ~isnan(Beta_VAR_SE_boot(h,b)) && ~isnan(v_lambda_b(h))
                v_lam = v_lambda_b(h);
                var_lp = Beta_LP_SE_boot(h,b)^2;
                var_var = Beta_VAR_SE_boot(h,b)^2;
                if ~isnan(cov_matrix_boot(1,2,h,b))
                    cov_lp_var = cov_matrix_boot(1,2,h,b);
                else
                    cov_lp_var = Variance_matrix_mean(1,2,h);
                end
                var_tlp = v_lam^2*var_lp + (1-v_lam)^2*var_var + 2*v_lam*(1-v_lam)*cov_lp_var;
                Beta_TLP_SE_boot(h,b) = sqrt(max(0, var_tlp));
            else
                Beta_TLP_SE_boot(h,b) = NaN;
            end
        end
    end
end

if use_second_level_tlp_variance == 1
    Beta_TLP_SE_boot_direct = nan(H_horizons, bootN);

    for b = 1:bootN
        for h = 1:H_horizons
            Beta_LP_second_h = Beta_LP_second_level_boot(h,:,b);
            Beta_VAR_second_h = Beta_VAR_second_level_boot(h,:,b);
            Beta_TLP_second_h = nan(1, B2);
            v_lambda_second_h = nan(1, B2);

            for b2 = 1:B2
                if ~isnan(Beta_LP_second_h(b2)) && ~isnan(Beta_VAR_second_h(b2))
                    v_lambda_second_h(b2) = compute_TLP_weight( ...
                        Beta_LP_second_h(b2), Beta_VAR_second_h(b2), Variance_matrix_mean, h);
                    if ~isnan(v_lambda_second_h(b2))
                        Beta_TLP_second_h(b2) = v_lambda_second_h(b2) * Beta_LP_second_h(b2) ...
                            + (1-v_lambda_second_h(b2)) * Beta_VAR_second_h(b2);
                    end
                end
            end

            valid_tlp = ~isnan(Beta_TLP_second_h);
            if sum(valid_tlp) > 10
                Beta_TLP_SE_boot_direct(h,b) = std(Beta_TLP_second_h(valid_tlp));
            end
        end
    end

    Beta_TLP_SE_boot_weighted = Beta_TLP_SE_boot;
    Beta_TLP_SE_boot = Beta_TLP_SE_boot_direct;
else
    Beta_TLP_SE_boot_direct = [];
    Beta_TLP_SE_boot_weighted = Beta_TLP_SE_boot;
end

%==========================================================================
% Z-scores for all methods
%==========================================================================
z_scores_full_LP = struct();
z_scores_full_TLP = struct();
z_scores_full_VAR = struct();
z_scores_full_SLP = struct();

methods_list = {'method1', 'method2', 'method3', 'method7', 'method8', 'method9'};
for m = 1:length(methods_list)
    z_scores_full_LP.(methods_list{m}) = zeros(H_horizons, bootN);
    z_scores_full_TLP.(methods_list{m}) = zeros(H_horizons, bootN);
    z_scores_full_VAR.(methods_list{m}) = zeros(H_horizons, bootN);
    z_scores_full_SLP.(methods_list{m}) = zeros(H_horizons, bootN);
end

Beta_LP_boot_mean = mean(Beta_LP_boot, 2);
Beta_TLP_boot_mean = mean(Beta_TLP_boot, 2);
Beta_VAR_boot_mean = mean(Beta_VAR_boot, 2);
Beta_SLP_boot_mean = nanmean(Beta_SLP_boot, 2);

for b = 1:bootN
    for h = 1:H_horizons
        % LP
        if ~isnan(Beta_LP_SE_boot(h,b)) && Beta_LP_SE_boot(h,b) > 1e-10
            z_scores_full_LP.method1(h,b) = (Beta_LP_boot(h,b) - Beta_LP_boot_mean(h)) / Beta_LP_SE_boot(h,b);
            z_scores_full_LP.method7(h,b) = z_scores_full_LP.method1(h,b);
            z_scores_full_LP.method2(h,b) = (Beta_LP_boot(h,b) - Beta_VAR_original(h)) / Beta_LP_SE_boot(h,b);
            z_scores_full_LP.method8(h,b) = z_scores_full_LP.method2(h,b);
            z_scores_full_LP.method3(h,b) = (Beta_LP_boot(h,b) - Beta_LP_original(h)) / Beta_LP_SE_boot(h,b);
            z_scores_full_LP.method9(h,b) = z_scores_full_LP.method3(h,b);
        else
            for m = 1:length(methods_list), z_scores_full_LP.(methods_list{m})(h,b) = NaN; end
        end
        
        % TLP
        if ~isnan(Beta_TLP_SE_boot(h,b)) && Beta_TLP_SE_boot(h,b) > 1e-10
            z_scores_full_TLP.method1(h,b) = (Beta_TLP_boot(h,b) - Beta_TLP_boot_mean(h)) / Beta_TLP_SE_boot(h,b);
            z_scores_full_TLP.method7(h,b) = z_scores_full_TLP.method1(h,b);
            z_scores_full_TLP.method2(h,b) = (Beta_TLP_boot(h,b) - Beta_VAR_original(h)) / Beta_TLP_SE_boot(h,b);
            z_scores_full_TLP.method8(h,b) = z_scores_full_TLP.method2(h,b);
            z_scores_full_TLP.method3(h,b) = (Beta_TLP_boot(h,b) - Beta_TLP_original(h)) / Beta_TLP_SE_boot(h,b);
            z_scores_full_TLP.method9(h,b) = z_scores_full_TLP.method3(h,b);
        else
            for m = 1:length(methods_list), z_scores_full_TLP.(methods_list{m})(h,b) = NaN; end
        end
        
        % VAR
        if ~isnan(Beta_VAR_SE_boot(h,b)) && Beta_VAR_SE_boot(h,b) > 1e-10
            z_scores_full_VAR.method1(h,b) = (Beta_VAR_boot(h,b) - Beta_VAR_boot_mean(h)) / Beta_VAR_SE_boot(h,b);
            z_scores_full_VAR.method7(h,b) = z_scores_full_VAR.method1(h,b);
            z_scores_full_VAR.method2(h,b) = (Beta_VAR_boot(h,b) - Beta_VAR_original(h)) / Beta_VAR_SE_boot(h,b);
            z_scores_full_VAR.method8(h,b) = z_scores_full_VAR.method2(h,b);
            z_scores_full_VAR.method3(h,b) = z_scores_full_VAR.method2(h,b);
            z_scores_full_VAR.method9(h,b) = z_scores_full_VAR.method2(h,b);
        else
            for m = 1:length(methods_list), z_scores_full_VAR.(methods_list{m})(h,b) = NaN; end
        end
        
        % SLP
        if compute_SLP && ~isnan(Beta_SLP_SE_boot(h,b)) && Beta_SLP_SE_boot(h,b) > 1e-10
            z_scores_full_SLP.method1(h,b) = (Beta_SLP_boot(h,b) - Beta_SLP_boot_mean(h)) / Beta_SLP_SE_boot(h,b);
            z_scores_full_SLP.method7(h,b) = z_scores_full_SLP.method1(h,b);
            z_scores_full_SLP.method2(h,b) = (Beta_SLP_boot(h,b) - Beta_VAR_original(h)) / Beta_SLP_SE_boot(h,b);
            z_scores_full_SLP.method8(h,b) = z_scores_full_SLP.method2(h,b);
            z_scores_full_SLP.method3(h,b) = (Beta_SLP_boot(h,b) - Beta_SLP_original(h)) / Beta_SLP_SE_boot(h,b);
            z_scores_full_SLP.method9(h,b) = z_scores_full_SLP.method3(h,b);
        else
            for m = 1:length(methods_list), z_scores_full_SLP.(methods_list{m})(h,b) = NaN; end
        end
    end
end

%==========================================================================
% Z-quantiles
%==========================================================================
z_quantiles_LP = struct();
z_quantiles_TLP = struct();
z_quantiles_VAR = struct();
z_quantiles_SLP = struct();

methods_equal_tailed = {'method1', 'method2', 'method3'};
for m = 1:3
    method = methods_equal_tailed{m};
    z_quantiles_LP.(method) = zeros(H_horizons, 2);
    z_quantiles_TLP.(method) = zeros(H_horizons, 2);
    z_quantiles_VAR.(method) = zeros(H_horizons, 2);
    z_quantiles_SLP.(method) = zeros(H_horizons, 2);
    
    for h = 1:H_horizons
        valid_z = z_scores_full_LP.(method)(h, ~isnan(z_scores_full_LP.(method)(h,:)));
        if length(valid_z) > 10
            z_quantiles_LP.(method)(h,:) = [quantile(valid_z, alpha/2), quantile(valid_z, 1-alpha/2)];
        else
            z_quantiles_LP.(method)(h,:) = [norminv(alpha/2), norminv(1-alpha/2)];
        end
        
        valid_z = z_scores_full_TLP.(method)(h, ~isnan(z_scores_full_TLP.(method)(h,:)));
        if length(valid_z) > 10
            z_quantiles_TLP.(method)(h,:) = [quantile(valid_z, alpha/2), quantile(valid_z, 1-alpha/2)];
        else
            z_quantiles_TLP.(method)(h,:) = [norminv(alpha/2), norminv(1-alpha/2)];
        end
        
        valid_z = z_scores_full_VAR.(method)(h, ~isnan(z_scores_full_VAR.(method)(h,:)));
        if length(valid_z) > 10
            z_quantiles_VAR.(method)(h,:) = [quantile(valid_z, alpha/2), quantile(valid_z, 1-alpha/2)];
        else
            z_quantiles_VAR.(method)(h,:) = [norminv(alpha/2), norminv(1-alpha/2)];
        end
        
        valid_z = z_scores_full_SLP.(method)(h, ~isnan(z_scores_full_SLP.(method)(h,:)));
        if length(valid_z) > 10
            z_quantiles_SLP.(method)(h,:) = [quantile(valid_z, alpha/2), quantile(valid_z, 1-alpha/2)];
        else
            z_quantiles_SLP.(method)(h,:) = [norminv(alpha/2), norminv(1-alpha/2)];
        end
    end
end

methods_symmetric = {'method7', 'method8', 'method9'};
source_methods = {'method1', 'method2', 'method3'};
for m = 1:3
    method = methods_symmetric{m};
    source = source_methods{m};
    z_quantiles_LP.(method) = zeros(H_horizons, 2);
    z_quantiles_TLP.(method) = zeros(H_horizons, 2);
    z_quantiles_VAR.(method) = zeros(H_horizons, 2);
    z_quantiles_SLP.(method) = zeros(H_horizons, 2);
    
    for h = 1:H_horizons
        valid_z = z_scores_full_LP.(source)(h, ~isnan(z_scores_full_LP.(source)(h,:)));
        if length(valid_z) > 10
            c = quantile(abs(valid_z), 1-alpha);
        else
            c = norminv(1-alpha/2);
        end
        z_quantiles_LP.(method)(h,:) = [-c, c];
        
        valid_z = z_scores_full_TLP.(source)(h, ~isnan(z_scores_full_TLP.(source)(h,:)));
        if length(valid_z) > 10
            c = quantile(abs(valid_z), 1-alpha);
        else
            c = norminv(1-alpha/2);
        end
        z_quantiles_TLP.(method)(h,:) = [-c, c];
        
        valid_z = z_scores_full_VAR.(source)(h, ~isnan(z_scores_full_VAR.(source)(h,:)));
        if length(valid_z) > 10
            c = quantile(abs(valid_z), 1-alpha);
        else
            c = norminv(1-alpha/2);
        end
        z_quantiles_VAR.(method)(h,:) = [-c, c];
        
        valid_z = z_scores_full_SLP.(source)(h, ~isnan(z_scores_full_SLP.(source)(h,:)));
        if length(valid_z) > 10
            c = quantile(abs(valid_z), 1-alpha);
        else
            c = norminv(1-alpha/2);
        end
        z_quantiles_SLP.(method)(h,:) = [-c, c];
    end
end

%==========================================================================
% Quantile bands
%==========================================================================
Quantile_bands = struct();
Quantile_bands.LP_lower = zeros(H_horizons, 1);
Quantile_bands.LP_upper = zeros(H_horizons, 1);
Quantile_bands.TLP_lower = zeros(H_horizons, 1);
Quantile_bands.TLP_upper = zeros(H_horizons, 1);
Quantile_bands.VAR_lower = zeros(H_horizons, 1);
Quantile_bands.VAR_upper = zeros(H_horizons, 1);
Quantile_bands.SLP_lower = zeros(H_horizons, 1);
Quantile_bands.SLP_upper = zeros(H_horizons, 1);

for h = 1:H_horizons
    valid = Beta_LP_boot(h, ~isnan(Beta_LP_boot(h,:)));
    if length(valid) > 10
        Quantile_bands.LP_lower(h) = quantile(valid, alpha/2);
        Quantile_bands.LP_upper(h) = quantile(valid, 1-alpha/2);
    end
    valid = Beta_TLP_boot(h, ~isnan(Beta_TLP_boot(h,:)));
    if length(valid) > 10
        Quantile_bands.TLP_lower(h) = quantile(valid, alpha/2);
        Quantile_bands.TLP_upper(h) = quantile(valid, 1-alpha/2);
    end
    valid = Beta_VAR_boot(h, ~isnan(Beta_VAR_boot(h,:)));
    if length(valid) > 10
        Quantile_bands.VAR_lower(h) = quantile(valid, alpha/2);
        Quantile_bands.VAR_upper(h) = quantile(valid, 1-alpha/2);
    end
    if compute_SLP
        valid = Beta_SLP_boot(h, ~isnan(Beta_SLP_boot(h,:)));
        if length(valid) > 10
            Quantile_bands.SLP_lower(h) = quantile(valid, alpha/2);
            Quantile_bands.SLP_upper(h) = quantile(valid, 1-alpha/2);
        end
    else
        Quantile_bands.SLP_lower(h) = NaN;
        Quantile_bands.SLP_upper(h) = NaN;
    end
end

%==========================================================================
% Covariance matrices
%==========================================================================
cov_matrix = nan(2,2,H_horizons);
cor_matrix = nan(2,2,H_horizons);
for h = 1:H_horizons
    vector = [Beta_LP_boot(h,:); Beta_VAR_boot(h,:)]';
    valid_rows = ~any(isnan(vector),2);
    if sum(valid_rows) > 10
        cov_matrix(:,:,h) = cov(vector(valid_rows,:));
        cor_matrix(:,:,h) = corr(vector(valid_rows,:));
    end
end

cov_matrix_avg = Variance_matrix_mean;
cor_matrix_avg = nanmean(cor_matrix_boot, 4);

%==========================================================================
% Output
%==========================================================================
output.Beta_LP_original = Beta_LP_original;
output.Beta_TLP_original = Beta_TLP_original;
output.Beta_VAR_original = Beta_VAR_original;
output.Beta_SLP_original = Beta_SLP_original;
output.v_lambda_original = v_lambda_original;
output.lambda_opt = lambda_opt;
output.Variance_matrix_original = Variance_matrix_mean;
output.compute_SLP = compute_SLP;
output.use_bootstrap_bias = use_bootstrap_bias;
output.use_second_level_tlp_variance = use_second_level_tlp_variance;

% Store bootstrap-specific weights (useful for diagnostics)
output.v_lambda_boot = v_lambda_boot;

output.z_scores_full_LP = z_scores_full_LP;
output.z_scores_full_TLP = z_scores_full_TLP;
output.z_scores_full_VAR = z_scores_full_VAR;
output.z_scores_full_SLP = z_scores_full_SLP;

output.z_quantiles_LP = z_quantiles_LP;
output.z_quantiles_TLP = z_quantiles_TLP;
output.z_quantiles_VAR = z_quantiles_VAR;
output.z_quantiles_SLP = z_quantiles_SLP;

output.cov_matrix = cov_matrix;
output.cor_matrix = cor_matrix;
output.cov_matrix_avg = cov_matrix_avg;
output.cor_matrix_avg = cor_matrix_avg;
output.cov_matrix2 = cov_matrix_avg;

output.Quantile_bands = Quantile_bands;

output.Beta_LP_boot = Beta_LP_boot;
output.Beta_TLP_boot = Beta_TLP_boot;
output.Beta_VAR_boot = Beta_VAR_boot;
output.Beta_SLP_boot = Beta_SLP_boot;

output.Beta_LP_SE_boot = Beta_LP_SE_boot;
output.Beta_TLP_SE_boot = Beta_TLP_SE_boot;
output.Beta_TLP_SE_boot_direct = Beta_TLP_SE_boot_direct;
output.Beta_TLP_SE_boot_weighted = Beta_TLP_SE_boot_weighted;
output.Beta_VAR_SE_boot = Beta_VAR_SE_boot;
output.Beta_SLP_SE_boot = Beta_SLP_SE_boot;

end

%==========================================================================
% Helper function
%==========================================================================
function [Y_boot, X_boot] = generate_bootstrap_sample(Y, X, A_hat, U, T, p, k, str, seed, use_iid_bootstrap, block_length_T)
    if use_iid_bootstrap == 1
        indx_boot2 = randi(T, T, 1);
    else
        number_of_blocks = ceil(T/block_length_T);
        temp = randi(number_of_blocks-2, number_of_blocks, 1) * block_length_T + 1;
        temp2 = temp + (0:block_length_T-1);
        indx_boot2 = reshape(temp2', block_length_T*number_of_blocks, 1);
        indx_boot2(T+1:end,:) = [];
    end
    
    U_boot = U(indx_boot2,:);
    U_boot_demean = U_boot - mean(U_boot);
    
    y_t = zeros(T, k);
    y_t(1:p,:) = U_boot_demean(1:p,:);
    for t = (p+1):T
        aaa = zeros(1, k);
        for i = 1:p
            aaa = aaa + (A_hat(k*i-k+1:k*i,:)' * y_t(t-i,:)')';
        end
        y_t(t,:) = aaa + U_boot_demean(t,:);
    end
    
    Y_boot = y_t;
    X_boot = lagmatrix(Y_boot, 1:str.P_VAR);
end
