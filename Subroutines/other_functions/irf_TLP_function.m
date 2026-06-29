function output = irf_TLP_function(data, str, use_TLP_boot, use_iid_boot, use_cov_from_var, centering_method, compute_SLP, use_bootstrap_bias, use_second_level_tlp_variance)
%==========================================================================
% irf_TLP_function.m  2025
%
% INPUTS:
%   data              - Struct with fields y_t, etc.
%   str               - Struct with estimation settings
%   use_TLP_boot      - 1 = run bootstrap, 0 = skip
%   use_iid_boot      - 1 = iid bootstrap, 0 = block bootstrap
%   use_cov_from_var  - 1 = use covariance from VAR bootstrap
%   centering_method  - (informational only; all methods computed)
%   compute_SLP       - FLAG: 1 = compute SLP, 0 = skip (default: 0)
%   use_bootstrap_bias - FLAG: 1 = compute TLP weights using bootstrap-layer
%                        bias estimates, 0 = use original sample (default: 0)
%   use_second_level_tlp_variance - FLAG: 1 = compute double-bootstrap TLP
%                        SEs directly from second-level TLP draws with
%                        draw-specific weights, 0 = weighted LP/VAR/COV
%                        variance formula (default: 1)
%==========================================================================

    if nargin < 3, use_TLP_boot = 1; end
    if nargin < 4, use_iid_boot = 0; end
    if nargin < 5, use_cov_from_var = 0; end
    if nargin < 6, centering_method = 1; end
    if nargin < 7, compute_SLP = 0; end
    if nargin < 8, use_bootstrap_bias = 0; end
    if nargin < 9
        if isfield(str, 'use_second_level_tlp_variance')
            use_second_level_tlp_variance = str.use_second_level_tlp_variance;
        else
            use_second_level_tlp_variance = 1;
        end
    end

    H_max = str.H_max;
    H_min = str.H_min;
    H_horizons = H_max - H_min + 1;

    VAR_est = estim_VAR(data, str.P_VAR, 0);
    IRF_est = irf_function(VAR_est.A_hat_f, VAR_est.B_0_inv2, str);
    irf_estim = reshape(permute(IRF_est, [2 1 3]), str.k^2, H_horizons);

    data_reshaped_lp = PrepareLP_VAR(data, H_max, H_min, str.P_LP, 0, str.which_irf_y, str.which_irf_x);

    Beta_OLS = nan(H_horizons, 1);
    Beta_VAR = nan(H_horizons, 1);

    for h = 1:H_horizons
        X = data_reshaped_lp.X{h};
        Y = data_reshaped_lp.Y{h};
        X1 = X(:, 1);
        X_mat = [ones(size(X, 1), 1) X(:, 2:end)];
        Xpr = X1 - X_mat * (X_mat \ X1);
        Ypr = Y - X_mat * (X_mat \ Y);
        Beta_OLS(h) = Xpr \ Ypr;
        Beta_VAR(h) = irf_estim(str.which_irf_var, h);
    end

    % Bootstrap
    if str.use_double_bootstrap == 1
        CI_boot = double_bootstrapTLP_parfor(VAR_est, str.BootN_TLP, str, use_iid_boot, compute_SLP, use_bootstrap_bias, use_second_level_tlp_variance);
        bootstrap_type = 'double';
    else
        CI_boot = bootstrapTLP_parfor(VAR_est, str.BootN_TLP, str, use_iid_boot);
        bootstrap_type = 'single';
    end

    if use_cov_from_var == 1 && isfield(CI_boot, 'cov_matrix2')
        Variance_matrix = CI_boot.cov_matrix2;
    else
        Variance_matrix = CI_boot.cov_matrix_avg;
    end

    % v_lambda and TLP
    v_of_lambda = zeros(H_horizons, 1);
    Beta_TLP = zeros(H_horizons, 1);

    for h = 1:H_horizons
        v_of_lambda(h) = compute_TLP_weight(Beta_OLS(h), Beta_VAR(h), Variance_matrix, h);
        Beta_TLP(h) = v_of_lambda(h)*Beta_OLS(h) + (1-v_of_lambda(h))*Beta_VAR(h);
    end

    % TLP variance
    var_TLP = zeros(H_horizons, 1);
    for h = 1:H_horizons
        var_LP = v_of_lambda(h)^2 * Variance_matrix(1,1,h);
        var_VAR = (1-v_of_lambda(h))^2 * Variance_matrix(2,2,h);
        Cov_LP_VAR = 2*v_of_lambda(h)*(1-v_of_lambda(h))*Variance_matrix(1,2,h);
        var_TLP(h) = var_LP + var_VAR + Cov_LP_VAR;
    end

    var_TLP_weighted = var_TLP;
    if strcmp(bootstrap_type, 'double') && use_second_level_tlp_variance == 1 ...
            && isfield(CI_boot, 'Beta_TLP_SE_boot')
        var_TLP = CI_boot.Beta_TLP_SE_boot(:,1).^2;
    end
    sd_TLP = sqrt(var_TLP);

    % Construct CIs
    if strcmp(bootstrap_type, 'double')
        Beta_LP_SE = CI_boot.Beta_LP_SE_boot(:,1);
        Beta_VAR_SE = CI_boot.Beta_VAR_SE_boot(:,1);
        Beta_TLP_SE = CI_boot.Beta_TLP_SE_boot(:,1);
        
        if compute_SLP
            Beta_SLP = CI_boot.Beta_SLP_original;
            Beta_SLP_SE = CI_boot.Beta_SLP_SE_boot(:,1);
        end

        methods_inverted = {'method1', 'method2', 'method3'};
        for m = 1:3
            mname = methods_inverted{m};
            output.LP.CI.Studentized.(mname)(:,1) = Beta_OLS - CI_boot.z_quantiles_LP.(mname)(:,2).*Beta_LP_SE;
            output.LP.CI.Studentized.(mname)(:,2) = Beta_OLS - CI_boot.z_quantiles_LP.(mname)(:,1).*Beta_LP_SE;
            output.VAR.CI.Studentized.(mname)(:,1) = Beta_VAR - CI_boot.z_quantiles_VAR.(mname)(:,2).*Beta_VAR_SE;
            output.VAR.CI.Studentized.(mname)(:,2) = Beta_VAR - CI_boot.z_quantiles_VAR.(mname)(:,1).*Beta_VAR_SE;
            output.TLP.CI.Studentized.(mname)(:,1) = Beta_TLP - CI_boot.z_quantiles_TLP.(mname)(:,2).*Beta_TLP_SE;
            output.TLP.CI.Studentized.(mname)(:,2) = Beta_TLP - CI_boot.z_quantiles_TLP.(mname)(:,1).*Beta_TLP_SE;
            if compute_SLP
                output.SLP.CI.Studentized.(mname)(:,1) = Beta_SLP - CI_boot.z_quantiles_SLP.(mname)(:,2).*Beta_SLP_SE;
                output.SLP.CI.Studentized.(mname)(:,2) = Beta_SLP - CI_boot.z_quantiles_SLP.(mname)(:,1).*Beta_SLP_SE;
            end
        end

        methods_direct = {'method4', 'method5', 'method6'};
        for m = 1:3
            mname = methods_direct{m};
            centering_m = methods_inverted{m};
            output.LP.CI.Studentized.(mname)(:,1) = Beta_OLS + CI_boot.z_quantiles_LP.(centering_m)(:,1).*Beta_LP_SE;
            output.LP.CI.Studentized.(mname)(:,2) = Beta_OLS + CI_boot.z_quantiles_LP.(centering_m)(:,2).*Beta_LP_SE;
            output.VAR.CI.Studentized.(mname)(:,1) = Beta_VAR + CI_boot.z_quantiles_VAR.(centering_m)(:,1).*Beta_VAR_SE;
            output.VAR.CI.Studentized.(mname)(:,2) = Beta_VAR + CI_boot.z_quantiles_VAR.(centering_m)(:,2).*Beta_VAR_SE;
            output.TLP.CI.Studentized.(mname)(:,1) = Beta_TLP + CI_boot.z_quantiles_TLP.(centering_m)(:,1).*Beta_TLP_SE;
            output.TLP.CI.Studentized.(mname)(:,2) = Beta_TLP + CI_boot.z_quantiles_TLP.(centering_m)(:,2).*Beta_TLP_SE;
            if compute_SLP
                output.SLP.CI.Studentized.(mname)(:,1) = Beta_SLP + CI_boot.z_quantiles_SLP.(centering_m)(:,1).*Beta_SLP_SE;
                output.SLP.CI.Studentized.(mname)(:,2) = Beta_SLP + CI_boot.z_quantiles_SLP.(centering_m)(:,2).*Beta_SLP_SE;
            end
        end

        methods_symmetric = {'method7', 'method8', 'method9'};
        for m = 1:3
            mname = methods_symmetric{m};
            c_LP = CI_boot.z_quantiles_LP.(mname)(:,2);
            c_VAR = CI_boot.z_quantiles_VAR.(mname)(:,2);
            c_TLP = CI_boot.z_quantiles_TLP.(mname)(:,2);
            output.LP.CI.Studentized.(mname)(:,1) = Beta_OLS - c_LP.*Beta_LP_SE;
            output.LP.CI.Studentized.(mname)(:,2) = Beta_OLS + c_LP.*Beta_LP_SE;
            output.VAR.CI.Studentized.(mname)(:,1) = Beta_VAR - c_VAR.*Beta_VAR_SE;
            output.VAR.CI.Studentized.(mname)(:,2) = Beta_VAR + c_VAR.*Beta_VAR_SE;
            output.TLP.CI.Studentized.(mname)(:,1) = Beta_TLP - c_TLP.*Beta_TLP_SE;
            output.TLP.CI.Studentized.(mname)(:,2) = Beta_TLP + c_TLP.*Beta_TLP_SE;
            if compute_SLP
                c_SLP = CI_boot.z_quantiles_SLP.(mname)(:,2);
                output.SLP.CI.Studentized.(mname)(:,1) = Beta_SLP - c_SLP.*Beta_SLP_SE;
                output.SLP.CI.Studentized.(mname)(:,2) = Beta_SLP + c_SLP.*Beta_SLP_SE;
            end
        end

        output.LP.z_scores_full = CI_boot.z_scores_full_LP;
        output.LP.z_quantiles = CI_boot.z_quantiles_LP;
        output.VAR.z_scores_full = CI_boot.z_scores_full_VAR;
        output.VAR.z_quantiles = CI_boot.z_quantiles_VAR;
        output.TLP.z_scores_full = CI_boot.z_scores_full_TLP;
        output.TLP.z_quantiles = CI_boot.z_quantiles_TLP;
        
        if compute_SLP
            output.SLP.z_scores_full = CI_boot.z_scores_full_SLP;
            output.SLP.z_quantiles = CI_boot.z_quantiles_SLP;
            output.SLP.irf = Beta_SLP;
            output.SLP.lambda_opt = CI_boot.lambda_opt;
            output.SLP.variance = var(CI_boot.Beta_SLP_boot, 0, 2);
            output.SLP.std_dev = sqrt(output.SLP.variance);
        end
    else
        % Single bootstrap (no SLP support in single bootstrap)
        Beta_LP_SE_analytical = CI_boot.Beta_LP_SE_analytical;
        methods_inverted = {'method1', 'method2', 'method3'};
        for m = 1:3
            mname = methods_inverted{m};
            output.LP.CI.Studentized.(mname)(:,1) = Beta_OLS - CI_boot.t_quantiles_LP.(mname)(:,2).*Beta_LP_SE_analytical;
            output.LP.CI.Studentized.(mname)(:,2) = Beta_OLS - CI_boot.t_quantiles_LP.(mname)(:,1).*Beta_LP_SE_analytical;
        end
        methods_direct = {'method4', 'method5', 'method6'};
        for m = 1:3
            mname = methods_direct{m};
            centering_m = methods_inverted{m};
            output.LP.CI.Studentized.(mname)(:,1) = Beta_OLS + CI_boot.t_quantiles_LP.(centering_m)(:,1).*Beta_LP_SE_analytical;
            output.LP.CI.Studentized.(mname)(:,2) = Beta_OLS + CI_boot.t_quantiles_LP.(centering_m)(:,2).*Beta_LP_SE_analytical;
        end
        output.LP.z_scores_full = CI_boot.t_scores_LP;
        output.LP.z_quantiles = CI_boot.t_quantiles_LP;
        output.LP.std_dev_analytical = Beta_LP_SE_analytical;
    end

    % Quantile bands
    output.LP.CI.Quantile = [CI_boot.Quantile_bands.LP_lower CI_boot.Quantile_bands.LP_upper];
    output.VAR.CI.Quantile = [CI_boot.Quantile_bands.VAR_lower CI_boot.Quantile_bands.VAR_upper];
    output.TLP.CI.Quantile = [CI_boot.Quantile_bands.TLP_lower CI_boot.Quantile_bands.TLP_upper];
    if compute_SLP && isfield(CI_boot.Quantile_bands, 'SLP_lower')
        output.SLP.CI.Quantile = [CI_boot.Quantile_bands.SLP_lower CI_boot.Quantile_bands.SLP_upper];
    end

    % Output organization
    output.LP.irf = Beta_OLS;
    output.LP.variance = squeeze(CI_boot.cov_matrix_avg(1,1,:));
    output.LP.std_dev = sqrt(output.LP.variance);

    output.VAR.irf = Beta_VAR;
    output.VAR.variance = squeeze(CI_boot.cov_matrix_avg(2,2,:));
    output.VAR.std_dev = sqrt(output.VAR.variance);

    output.TLP.irf = Beta_TLP;
    output.TLP.v_lambda = v_of_lambda;
    output.TLP.variance = var_TLP;
    output.TLP.variance_weighted_formula = var_TLP_weighted;
    output.TLP.std_dev = sd_TLP;

    output.covariance.matrix = Variance_matrix;
    output.covariance.matrix_avg = CI_boot.cov_matrix_avg;
    output.covariance.correlation = CI_boot.cor_matrix;

    output.bootstrap.type = bootstrap_type;
    output.bootstrap.n_boot = str.BootN_TLP;
    output.bootstrap.iid_boot = use_iid_boot;
    output.bootstrap.alpha = str.alpha;
    output.bootstrap.compute_SLP = compute_SLP;
    output.bootstrap.use_bootstrap_bias = use_bootstrap_bias;
    output.bootstrap.use_second_level_tlp_variance = use_second_level_tlp_variance;
    
    % Store bootstrap-specific weights if available
    if isfield(CI_boot, 'v_lambda_boot')
        output.TLP.v_lambda_boot = CI_boot.v_lambda_boot;
    end
end
