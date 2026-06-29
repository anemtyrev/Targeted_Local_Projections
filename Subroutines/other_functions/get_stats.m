function stats = get_stats(cell2, irf_true, method_name, str)
    if nargin <= 3
        str.P_VAR = [];
    end
    
    h_max = length(irf_true);
    sim = size(cell2.CI, 3);
    P = size(cell2.CI, 4);
    
    for p = 1:P
        for h = 1:h_max
            % Extract vectors across simulations
            IRF_vec = squeeze(cell2.irf(h, :, p));      % [1 x sim]
            CI_low  = squeeze(cell2.CI(h, 1, :, p));    % [sim x 1]
            CI_high = squeeze(cell2.CI(h, 2, :, p));    % [sim x 1]
            
            % CI statistics
            mean_coverage(1, h)  = mean((irf_true(h) > CI_low) & (irf_true(h) < CI_high));
            mean_length_CI(1, h) = mean(CI_high - CI_low);
            mean_below_CI(1, h)  = mean(irf_true(h) < CI_low);
            mean_above_CI(1, h)  = mean(irf_true(h) > CI_high);
            
            % MSE decomposition
            mean_bias(1, h)      = mean(IRF_vec) - irf_true(h);  % E[θ̂] - θ
            bias_sq(1, h)        = mean_bias(1, h)^2;
            mean_variance(1, h)  = var(IRF_vec, 1);              % population variance (divide by n)
            MSE(1, h)            = mean_variance(1, h) + bias_sq(1, h);
        end
        
        % take averages across variances 

        % === Table: Raw stats ===
        table_out = array2table(real(mean([ ...
            mean_coverage; ...
            mean_below_CI; ...
            mean_above_CI; ...
            mean_length_CI; ...
            bias_sq; ...
            mean_variance; ...
            MSE], 2)));
        table_out.Properties.RowNames = {'Coverage', 'Below_CI', 'Above_CI', ...
            'Length', 'Bias_sq', 'Variance', 'MSE'};
        
        if isempty(str.P_VAR)
            table_out.Properties.VariableNames = method_name;
        else
            table_out.Properties.VariableNames = append(method_name, "_pV=", num2str(str.P_VAR(p)));
        end
        
        % === Table: Root stats ===
        table_out2 = array2table(real(mean([ ...
            mean_coverage; ...
            mean_below_CI; ...
            mean_above_CI; ...
            mean_length_CI; ...
            mean_bias; ...
            sqrt(mean_variance); ...
            sqrt(MSE)], 2)));
        table_out2.Properties.RowNames = {'Coverage', 'Below_CI', 'Above_CI', ...
            'Length', 'Bias', 'Std', 'RMSE'};
        
        if isempty(str.P_VAR)
            table_out2.Properties.VariableNames = method_name;
        else
            table_out2.Properties.VariableNames = append(method_name, "_pV=", num2str(str.P_VAR(p)));
        end
        
        % === Store everything into output struct ===
        stats.coverage(:, :, p)    = mean_coverage;
        stats.below_CI(:, :, p)    = mean_below_CI;
        stats.above_CI(:, :, p)    = mean_above_CI;
        stats.length(:, :, p)      = mean_length_CI;
        stats.bias(:, :, p)        = mean_bias;
        stats.bias_sq(:, :, p)     = bias_sq;
        stats.variance(:, :, p)    = mean_variance;
        stats.MSE(:, :, p)         = MSE;
        stats.bias_abs(:, :, p)    = abs(mean_bias);
        stats.std(:, :, p)         = sqrt(mean_variance);
        stats.RMSE(:, :, p)        = sqrt(MSE);
        stats.table{p}             = round(table_out, 3);
        stats.table2{p}            = round(table_out2, 3);
        stats.method_name          = method_name;
        stats.irf                  = cell2.irf;
        stats.CI                   = cell2.CI;
    end
end