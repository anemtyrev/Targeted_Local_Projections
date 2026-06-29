function output = unpack_cell(cell, str)
%==========================================================================
% unpack_cell.m - UPDATED VERSION FOR SINGLE AND DOUBLE BOOTSTRAP
%
% Purpose:
%   Unpacks simulation results from cell array into organized output structure.
%   Supports legacy LP/SLP/VAR/BLP outputs, plus new TLP2 outputs with
%   multiple centering methods for DoubleBoot and combined z-quantiles.
%
% UPDATE: Now properly handles single bootstrap case where DoubleBoot and
%         CombinedQuantile fields don't exist for VAR and TLP.
%         Only creates output fields when data actually exists.
%
% UPDATE 2: Now directly accesses 'variance' field from bootstrap output
%           instead of computing variance from beta_sort.
%
% UPDATE 3: Added support for symmetric bootstrap methods 7-9 (Hall 1988)
%
% UPDATE 4: Added SLP2 outputs with all 9 centering methods
%
% Inputs:
%   cell - cell array of simulation results
%   str  - structure containing model specifications
%
% Outputs:
%   output - organized structure with all estimation results
%     Legacy fields:
%       .LP, .SLP, .VAR, .BLP
%     New fields:
%       .LP2, .VAR2, .TLP2, .SLP2
%==========================================================================
    simulations = length(cell);
    output = struct();
    
    % Include symmetric methods 7-9
    methods = {'method1','method2','method3','method4','method5','method6',...
               'method7','method8','method9'};

    for i = 1:simulations

        %% ================= Legacy LP =================
        if isfield(cell{i}, 'irf_LP') && ~isempty(cell{i}.irf_LP)
            if isfield(cell{i}.irf_LP, 'irf')
                output.LP.B.irf(:,i) = safe_get(cell{i}, 'irf_LP', 'irf', []);
                output.LP.H.irf(:,i) = safe_get(cell{i}, 'irf_LP', 'irf', []);
            end
            if isfield(cell{i}.irf_LP, 'CI_Boot')
                output.LP.B.CI(:,:,i) = safe_get(cell{i}, 'irf_LP', 'CI_Boot', 'conf_int', []);
                output.LP.B.variance(:,:,i) = safe_var(cell{i}, 'irf_LP', 'CI_Boot');
            end
            if isfield(cell{i}.irf_LP, 'CI_HAC')
                output.LP.H.CI(:,:,i) = safe_get(cell{i}, 'irf_LP', 'CI_HAC', 'conf_int', []);
                output.LP.H.variance(:,:,i) = safe_diag(cell{i}, 'irf_LP', 'CI_HAC', 'VAR');
            end
        end

        %% ================= Legacy SLP =================
        if isfield(cell{i}, 'irf_SLP') && ~isempty(cell{i}.irf_SLP)
            if isfield(cell{i}.irf_SLP, 'irf')
                output.SLP.B.irf(:,i) = safe_get(cell{i}, 'irf_SLP', 'irf', []);
                output.SLP.H.irf(:,i) = safe_get(cell{i}, 'irf_SLP', 'irf', []);
                output.SLP.H_US.irf(:,i) = safe_get(cell{i}, 'irf_SLP', 'irf', []);
            end
            if isfield(cell{i}.irf_SLP, 'CI_Boot')
                output.SLP.B.CI(:,:,i) = safe_get(cell{i}, 'irf_SLP', 'CI_Boot', 'conf_int', []);
                output.SLP.B.variance(:,:,i) = safe_var(cell{i}, 'irf_SLP', 'CI_Boot');
            end
            if isfield(cell{i}.irf_SLP, 'CI_HAC')
                output.SLP.H.CI(:,:,i) = safe_get(cell{i}, 'irf_SLP', 'CI_HAC', 'conf_int', []);
                output.SLP.H.variance(:,:,i) = safe_diag(cell{i}, 'irf_SLP', 'CI_HAC', 'VAR');
            end
            if isfield(cell{i}.irf_SLP, 'CI_HAC_UNDERSMOOTH')
                output.SLP.H_US.CI(:,:,i) = safe_get(cell{i}, 'irf_SLP', 'CI_HAC_UNDERSMOOTH', 'conf_int', []);
                output.SLP.H_US.variance(:,:,i) = safe_diag(cell{i}, 'irf_SLP', 'CI_HAC_UNDERSMOOTH', 'VAR');
            end
        end

        %% ================= Legacy VAR =================
        if isfield(cell{i}, 'irf_VAR') && ~isempty(cell{i}.irf_VAR)
            for p = 1:length(cell{i}.irf_VAR)
                if isfield(cell{i}.irf_VAR{p}, 'irf')
                    output.VAR.B.irf(:, i, p) = safe_get(cell{i}, 'irf_VAR', p, 'irf', []);
                end
                if isfield(cell{i}.irf_VAR{p}, 'CI')
                    output.VAR.B.CI(:, :, i, p) = safe_get(cell{i}, 'irf_VAR', p, 'CI', 'conf_int', []);
                    output.VAR.B.variance(:, :, i, p) = safe_stddev(cell{i}, 'irf_VAR', p, 'CI', 'standard_dev_est', str);
                end
            end
        end

        %% ================= Legacy BLP =================
        if isfield(cell{i}, 'irf_BLP') && ~isempty(cell{i}.irf_BLP)
            for p = 1:length(cell{i}.irf_BLP)
                output.BLP.irf(:, i, p) = safe_get(cell{i}, 'irf_BLP', p, 'irf', []);
                output.BLP.CI(:, :, i, p) = safe_get(cell{i}, 'irf_BLP', p, 'CI', []);
                output.BLP.variance(:, :, i, p) = safe_get(cell{i}, 'irf_BLP', p, 'var', []);
            end
        end

        %% ================= New TLP2-style outputs =================
        if isfield(cell{i}, 'irf_TLP2') && ~isempty(cell{i}.irf_TLP2)
            for p = 1:length(str.P_VAR)

                % ---------------- LP2 ----------------
                % Studentized (check if field exists and has data)
                for m = 1:length(methods)
                    method_name = methods{m};
                    studentized_ci = safe_get(cell{i}, 'irf_TLP2', p, 'LP', 'CI', 'Studentized', method_name, []);
                    if ~isempty(studentized_ci)
                        output.LP2.Studentized.(method_name).irf(:,i,p) = safe_get(cell{i}, 'irf_TLP2', p, 'LP', 'irf', []);
                        output.LP2.Studentized.(method_name).CI(:, :, i, p) = studentized_ci;
                        output.LP2.Studentized.(method_name).variance(:, :, i, p) = safe_get(cell{i}, 'irf_TLP2', p, 'LP', 'variance', [])';
                    end
                end
                % Quantile (always present)
                quantile_ci = safe_get(cell{i}, 'irf_TLP2', p, 'LP', 'CI', 'Quantile', []);
                if ~isempty(quantile_ci)
                    output.LP2.Quantile.irf(:,i,p) = safe_get(cell{i}, 'irf_TLP2', p, 'LP', 'irf', []);
                    output.LP2.Quantile.CI(:, :, i, p) = quantile_ci;
                    output.LP2.Quantile.variance(:, :, i, p) = safe_get(cell{i}, 'irf_TLP2', p, 'LP', 'variance', [])';
                end

                % ---------------- VAR2 ----------------
                % Studentized (check if field exists and has data)
                for m = 1:length(methods)
                    method_name = methods{m};
                    studentized_ci = safe_get(cell{i}, 'irf_TLP2', p, 'VAR', 'CI', 'Studentized', method_name, []);
                    if ~isempty(studentized_ci)
                        output.VAR2.Studentized.(method_name).irf(:,i,p) = safe_get(cell{i}, 'irf_TLP2', p, 'VAR', 'irf', []);
                        output.VAR2.Studentized.(method_name).CI(:, :, i, p) = studentized_ci;
                        output.VAR2.Studentized.(method_name).variance(:, :, i, p) = safe_get(cell{i}, 'irf_TLP2', p, 'VAR', 'variance', [])';
                    end
                end
                % Quantile (always present)
                quantile_ci = safe_get(cell{i}, 'irf_TLP2', p, 'VAR', 'CI', 'Quantile', []);
                if ~isempty(quantile_ci)
                    output.VAR2.Quantile.irf(:,i,p) = safe_get(cell{i}, 'irf_TLP2', p, 'VAR', 'irf', []);
                    output.VAR2.Quantile.CI(:, :, i, p) = quantile_ci;
                    output.VAR2.Quantile.variance(:, :, i, p) = safe_get(cell{i}, 'irf_TLP2', p, 'VAR', 'variance', [])';
                end

                % ---------------- TLP2 ----------------
                % Studentized (check if field exists and has data)
                for m = 1:length(methods)
                    method_name = methods{m};
                    
                    studentized_ci = safe_get(cell{i}, 'irf_TLP2', p, 'TLP', 'CI', 'Studentized', method_name, []);
                    if ~isempty(studentized_ci)
                        output.TLP2.Studentized.(method_name).irf(:,i,p) = safe_get(cell{i}, 'irf_TLP2', p, 'TLP', 'irf', []);
                        output.TLP2.Studentized.(method_name).CI(:, :, i, p) = studentized_ci;
                        output.TLP2.Studentized.(method_name).variance(:, :, i, p) = safe_get(cell{i}, 'irf_TLP2', p, 'TLP', 'variance', [])';
                        variance_weighted = safe_get(cell{i}, 'irf_TLP2', p, 'TLP', 'variance_weighted_formula', []);
                        if ~isempty(variance_weighted)
                            output.TLP2.Studentized.(method_name).variance_weighted(:, :, i, p) = variance_weighted';
                        end
                        output.TLP2.Studentized.(method_name).v_lambda(:,i,p) = safe_get(cell{i}, 'irf_TLP2', p, 'TLP', 'v_lambda', []);
                    end
                    
                    % CombinedQuantile CI (check if exists and has data)
                    combined_ci = safe_get(cell{i}, 'irf_TLP2', p, 'TLP', 'CI', 'CombinedQuantile', method_name, []);
                    if ~isempty(combined_ci)
                        output.TLP2.CombinedQuantile.(method_name).CI(:, :, i, p) = combined_ci;
                        output.TLP2.CombinedQuantile.(method_name).irf(:, i, p) = safe_get(cell{i}, 'irf_TLP2', p, 'TLP', 'irf', []);
                        output.TLP2.CombinedQuantile.(method_name).variance(:, :, i, p) = safe_get(cell{i}, 'irf_TLP2', p, 'TLP', 'variance', [])';
                        output.TLP2.CombinedQuantile.(method_name).v_lambda(:, i, p) = safe_get(cell{i}, 'irf_TLP2', p, 'TLP', 'v_lambda', []);
                    end
                end
                
                % Quantile bands (always present)
                quantile_ci = safe_get(cell{i}, 'irf_TLP2', p, 'TLP', 'CI', 'Quantile', []);
                if ~isempty(quantile_ci)
                    output.TLP2.Quantile.irf(:,i,p) = safe_get(cell{i}, 'irf_TLP2', p, 'TLP', 'irf', []);
                    output.TLP2.Quantile.CI(:, :, i, p) = quantile_ci;
                    output.TLP2.Quantile.variance(:, :, i, p) = safe_get(cell{i}, 'irf_TLP2', p, 'TLP', 'variance', [])';
                    output.TLP2.Quantile.v_lambda(:,i,p) = safe_get(cell{i}, 'irf_TLP2', p, 'TLP', 'v_lambda', []);
                end

                % ---------------- SLP2 (NEW) ----------------
                % Studentized (check if field exists and has data)
                for m = 1:length(methods)
                    method_name = methods{m};
                    studentized_ci = safe_get(cell{i}, 'irf_TLP2', p, 'SLP', 'CI', 'Studentized', method_name, []);
                    if ~isempty(studentized_ci)
                        output.SLP2.Studentized.(method_name).irf(:,i,p) = safe_get(cell{i}, 'irf_TLP2', p, 'SLP', 'irf', []);
                        output.SLP2.Studentized.(method_name).CI(:, :, i, p) = studentized_ci;
                        output.SLP2.Studentized.(method_name).variance(:, :, i, p) = safe_get(cell{i}, 'irf_TLP2', p, 'SLP', 'variance', [])';
                        output.SLP2.Studentized.(method_name).lambda_opt(:,i,p) = safe_get(cell{i}, 'irf_TLP2', p, 'SLP', 'lambda_opt', []);
                    end
                end
                
                % Quantile bands for SLP
                quantile_ci = safe_get(cell{i}, 'irf_TLP2', p, 'SLP', 'CI', 'Quantile', []);
                if ~isempty(quantile_ci)
                    output.SLP2.Quantile.irf(:,i,p) = safe_get(cell{i}, 'irf_TLP2', p, 'SLP', 'irf', []);
                    output.SLP2.Quantile.CI(:, :, i, p) = quantile_ci;
                    output.SLP2.Quantile.variance(:, :, i, p) = safe_get(cell{i}, 'irf_TLP2', p, 'SLP', 'variance', [])';
                    output.SLP2.Quantile.lambda_opt(:,i,p) = safe_get(cell{i}, 'irf_TLP2', p, 'SLP', 'lambda_opt', []);
                end

            end
        end

    end
end

%% ================= Helper Functions =================
function val = safe_get(s, varargin)
    % Safe recursive field/cell extraction with default fallback
    default = varargin{end};
    fields = varargin(1:end-1);
    for i = 1:length(fields)
        f = fields{i};
        if ischar(f)
            if isfield(s, f)
                s = s.(f);
            else
                val = default; 
                return;
            end
        elseif isnumeric(f) && iscell(s)
            if f > 0 && f <= length(s)
                s = s{f};
            else
                val = default; 
                return;
            end
        else
            val = default; 
            return;
        end
    end
    val = s;
end

function val = safe_var(s, field1, field2)
    % Updated: Now directly accesses 'variance' field instead of computing from beta_sort
    if isfield(s, field1) && isfield(s.(field1), field2)
        if isfield(s.(field1).(field2), 'variance')
            val = s.(field1).(field2).variance;
        elseif isfield(s.(field1).(field2), 'beta_sort')
            % Fallback for legacy data that still has beta_sort
            val = var(s.(field1).(field2).beta_sort);
        else
            val = [];
        end
    else
        val = [];
    end
end

function val = safe_diag(s, field1, field2, field3)
    if isfield(s, field1) && isfield(s.(field1), field2) && isfield(s.(field1).(field2), field3)
        val = diag(s.(field1).(field2).(field3))';
    else
        val = [];
    end
end

function val = safe_stddev(s, field1, index, field2, field3, str)
    if isfield(str, 'which_irf_var') && isfield(s, field1)
        sub = s.(field1);
        if iscell(sub) && index <= length(sub) && isstruct(sub{index}) && ...
           isfield(sub{index}, field2) && isfield(sub{index}.(field2), field3)
            data = sub{index}.(field2).(field3);
            if str.which_irf_var > 0 && str.which_irf_var <= size(data, 1)
                val = data(str.which_irf_var, :).^2;
                return;
            end
        end
    end
    val = [];
end
