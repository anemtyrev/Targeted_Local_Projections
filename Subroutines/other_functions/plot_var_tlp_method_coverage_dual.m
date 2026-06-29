function plot_var_tlp_method_coverage_dual(mat_files, p_value, varargin)
% PLOT_VAR_TLP_METHOD_COVERAGE_DUAL  4x2 method-by-T metric comparison
%
% Rows are the four studentized CI variants:
%   method2: Subtract VAR
%   method1: Subtract Bootstrap Mean
%   method8: Symmetric around VAR
%   method7: Subtract Bootstrap Mean + Symmetric
%
% Columns are the two sample sizes. Each panel overlays selected estimators
% for a selected metric and a single VAR lag order, by default pV = 8.
%
% Example:
%   plot_var_tlp_method_coverage_dual( ...
%       {'DiffEta_T200_sim1000_pL10_pV8_eta1_DGP_OleaSW_P1.mat', ...
%        'DiffEta_T800_sim1000_pL10_pV8_eta1_DGP_OleaSW_P1.mat'}, ...
%       8, 'coverage', 1, false);
%
% With extra estimators:
%   plot_var_tlp_method_coverage_dual(mat_files, 8, 'coverage', ...
%       {'LP', 'SLP', 'VAR', 'TLP'}, 1, false);
%
% Backward-compatible old call:
%   plot_var_tlp_method_coverage_dual(mat_files, 8, 1, false);

    if nargin < 2 || isempty(p_value), p_value = 8; end
    [metric, estimators, save_picture, close_picture] = parse_plot_args(varargin{:});

    if ~iscell(mat_files) || numel(mat_files) ~= 2
        error('mat_files must be a 1x2 cell array: {T200_file, T800_file}.');
    end

    method_fields = {'method2', 'method1', 'method8', 'method7'};
    row_labels = { ...
        'Subtract VAR', ...
        'Subtract Bootstrap Mean', ...
        'Symmetric around VAR', ...
        'Subtract Boot. Mean + Symmetric'};

    data = cell(1, 2);
    for col = 1:2
        data{col} = load_plot_data(mat_files{col}, p_value);
    end

    metric_lower = lower(char(metric));
    [y_label, y_lim, y_ticks, ref_line, ref_label, ref_line_style] = ...
        get_metric_settings(metric_lower);

    if isempty(y_lim)
        y_lim = compute_auto_y_limits(data, method_fields, metric_lower, estimators);
    end

    [estimator_colors, estimator_line_styles, estimator_line_widths] = get_estimator_styles(estimators);
    ref_color = [0.5, 0.5, 0.5];

    line_width_ref = 3;

    axis_font = 18;
    label_font = 20;
    title_font = 22;
    legend_font = 14;
    row_label_font = 16;

    figure;
    fig_width = 0.37;
    fig_height = min(0.95, 0.14 * 4 + 0.40);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05 0.025 fig_width fig_height]);

    tiledlayout(4, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    legend_built = false;
    legend_handles = [];
    legend_entries = {};

    for row = 1:4
        for col = 1:2
            nexttile;
            hold on;

            h_ref = [];
            if ~isempty(ref_line)
                h_ref = yline(ref_line, 'LineWidth', line_width_ref, ...
                    'Color', ref_color, 'LineStyle', ref_line_style);
                uistack(h_ref, 'bottom');
            end

            max_horizon = 0;
            estimator_handles = [];
            for est_idx = 1:length(estimators)
                estimator = estimators{est_idx};
                method_data = resolve_field_path(data{col}.stats, ...
                    [estimator '.Studentized.' method_fields{row}]);
                [horizons, values] = extract_metric_values(method_data, metric_lower, data{col}.p_idx);
                max_horizon = max(max_horizon, max(horizons));

                h_est = plot(horizons, values, ...
                    'Color', estimator_colors(est_idx, :), ...
                    'LineWidth', estimator_line_widths(est_idx), ...
                    'LineStyle', estimator_line_styles{est_idx});
                estimator_handles = [estimator_handles, h_est]; %#ok<AGROW>
            end

            hold off;

            xlim([0 max_horizon]);
            if ~isempty(y_lim)
                ylim(y_lim);
            end
            xticks(unique([0:5:max_horizon, max_horizon]));
            if ~isempty(y_ticks)
                yticks(y_ticks);
            end

            ax = gca;
            ax.YAxis.TickLabelFormat = '%.1f';
            ax.Layer = 'bottom';
            set(gca, 'FontSize', axis_font);
            grid on;
            box on;
            if ~isempty(h_ref)
                uistack(h_ref, 'bottom');
            end
            if ~isempty(estimator_handles)
                uistack(estimator_handles, 'top');
            end

            current_ylim = ylim;
            y_pos = current_ylim(1) + 0.03 * (current_ylim(2) - current_ylim(1));
            text(1, y_pos, row_labels{row}, ...
                'FontSize', row_label_font, ...
                'FontName', get(gca, 'FontName'), ...
                'Color', [0.2 0.2 0.2], ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'bottom', ...
                'BackgroundColor', [1 1 1 0.8], ...
                'EdgeColor', [0.6 0.6 0.6], ...
                'Margin', 3);

            if row == 1
                title(sprintf('T = %d', round(data{col}.T)), ...
                    'FontSize', title_font, 'FontWeight', 'bold');
            end
            if col == 1
                ylabel(y_label, 'FontSize', label_font);
            end
            if row == 4
                xlabel('Horizon', 'FontSize', label_font);
            end

            if ~legend_built
                legend_handles = estimator_handles;
                legend_entries = estimators;
                if ~isempty(ref_label) && ~isempty(h_ref)
                    legend_handles = [legend_handles, h_ref];
                    legend_entries = [legend_entries, {ref_label}];
                end
                legend_built = true;
            end
        end
    end

    if ~isempty(legend_handles)
        lgd = legend(legend_handles, legend_entries, ...
            'Orientation', 'horizontal', ...
            'FontSize', legend_font, ...
            'NumColumns', 3);
        lgd.Layout.Tile = 'north';
        lgd.ItemTokenSize = [40, 14];
    end

    set(gcf, 'Color', [1 1 1]);

    align_columns_like_vertical_dual_plot();

    if should_save_figure(save_picture)
        save_var_tlp_method_coverage_figure(data, p_value, save_picture, metric_lower, estimators);
    end

    if close_picture
        close(gcf);
    end
end

function [metric, estimators, save_picture, close_picture] = parse_plot_args(varargin)
    metric = 'coverage';
    estimators = {'VAR', 'TLP'};
    save_picture = 0;
    close_picture = false;

    if isempty(varargin)
        return;
    end

    first_arg = varargin{1};
    if ischar(first_arg) || (isstring(first_arg) && isscalar(first_arg))
        metric = char(first_arg);
        remaining_args = varargin(2:end);
        [estimators, save_picture, close_picture] = parse_remaining_plot_args(remaining_args, estimators);
        return;
    end

    if is_estimator_spec(first_arg)
        estimators = normalize_estimators(first_arg);
        remaining_args = varargin(2:end);
        [estimators, save_picture, close_picture] = parse_remaining_plot_args(remaining_args, estimators);
        return;
    end

    save_picture = first_arg;
    if length(varargin) >= 2
        close_picture = varargin{2};
    end
    if length(varargin) >= 3 && (ischar(varargin{3}) || (isstring(varargin{3}) && isscalar(varargin{3})))
        metric = char(varargin{3});
    end
    if length(varargin) >= 4 && is_estimator_spec(varargin{4})
        estimators = normalize_estimators(varargin{4});
    end
end

function [estimators, save_picture, close_picture] = parse_remaining_plot_args(args, default_estimators)
    estimators = default_estimators;
    save_picture = 0;
    close_picture = false;

    if isempty(args)
        return;
    end

    if is_estimator_spec(args{1})
        estimators = normalize_estimators(args{1});
        args = args(2:end);
    end

    if ~isempty(args)
        save_picture = args{1};
    end
    if length(args) >= 2
        close_picture = args{2};
    end
end

function tf = is_estimator_spec(value)
    tf = false;
    if iscell(value)
        tf = all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), value));
    elseif isstring(value)
        tf = numel(value) >= 1 && all(ismember(upper(cellstr(value)), {'LP', 'SLP', 'VAR', 'TLP', 'BLP'}));
    elseif ischar(value)
        tf = ismember(upper(value), {'LP', 'SLP', 'VAR', 'TLP', 'BLP'});
    end
end

function estimators = normalize_estimators(value)
    if ischar(value) || (isstring(value) && isscalar(value))
        estimators = {upper(char(value))};
    elseif isstring(value)
        estimators = cellstr(upper(value));
    elseif iscell(value)
        estimators = cellfun(@(x) upper(char(x)), value, 'UniformOutput', false);
    else
        error('Estimator specification must be a string or cell array of strings.');
    end

    allowed_estimators = {'LP', 'SLP', 'VAR', 'TLP', 'BLP'};
    for i = 1:length(estimators)
        if ~ismember(estimators{i}, allowed_estimators)
            error('Unsupported estimator "%s". Use LP, SLP, VAR, TLP, or BLP.', estimators{i});
        end
    end
end

function [y_label, y_lim, y_ticks, ref_line, ref_label, ref_line_style] = get_metric_settings(metric)
    y_lim = [];
    y_ticks = [];
    ref_line = [];
    ref_label = '';
    ref_line_style = '-';

    switch metric
        case 'coverage'
            y_label = 'Coverage';
            y_lim = [0.5 1];
            y_ticks = 0.5:0.1:1.0;
            ref_line = 0.9;
            ref_label = '90 %';
        case 'length'
            y_label = 'Length';
        case 'bias'
            y_label = 'Bias';
            ref_line = 0;
            ref_label = 'Zero';
        case 'rmse'
            y_label = 'RMSE';
        case 'mse'
            y_label = 'MSE';
        case 'variance'
            y_label = 'Variance';
        case 'bias_sq'
            y_label = 'Bias Squared';
        case 'std'
            y_label = 'Std';
        case 'below_ci'
            y_label = 'Below CI';
            ref_line = 0.05;
            ref_label = '5 %';
            ref_line_style = '--';
        case 'above_ci'
            y_label = 'Above CI';
            ref_line = 0.05;
            ref_label = '5 %';
            ref_line_style = '--';
        otherwise
            y_label = metric;
    end
end

function [estimator_colors, estimator_line_styles, estimator_line_widths] = get_estimator_styles(estimators)
    base_colors = containers.Map();
    base_colors('TLP') = [0.9, 0.5, 0.1];
    base_colors('LP')  = [0.2, 0.6, 0.4];
    base_colors('VAR') = [0.8, 0.4, 0.8];
    base_colors('SLP') = [0.1, 0.6, 0.8];
    base_colors('BLP') = [0.2, 0.4, 0.6];

    style_map = containers.Map();
    style_map('TLP') = '-';
    style_map('LP')  = ':';
    style_map('VAR') = '--';
    style_map('SLP') = '-.';
    style_map('BLP') = '-';

    estimator_colors = zeros(length(estimators), 3);
    estimator_line_styles = cell(1, length(estimators));
    estimator_line_widths = 3 * ones(1, length(estimators));

    for i = 1:length(estimators)
        estimator = estimators{i};
        if ~isKey(base_colors, estimator)
            error('Unsupported estimator "%s".', estimator);
        end
        estimator_colors(i, :) = base_colors(estimator);
        estimator_line_styles{i} = style_map(estimator);
        if strcmp(estimator, 'TLP')
            estimator_line_widths(i) = 5;
        end
    end
end

function y_lim = compute_auto_y_limits(data, method_fields, metric, estimators)
    all_values = [];
    for col = 1:2
        for row = 1:length(method_fields)
            for est_idx = 1:length(estimators)
                method_data = resolve_field_path(data{col}.stats, ...
                    [estimators{est_idx} '.Studentized.' method_fields{row}]);
                [~, values] = extract_metric_values(method_data, metric, data{col}.p_idx);
                all_values = [all_values, values]; %#ok<AGROW>
            end
        end
    end

    if isempty(all_values)
        y_lim = [];
        return;
    end

    y_min = min(all_values);
    y_max = max(all_values);
    y_range = y_max - y_min;
    if y_range <= 0
        y_range = max(abs(y_max), 1);
    end

    y_lim = [y_min - 0.05*y_range, y_max + 0.05*y_range];

    if ismember(metric, {'length', 'rmse', 'mse', 'variance', 'bias_sq', 'std'})
        y_lim(1) = max(0, y_lim(1));
    elseif strcmp(metric, 'bias')
        y_lim(1) = min(y_lim(1), -0.05*y_range);
        y_lim(2) = max(y_lim(2), 0.05*y_range);
    end
end

function align_columns_like_vertical_dual_plot()
    all_axes = flipud(findobj(gcf, 'Type', 'axes'));
    real_axes = [];
    for i = 1:length(all_axes)
        if isempty(all_axes(i).Tag) || ~strcmp(all_axes(i).Tag, 'legend')
            real_axes = [real_axes; all_axes(i)]; %#ok<AGROW>
        end
    end

    left_axes = [];
    right_axes = [];
    for i = 1:length(real_axes)
        pos = real_axes(i).Position;
        if pos(1) < 0.5
            left_axes = [left_axes; real_axes(i)]; %#ok<AGROW>
        else
            right_axes = [right_axes; real_axes(i)]; %#ok<AGROW>
        end
    end

    if ~isempty(left_axes)
        max_left = max(arrayfun(@(ax) ax.Position(1), left_axes));
        min_width_L = min(arrayfun(@(ax) ax.Position(3), left_axes));
        for i = 1:length(left_axes)
            pos = left_axes(i).Position;
            left_axes(i).Position = [max_left, pos(2), min_width_L, pos(4)];
        end
    end

    if ~isempty(right_axes)
        max_right = max(arrayfun(@(ax) ax.Position(1), right_axes));
        min_width_R = min(arrayfun(@(ax) ax.Position(3), right_axes));
        for i = 1:length(right_axes)
            pos = right_axes(i).Position;
            right_axes(i).Position = [max_right, pos(2), min_width_R, pos(4)];
        end
    end
end

function data = load_plot_data(mat_file, p_value)
    resolved_file = resolve_mat_file(mat_file);
    fprintf('Loading: %s\n', resolved_file);

    loaded = load(resolved_file, 'stats', 'str');
    if ~isfield(loaded, 'stats') || ~isfield(loaded, 'str')
        error('File "%s" must contain variables "stats" and "str".', resolved_file);
    end

    p_values = get_numeric_vector_field(loaded.str, 'P_VAR');
    if isempty(p_values)
        error('File "%s" does not contain a valid str.P_VAR.', resolved_file);
    end

    p_idx = find(abs(p_values - p_value) < 1e-10, 1);
    if isempty(p_idx)
        error('Requested pV=%g, but file "%s" has P_VAR=%s.', ...
            p_value, resolved_file, mat2str(p_values));
    end

    T_value = get_scalar_numeric_field(loaded.str, 'T');
    if isnan(T_value)
        T_value = extract_T_from_filename(resolved_file);
    end

    data = struct();
    data.stats = loaded.stats;
    data.str = loaded.str;
    data.file = resolved_file;
    data.p_values = p_values;
    data.p_idx = p_idx;
    data.T = T_value;
end

function resolved_file = resolve_mat_file(mat_file)
    if isstring(mat_file) && isscalar(mat_file)
        mat_file = char(mat_file);
    end
    if ~ischar(mat_file)
        error('Each mat_files entry must be a file path string.');
    end

    if exist(mat_file, 'file')
        resolved_file = mat_file;
        return;
    end

    function_dir = fileparts(mfilename('fullpath'));
    candidate = fullfile(function_dir, mat_file);
    if exist(candidate, 'file')
        resolved_file = candidate;
        return;
    end

    error('Could not find MAT file "%s".', mat_file);
end

function result = resolve_field_path(s, path)
    result = s;
    parts = strsplit(path, '.');
    for i = 1:length(parts)
        if isstruct(result) && isfield(result, parts{i})
            result = result.(parts{i});
        else
            error('Field path "%s" not found; failed at "%s".', path, parts{i});
        end
    end
end

function [horizons, values] = extract_metric_values(method, metric, p_idx)
    field_name = resolve_metric_field(method, metric);
    if isempty(field_name)
        error('The selected method does not contain metric "%s".', metric);
    end

    metric_data = method.(field_name);

    if ndims(metric_data) == 3
        if p_idx > size(metric_data, 3)
            error('Requested p index %d exceeds %s third dimension %d.', ...
                p_idx, field_name, size(metric_data, 3));
        end
        values = squeeze(metric_data(:, :, p_idx));
        values = values(:)';
    elseif isvector(metric_data)
        values = metric_data(:)';
    elseif ismatrix(metric_data)
        [nrows, ncols] = size(metric_data);
        if nrows == 1
            values = metric_data(1, :);
        elseif ncols == 1
            values = metric_data(:, 1)';
        elseif p_idx <= ncols
            values = metric_data(:, p_idx)';
        else
            error('Could not extract p index %d from %s matrix of size %dx%d.', ...
                p_idx, field_name, nrows, ncols);
        end
    else
        error('Unsupported %s array shape.', field_name);
    end

    values = real(values);
    horizons = 0:(length(values) - 1);
end

function field_name = resolve_metric_field(method, metric)
    if isempty(metric)
        field_name = '';
        return;
    end

    switch lower(metric)
        case 'rmse'
            field_variants = {'RMSE', 'rmse', 'Rmse'};
        case 'mse'
            field_variants = {'MSE', 'mse', 'Mse'};
        case 'bias_sq'
            field_variants = {'bias_sq', 'Bias_sq', 'BIAS_SQ', 'bias_squared'};
        case 'bias_abs'
            field_variants = {'bias_abs', 'Bias_abs', 'BIAS_ABS'};
        case 'below_ci'
            field_variants = {'below_CI', 'below_ci', 'belowCI', 'BELOW_CI'};
        case 'above_ci'
            field_variants = {'above_CI', 'above_ci', 'aboveCI', 'ABOVE_CI'};
        otherwise
            metric = char(metric);
            field_variants = {metric, lower(metric), upper(metric), ...
                [upper(metric(1)) lower(metric(2:end))]};
    end

    field_name = '';
    for i = 1:length(field_variants)
        if isfield(method, field_variants{i})
            field_name = field_variants{i};
            return;
        end
    end
end

function values = get_numeric_vector_field(s, field_name)
    values = [];
    if isstruct(s) && isfield(s, field_name) && isnumeric(s.(field_name))
        raw = s.(field_name);
        values = unique(raw(:)', 'sorted');
    end
end

function value = get_scalar_numeric_field(s, field_name)
    value = NaN;
    if isstruct(s) && isfield(s, field_name) && isnumeric(s.(field_name)) && ~isempty(s.(field_name))
        value = s.(field_name)(1);
    end
end

function T_value = extract_T_from_filename(file_path)
    T_value = NaN;
    [~, base_name, ~] = fileparts(file_path);
    tokens = regexp(base_name, '(?:^|[_-])T(\d+)(?:[_-]|$)', 'tokens', 'once');
    if ~isempty(tokens)
        T_value = str2double(tokens{1});
    end
end

function value = should_save_figure(save_picture)
    value = false;
    if ischar(save_picture) || (isstring(save_picture) && isscalar(save_picture))
        value = true;
    elseif (isnumeric(save_picture) || islogical(save_picture)) && any(save_picture)
        value = true;
    end
end

function save_var_tlp_method_coverage_figure(data, p_value, save_picture, metric, estimators)
    output_dir = 'TablesAndPlots';
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    T1 = round(data{1}.T);
    T2 = round(data{2}.T);
    estimator_tag = strjoin(estimators, '_');
    base_name = sprintf('%s_METHOD_%s_t%d_t%d_pV%g', ...
        estimator_tag, upper(metric), T1, T2, p_value);

    if ischar(save_picture) || (isstring(save_picture) && isscalar(save_picture))
        base_name = sprintf('%s_%s', char(save_picture), base_name);
    end

    filename = fullfile(output_dir, base_name);
    exportgraphics(gcf, [filename '.pdf'], 'ContentType', 'vector');
    fprintf('Saved: %s.pdf\n', filename);
end
