function plot_var_tlp_method_metric_2x2(mat_files, p_value, metric, estimators, save_picture, close_picture)
% PLOT_VAR_TLP_METHOD_METRIC_2X2  Estimator-by-T method comparison
%
% Layout:
%   rows    = selected estimators
%   columns = sample sizes, normally T = 200 and T = 800
%   lines   = the four confidence-band constructions
%
% Example:
%   plot_var_tlp_method_metric_2x2( ...
%       {'Final15May_T200_sim1000_pL10_pV8_eta1_DGP_OleaSW_P1.mat', ...
%        'Final15May_T800_sim1000_pL10_pV8_eta1_DGP_OleaSW_P1.mat'}, ...
%       8, 'length', {'VAR', 'TLP'}, 1, false);

    if nargin < 6, close_picture = false; end
    if nargin < 5, save_picture = 0; end
    if nargin < 4 || isempty(estimators), estimators = {'VAR', 'TLP'}; end
    if nargin < 3 || isempty(metric), metric = 'length'; end
    if nargin < 2 || isempty(p_value), p_value = 8; end

    if ischar(estimators) || (isstring(estimators) && isscalar(estimators))
        estimators = {char(estimators)};
    elseif isstring(estimators)
        estimators = cellstr(estimators);
    elseif ~iscell(estimators)
        error('estimators must be a string or cell array of strings.');
    end
    estimators = cellfun(@(x) upper(char(x)), estimators, 'UniformOutput', false);

    if ~iscell(mat_files) || numel(mat_files) ~= 2
        error('mat_files must be a 1x2 cell array: {T200_file, T800_file}.');
    end
    if isempty(estimators)
        error('Provide at least one estimator, e.g. {''VAR''} or {''LP'', ''SLP'', ''VAR'', ''TLP''}.');
    end
    for i = 1:length(estimators)
        validate_estimator(estimators{i});
    end

    method_fields = {'method2', 'method1', 'method8', 'method7'};
    method_labels = { ...
        'Subtract VAR', ...
        'Subtract Bootstrap Mean', ...
        'Symmetric around VAR', ...
        'Subtract Boot. Mean + Symmetric'};

    data = cell(1, 2);
    for col = 1:2
        data{col} = load_plot_data(mat_files{col}, p_value);
    end

    metric_lower = lower(char(metric));
    [y_label, base_y_lim, y_ticks, ref_line, ref_label, ref_line_style] = ...
        get_metric_settings(metric_lower);

    y_lims_by_column = cell(1, 2);
    for col = 1:2
        if isempty(base_y_lim)
            y_lims_by_column{col} = compute_auto_y_limits(data, method_fields, ...
                metric_lower, estimators, col);
        else
            y_lims_by_column{col} = base_y_lim;
        end
    end

    [method_colors, method_line_styles, method_line_widths] = get_method_styles();
    ref_color = [0.5, 0.5, 0.5];
    line_width_ref = 3;

    axis_font = 16;
    label_font = 18;
    title_font = 20;
    legend_font = 12;
    panel_label_font = 16;

    n_rows = length(estimators);

    figure;
    fig_width = 0.37;
    fig_height = min(0.95, 0.18 * n_rows + 0.12);
    fig_bottom = max(0.025, 0.50 - fig_height / 2);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05 fig_bottom fig_width fig_height]);

    tiledlayout(n_rows, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    legend_built = false;
    legend_handles = [];
    legend_entries = {};

    for row = 1:n_rows
        estimator = estimators{row};

        for col = 1:2
            y_lim = y_lims_by_column{col};
            nexttile;
            hold on;

            h_ref = [];
            if ~isempty(ref_line)
                h_ref = yline(ref_line, 'LineWidth', line_width_ref, ...
                    'Color', ref_color, 'LineStyle', ref_line_style);
                uistack(h_ref, 'bottom');
            end

            line_handles_this_tile = [];
            max_horizon = 0;

            for method_idx = 1:length(method_fields)
                method = resolve_field_path(data{col}.stats, ...
                    [estimator '.Studentized.' method_fields{method_idx}]);
                [horizons, values] = extract_metric_values(method, metric_lower, data{col}.p_idx);
                max_horizon = max(max_horizon, max(horizons));

                h_line = plot(horizons, values, ...
                    'Color', method_colors(method_idx, :), ...
                    'LineWidth', method_line_widths(method_idx), ...
                    'LineStyle', method_line_styles{method_idx});
                line_handles_this_tile = [line_handles_this_tile, h_line]; %#ok<AGROW>
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
            set(gca, 'FontSize', axis_font);
            grid on;
            box on;

            current_ylim = ylim;
            y_pos = current_ylim(1) + 0.03 * (current_ylim(2) - current_ylim(1));
            text(1, y_pos, estimator, ...
                'FontSize', panel_label_font, ...
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
            if row == n_rows
                xlabel('Horizon', 'FontSize', label_font);
            end

            if ~legend_built
                legend_handles = line_handles_this_tile;
                legend_entries = method_labels;
                if ~isempty(ref_label) && ~isempty(h_ref)
                    legend_handles = [legend_handles, h_ref];
                    legend_entries = [legend_entries, {ref_label}];
                end
                legend_built = true;
            end
        end
    end

    if ~isempty(legend_handles)
        if strcmp(metric_lower, 'coverage')
            legend_num_columns = 3;
            if length(legend_handles) == 5
                hold on;
                h_blank = plot(NaN, NaN, 'LineStyle', 'none', 'Marker', 'none');
                hold off;
                legend_handles = [legend_handles(1:2), h_blank, legend_handles(3:4), legend_handles(5)];
                legend_entries = [legend_entries(1:2), {''}, legend_entries(3:4), legend_entries(5)];
            end
        else
            legend_num_columns = 2;
        end

        lgd = legend(legend_handles, legend_entries, ...
            'Orientation', 'horizontal', ...
            'FontSize', legend_font, ...
            'NumColumns', legend_num_columns);
        lgd.Layout.Tile = 'north';
        lgd.ItemTokenSize = [46, 14];
    end

    set(gcf, 'Color', [1 1 1]);
    align_columns_like_vertical_dual_plot();

    if should_save_figure(save_picture)
        save_metric_2x2_figure(data, p_value, save_picture, metric_lower, estimators);
    end

    if close_picture
        close(gcf);
    end
end

function [method_colors, method_line_styles, method_line_widths] = get_method_styles()
    % Colors identify CI construction, not estimator.
    % Purple/red match the original coverage-comparison palette exactly.
    method_colors = [
        0.55, 0.40, 0.60;  % Subtract VAR
        0.45, 0.45, 0.00;  % Subtract Bootstrap Mean
        0.00, 0.50, 0.50;  % Symmetric around VAR
        0.85, 0.35, 0.35;  % Subtract Boot. Mean + Symmetric
    ];

    method_line_styles = {'--', ':', '-.', '-'};
    method_line_widths = [4.2, 3.5, 3.5, 4.2];
end

function validate_estimator(estimator)
    if ~ismember(upper(estimator), {'LP', 'SLP', 'VAR', 'TLP'})
        error('Estimator must be LP, SLP, VAR, or TLP. Got "%s".', estimator);
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

function y_lim = compute_auto_y_limits(data, method_fields, metric, estimators, col)
    all_values = [];
    for row = 1:length(estimators)
        estimator = estimators{row};
        for method_idx = 1:length(method_fields)
            method = resolve_field_path(data{col}.stats, ...
                [estimator '.Studentized.' method_fields{method_idx}]);
            [~, values] = extract_metric_values(method, metric, data{col}.p_idx);
            all_values = [all_values, values]; %#ok<AGROW>
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

function save_metric_2x2_figure(data, p_value, save_picture, metric, estimators)
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
