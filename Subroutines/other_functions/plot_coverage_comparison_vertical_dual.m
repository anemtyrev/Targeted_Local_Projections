function plot_coverage_comparison_vertical_dual(method_pair_fields, mat_files, p_picture, ...
    save_picture, close_picture, metric)
% PLOT_COVERAGE_COMPARISON_VERTICAL_DUAL  Publication-ready 4x2 method comparison
%
%==========================================================================
% PURPOSE:
%   Creates a rowsx2 tiled figure comparing two inference methods within
%   each estimator (LP, SLP, VAR, TLP), with columns for different sample
%   sizes (e.g. T=200 vs T=800). Single shared legend at the top.
%
%   This function supports both:
%   1) Legacy combined files (one file per T, with multiple p values inside)
%   2) Split files (separate file per p value, e.g. pV1 and pV8 for each T)
%
%==========================================================================
% INPUTS:
%   method_pair_fields - Cell array of method pairs using dot-notation
%                        strings:
%                        {{'LP.Studentized.method7','LP.Studentized.method2'}, ...}
%
%   mat_files         - Flexible input:
%                       Legacy:
%                         {'T200_combined.mat','T800_combined.mat'}
%                       Split (flat):
%                         {'T200_pV1.mat','T200_pV8.mat',...
%                          'T800_pV1.mat','T800_pV8.mat'}
%                       Split (nested):
%                         {{'T200_pV1.mat','T200_pV8.mat'}, ...
%                          {'T800_pV1.mat','T800_pV8.mat'}}
%
%   p_picture         - p selectors. Can be:
%                       - Index selectors (legacy): [1 2]
%                       - Actual p values: [1 8]
%                       Single value behaves as before.
%
%   save_picture      - Save control:
%                       0/false: do not save
%                       1/true: save with auto-generated filename
%                       'prefix': save with custom prefix
%
%   close_picture     - Logical/numeric flag to close figure after saving
%
%   metric            - (Optional) 'coverage', 'length', 'bias', 'rmse', etc.
%                       Default: 'coverage'
%
%==========================================================================
% SEE ALSO:
%   plot_coverage_comparison_vertical, plot_stat_by_horizon_dual
%==========================================================================

    if nargin < 6 || isempty(metric), metric = 'coverage'; end
    if nargin < 5, close_picture = false; end
    if nargin < 4, save_picture = 0; end
    if nargin < 3, p_picture = 2; end

    %% Load input files into two plotting columns (T=left, T=right)
    column_data = prepare_column_data(mat_files);

    n_rows = length(method_pair_fields);

    %% Resolve p selectors to actual p values per column
    p_values_by_col = cell(1, 2);
    for col = 1:2
        p_values_by_col{col} = resolve_p_selectors_to_values(p_picture, column_data(col).available_p_values);
    end

    %% Settings
    p_compare_mode = (length(p_picture) > 1);

    % Row labels
    row_labels = {'LP', 'SLP', 'VAR', 'TLP'};
    if n_rows > length(row_labels)
        for i = length(row_labels)+1:n_rows
            row_labels{i} = sprintf('Method%d', i);
        end
    end

    % Which rows get 4 lines in p_compare_mode (VAR and TLP)
    four_line_rows = [3, 4];

    %% Metric-specific settings
    metric_lower = lower(metric);
    switch metric_lower
        case 'coverage'
            y_label = 'Coverage';
            y_lim = [0.5 1];
            y_ticks = [0.5 0.6 0.7 0.8 0.9 1.0];
            ref_line = 0.9;
            ref_label = '90 %';
        case 'length'
            y_label = 'Length';
            y_lim = [];
            y_ticks = [];
            ref_line = [];
            ref_label = '';
        case 'bias'
            y_label = 'Bias';
            y_lim = [];
            y_ticks = [];
            ref_line = 0;
            ref_label = 'Zero';
        case 'rmse'
            y_label = 'RMSE';
            y_lim = [];
            y_ticks = [];
            ref_line = [];
            ref_label = '';
        case 'mse'
            y_label = 'MSE';
            y_lim = [];
            y_ticks = [];
            ref_line = [];
            ref_label = '';
        case 'variance'
            y_label = 'Variance';
            y_lim = [];
            y_ticks = [];
            ref_line = [];
            ref_label = '';
        case 'bias_sq'
            y_label = 'Bias Squared';
            y_lim = [];
            y_ticks = [];
            ref_line = [];
            ref_label = '';
        otherwise
            y_label = metric;
            y_lim = [];
            y_ticks = [];
            ref_line = [];
            ref_label = '';
    end

    %% For auto y-limits, compute global range across all lines
    if isempty(y_lim)
        all_values = [];
        for col = 1:2
            for row = 1:n_rows
                p_vals_for_ylim = p_values_by_col{col};
                for pp = 1:length(p_vals_for_ylim)
                    [method1, method2, p_idx] = get_method_pair_for_p( ...
                        column_data(col), method_pair_fields{row}, p_vals_for_ylim(pp));

                    [~, vals1] = extract_metric_data_local(method1, metric_lower, p_idx);
                    [~, vals2] = extract_metric_data_local(method2, metric_lower, p_idx);
                    all_values = [all_values, vals1, vals2]; %#ok<AGROW>
                end
            end
        end
        if ~isempty(all_values)
            y_min = min(all_values);
            y_max = max(all_values);
            y_range = y_max - y_min;
            y_lim = [y_min - 0.05*y_range, y_max + 0.05*y_range];
            if strcmp(metric_lower, 'bias')
                y_lim(1) = min(y_lim(1), -0.05*y_range);
                y_lim(2) = max(y_lim(2), 0.05*y_range);
            end
        end
    end

    %% Colors (colorblind-safe, distinct from plot_stat_by_horizon palette)
    color_main = [0.85, 0.35, 0.35];
    color_alt  = [0.55, 0.40, 0.60];
    ref_color  = [0.5, 0.5, 0.5];

    line_width_solid = 5;
    line_width_alt = 3;
    line_width_dashed = 2.5;
    line_width_alt_dashed = 2.5;
    line_width_ref = 3;

    %% Font sizes
    axis_font = 18;
    label_font = 20;
    title_font = 22;
    legend_font = 14;

    %% Figure setup
    figure;
    fig_width = 0.37;
    fig_height = min(0.95, 0.14 * n_rows + 0.40);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05 0.025 fig_width fig_height]);

    tiledlayout(n_rows, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    %% Column headers
    T1_label = sprintf('T = %d', round(column_data(1).T));
    T2_label = sprintf('T = %d', round(column_data(2).T));

    %% p values for legend (use column 1)
    if p_compare_mode
        p_low = p_values_by_col{1}(1);
        p_high = p_values_by_col{1}(end);
    else
        p_single = p_values_by_col{1}(1);
    end

    %% Legend tracking
    legend_built = false;
    legend_handles = [];
    legend_entries = {};

    %% Plot loop: rows x columns
    for row = 1:n_rows
        for col = 1:2
            nexttile;
            hold on;

            % Determine if this row gets 4 lines
            is_four_line_row = p_compare_mode && ismember(row, four_line_rows);

            if is_four_line_row
                p_values_to_plot = p_values_by_col{col};
            else
                p_values_to_plot = p_values_by_col{col}(end);
            end

            max_horizon = 0;

            % Reference line
            h_ref = [];
            if ~isempty(ref_line)
                h_ref = yline(ref_line, 'LineWidth', line_width_ref, 'Color', ref_color);
                uistack(h_ref, 'bottom');
            end

            % Initialize handles
            h_main_high = [];
            h_main_low = [];
            h_alt_high = [];
            h_alt_low = [];

            % Plot lines for each p value
            for pp_idx = 1:length(p_values_to_plot)
                p_value = p_values_to_plot(pp_idx);

                [method1, method2, p_idx] = get_method_pair_for_p( ...
                    column_data(col), method_pair_fields{row}, p_value);

                [horizons1, values1] = extract_metric_data_local(method1, metric_lower, p_idx);
                [horizons2, values2] = extract_metric_data_local(method2, metric_lower, p_idx);

                if ~isempty(horizons1)
                    max_horizon = max(max_horizon, max(horizons1));
                end
                if ~isempty(horizons2)
                    max_horizon = max(max_horizon, max(horizons2));
                end

                % Line style based on low/high position in requested list
                if length(p_values_to_plot) == 1 || pp_idx == length(p_values_to_plot)
                    style_main = '-';
                    style_alt = '-';
                    lw_main = line_width_solid;
                    lw_alt = line_width_alt;
                else
                    style_main = '--';
                    style_alt = '--';
                    lw_main = line_width_dashed;
                    lw_alt = line_width_alt_dashed;
                end

                % Main method
                if ~isempty(values1)
                    h_tmp = plot(horizons1, values1, style_main, 'Color', color_main, 'LineWidth', lw_main);
                    if pp_idx == length(p_values_to_plot)
                        h_main_high = h_tmp;
                    else
                        h_main_low = h_tmp;
                    end
                end

                % Alternative method
                if ~isempty(values2)
                    h_tmp = plot(horizons2, values2, style_alt, 'Color', color_alt, 'LineWidth', lw_alt);
                    if pp_idx == length(p_values_to_plot)
                        h_alt_high = h_tmp;
                    else
                        h_alt_low = h_tmp;
                    end
                end
            end

            hold off;

            %% Formatting
            xlim([0 max_horizon]);
            if ~isempty(y_lim)
                ylim(y_lim);
            end

            base_ticks = 0:5:max_horizon;
            if ~isempty(base_ticks) && base_ticks(end) ~= max_horizon
                base_ticks = [base_ticks, max_horizon];
            end
            if ~isempty(base_ticks)
                xticks(base_ticks);
            end
            if ~isempty(y_ticks)
                yticks(y_ticks);
            end

            ax = gca;
            ax.YAxis.TickLabelFormat = '%.1f';
            set(gca, 'FontSize', axis_font);
            grid on;
            box on;

            % Estimator label in bottom left
            if row <= length(row_labels)
                label_text = sprintf('%-3s', row_labels{row});
                current_ylim = ylim;
                y_pos = current_ylim(1) + 0.03 * (current_ylim(2) - current_ylim(1));
                text(1, y_pos, label_text, 'FontSize', label_font + 2, 'FontWeight', 'normal', ...
                    'FontName', get(gca, 'FontName'), 'Color', [0.2 0.2 0.2], ...
                    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
                    'BackgroundColor', [1 1 1 0.8], 'EdgeColor', [0.6 0.6 0.6], ...
                    'Margin', 3);
            end

            % Column title on top row
            if row == 1
                if col == 1
                    title(T1_label, 'FontSize', title_font, 'FontWeight', 'bold');
                else
                    title(T2_label, 'FontSize', title_font, 'FontWeight', 'bold');
                end
            end

            % Y-label only on left column
            if col == 1
                ylabel(y_label, 'FontSize', label_font);
            end

            % X-label only on bottom row
            if row == n_rows
                xlabel('Horizon', 'FontSize', label_font);
            end

            %% Build legend from first tile only
            if ~legend_built
                if p_compare_mode
                    legend_handles = [];
                    legend_entries = {};

                    hold on;
                    if ~isempty(h_main_high)
                        legend_entries{end+1} = sprintf('Mean Symmetric Double Bootstrap (q=%g)', p_high); %#ok<AGROW>
                        legend_handles = [legend_handles, h_main_high]; %#ok<AGROW>
                    end
                    if ~isempty(h_alt_high)
                        legend_entries{end+1} = sprintf('Alternative (q=%g)', p_high); %#ok<AGROW>
                        legend_handles = [legend_handles, h_alt_high]; %#ok<AGROW>
                    end
                    if ~isempty(ref_label) && ~isempty(h_ref)
                        legend_entries{end+1} = ref_label; %#ok<AGROW>
                        legend_handles = [legend_handles, h_ref]; %#ok<AGROW>
                    end

                    h_dummy_main_low = plot(NaN, NaN, '--', 'Color', color_main, 'LineWidth', line_width_dashed);
                    legend_entries{end+1} = sprintf('Mean Symmetric Double Bootstrap (q=%g)', p_low); %#ok<AGROW>
                    legend_handles = [legend_handles, h_dummy_main_low]; %#ok<AGROW>

                    h_dummy_alt_low = plot(NaN, NaN, '--', 'Color', color_alt, 'LineWidth', line_width_dashed);
                    legend_entries{end+1} = sprintf('Alternative (q=%g)', p_low); %#ok<AGROW>
                    legend_handles = [legend_handles, h_dummy_alt_low]; %#ok<AGROW>
                    hold off;
                else
                    if ~isempty(h_main_high)
                        legend_entries{end+1} = sprintf('Mean Symmetric Double Bootstrap (q=%g)', p_single); %#ok<AGROW>
                        legend_handles = [legend_handles, h_main_high]; %#ok<AGROW>
                    end
                    if ~isempty(h_alt_high)
                        legend_entries{end+1} = sprintf('Alternative (q=%g)', p_single); %#ok<AGROW>
                        legend_handles = [legend_handles, h_alt_high]; %#ok<AGROW>
                    end
                    if ~isempty(ref_label) && ~isempty(h_ref)
                        legend_entries{end+1} = ref_label; %#ok<AGROW>
                        legend_handles = [legend_handles, h_ref]; %#ok<AGROW>
                    end
                end

                legend_built = true;
            end
        end
    end

    %% Shared legend at top
    if ~isempty(legend_handles)
        lgd = legend(legend_handles, legend_entries, ...
            'Orientation', 'horizontal', ...
            'FontSize', legend_font, ...
            'NumColumns', 3);
        lgd.Layout.Tile = 'north';
        lgd.ItemTokenSize = [40, 14];
    end

    set(gcf, 'Color', [1 1 1]);

    %% Force consistent column alignment
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

    %% Save
    if (ischar(save_picture) || isstring(save_picture)) || (isnumeric(save_picture) && save_picture)
        save_vertical_dual_figure(metric_lower, p_values_by_col{1}, ...
            column_data(1).representative_str, column_data(2).representative_str, save_picture);
    end

    if close_picture
        close(gcf);
    end
end

%% ========================================================================
%  HELPER FUNCTIONS
%% ========================================================================

function column_data = prepare_column_data(mat_files)
    file_groups = group_files_into_columns(mat_files);

    column_data = repmat(struct( ...
        'datasets', [], ...
        'available_p_values', [], ...
        'T', NaN, ...
        'representative_str', struct()), 1, 2);

    for col = 1:2
        files_for_col = file_groups{col};
        if isempty(files_for_col)
            error('No files provided for plotting column %d.', col);
        end

        datasets = repmat(struct('stats', [], 'str', [], 'file', '', 'p_values', []), 1, length(files_for_col));
        for i = 1:length(files_for_col)
            data = load_required_fields(files_for_col{i});
            datasets(i).stats = data.stats;
            datasets(i).str = data.str;
            datasets(i).file = files_for_col{i};
            datasets(i).p_values = get_numeric_vector_field(data.str, 'P_VAR');

            if isempty(datasets(i).p_values)
                error('File "%s" does not contain a valid str.P_VAR.', files_for_col{i});
            end
        end

        sort_keys = arrayfun(@(d) min(d.p_values), datasets);
        [~, sort_order] = sort(sort_keys);
        datasets = datasets(sort_order);

        all_p = [];
        for i = 1:length(datasets)
            all_p = [all_p, datasets(i).p_values]; %#ok<AGROW>
        end
        available_p = unique(all_p, 'sorted');

        t_candidates = NaN(1, length(datasets));
        for i = 1:length(datasets)
            t_candidates(i) = get_scalar_numeric_field(datasets(i).str, 'T');
        end
        T_val = t_candidates(find(~isnan(t_candidates), 1, 'first'));
        if isempty(T_val)
            T_val = extract_T_from_filename(files_for_col{1});
        end

        [~, idx_rep] = max(arrayfun(@(d) max(d.p_values), datasets));

        column_data(col).datasets = datasets;
        column_data(col).available_p_values = available_p;
        column_data(col).T = T_val;
        column_data(col).representative_str = datasets(idx_rep).str;
    end
end

function file_groups = group_files_into_columns(mat_files)
    file_list = flatten_file_list(mat_files);
    if length(file_list) < 2
        error('mat_files must contain at least two file paths.');
    end

    if length(file_list) == 2
        file_groups = {{file_list{1}}, {file_list{2}}};
        return;
    end

    T_vals = NaN(1, length(file_list));
    for i = 1:length(file_list)
        T_vals(i) = extract_T_from_filename(file_list{i});
    end

    if all(~isnan(T_vals))
        unique_T = unique(T_vals, 'sorted');
        if length(unique_T) ~= 2
            error(['Expected exactly two distinct T values in mat_files for a dual plot. ', ...
                   'Found T values: %s'], mat2str(unique_T));
        end

        file_groups = cell(1, 2);
        for col = 1:2
            idx = find(T_vals == unique_T(col));
            file_groups{col} = file_list(idx);
        end
        return;
    end

    if length(file_list) == 4
        file_groups = {file_list(1:2), file_list(3:4)};
        return;
    end

    error(['Could not infer T-grouping from mat_files. Provide either:', newline, ...
           '  1) Two files {T_left, T_right}, or', newline, ...
           '  2) Four files grouped by T, or', newline, ...
           '  3) Nested groups {{T_left_files...}, {T_right_files...}}']);
end

function file_list = flatten_file_list(input_value)
    file_list = {};
    if iscell(input_value)
        for i = 1:numel(input_value)
            nested = flatten_file_list(input_value{i});
            file_list = [file_list, nested]; %#ok<AGROW>
        end
    elseif ischar(input_value) || (isstring(input_value) && isscalar(input_value))
        file_list = {char(input_value)};
    else
        error('mat_files entries must be file path strings or nested cell arrays.');
    end
end

function T_val = extract_T_from_filename(file_path)
    T_val = NaN;
    [~, base_name, ~] = fileparts(file_path);
    tokens = regexp(base_name, '(?:^|[_-])T(\d+)(?:[_-]|$)', 'tokens', 'once');
    if ~isempty(tokens)
        T_val = str2double(tokens{1});
    end
end

function data = load_required_fields(file_path)
    fprintf('Loading: %s\n', file_path);
    data = load(file_path);
    if ~isfield(data, 'stats') || ~isfield(data, 'str')
        error('File "%s" must contain variables "stats" and "str".', file_path);
    end
end

function p_values = resolve_p_selectors_to_values(p_selectors, available_p_values)
    if isempty(available_p_values)
        error('No available p values were found in the provided files.');
    end

    selectors = p_selectors(:)';
    available = available_p_values(:)';

    if length(available) == 1
        only_p = available(1);

        is_index = all(abs(selectors - round(selectors)) < 1e-10) && all(selectors == 1);
        is_actual = all(abs(selectors - only_p) < 1e-10);

        if is_index || is_actual
            p_values = repmat(only_p, size(selectors));
            return;
        end

        error(['Requested p selector(s) %s, but available P_VAR is only %s. ', ...
               'For a single-p file, use selector 1 or actual value %g.'], ...
              mat2str(selectors), mat2str(available), only_p);
    end

    is_actual = all(arrayfun(@(x) any(abs(available - x) < 1e-10), selectors));
    is_index = all(abs(selectors - round(selectors)) < 1e-10) && ...
               all(selectors >= 1 & selectors <= length(available));

    if is_actual && ~is_index
        p_values = selectors;
        return;
    end
    if is_index && ~is_actual
        p_values = available(selectors);
        return;
    end
    if is_actual && is_index
        % Ambiguous case (e.g. selector=1 and available starts at 1):
        % keep legacy behavior and interpret as index.
        p_values = available(selectors);
        return;
    end

    error(['Could not resolve p selectors %s against available P_VAR %s.', newline, ...
           'Use index selectors (e.g. [1 2]) or actual p values (e.g. [1 8]).'], ...
          mat2str(selectors), mat2str(available));
end

function [method1, method2, p_idx] = get_method_pair_for_p(column_entry, method_pair_field, p_value)
    dataset_idx = [];
    p_idx = [];

    for i = 1:length(column_entry.datasets)
        p_vals = column_entry.datasets(i).p_values;
        idx = find(abs(p_vals - p_value) < 1e-10, 1);
        if ~isempty(idx)
            dataset_idx = i;
            p_idx = idx;
            break;
        end
    end

    if isempty(dataset_idx)
        error('Requested p=%g not found. Available P_VAR values: %s', ...
            p_value, mat2str(column_entry.available_p_values));
    end

    stats_struct = column_entry.datasets(dataset_idx).stats;
    method1 = resolve_field_path(stats_struct, method_pair_field{1});
    method2 = resolve_field_path(stats_struct, method_pair_field{2});
end

function values = get_numeric_vector_field(s, field_name)
    values = [];
    if ~isstruct(s) || ~isfield(s, field_name)
        return;
    end

    raw = s.(field_name);
    if isnumeric(raw)
        values = unique(raw(:)', 'sorted');
    end
end

function value = get_scalar_numeric_field(s, field_name)
    value = NaN;
    if ~isstruct(s) || ~isfield(s, field_name)
        return;
    end
    raw = s.(field_name);
    if isnumeric(raw) && ~isempty(raw)
        value = raw(1);
    end
end

function result = resolve_field_path(s, path)
    result = s;
    parts = strsplit(path, '.');
    for i = 1:length(parts)
        if isfield(result, parts{i})
            result = result.(parts{i});
        else
            warning('Field path "%s" not found (failed at "%s").', path, parts{i});
            result = [];
            return;
        end
    end
end

function [horizons, values] = extract_metric_data_local(method, metric, p_idx)
    horizons = [];
    values = [];

    field_variants = {metric, lower(metric), upper(metric), ...
                      [upper(metric(1)) lower(metric(2:end))]};

    if strcmpi(metric, 'rmse')
        field_variants = {'rmse', 'RMSE', 'Rmse'};
    elseif strcmpi(metric, 'mse')
        field_variants = {'mse', 'MSE', 'Mse'};
    elseif strcmpi(metric, 'bias_sq')
        field_variants = {'bias_sq', 'Bias_sq', 'BIAS_SQ', 'bias_squared'};
    elseif strcmpi(metric, 'below_ci')
        field_variants = {'below_CI', 'below_ci', 'belowCI', 'BELOW_CI'};
    elseif strcmpi(metric, 'above_ci')
        field_variants = {'above_CI', 'above_ci', 'aboveCI', 'ABOVE_CI'};
    end

    found = false;
    field_name = '';
    for v = 1:length(field_variants)
        if isfield(method, field_variants{v})
            field_name = field_variants{v};
            found = true;
            break;
        end
    end

    if ~found && isfield(method, 'B')
        for v = 1:length(field_variants)
            if isfield(method.B, field_variants{v})
                method = method.B;
                field_name = field_variants{v};
                found = true;
                break;
            end
        end
    end

    if ~found
        return;
    end

    data_all = method.(field_name);

    if ndims(data_all) == 3
        if p_idx <= size(data_all, 3)
            data = squeeze(data_all(:, :, p_idx));
            values = data(:)';
        end
    elseif ismatrix(data_all)
        [nrows, ncols] = size(data_all);
        if nrows == 1
            values = data_all;
        elseif ncols == 1
            values = data_all(:)';
        elseif ncols >= p_idx
            values = data_all(:, p_idx)';
        end
    end

    if ~isempty(values)
        horizons = 0:(length(values) - 1);
    end
end

function save_vertical_dual_figure(metric, p_values_plot, str1, str2, save_picture)
    T1 = str1.T;
    T2 = str2.T;

    if isfield(str1, 'OleaDGP') && str1.OleaDGP == 1
        dgp = 3;
    elseif isfield(str1, 'misspec_VARMA') && str1.misspec_VARMA == 1
        dgp = 2;
    else
        dgp = 1;
    end

    if length(p_values_plot) > 1
        p_str = 'p_COMPARE';
    else
        p_str = sprintf('p%g', p_values_plot(1));
    end

    if ~exist('TablesAndPlots', 'dir')
        mkdir('TablesAndPlots');
    end

    base_name = sprintf('dgp_%d_DUAL_%s_VERTICAL_t%d_t%d_%s', ...
        dgp, upper(metric), T1, T2, p_str);
    if (ischar(save_picture) || isstring(save_picture))
        base_name = sprintf('%s_%s', char(save_picture), base_name);
    end

    filename = fullfile('TablesAndPlots', base_name);
    exportgraphics(gcf, [filename '.pdf'], 'ContentType', 'vector');
    fprintf('Saved: %s.pdf\n', filename);
end
