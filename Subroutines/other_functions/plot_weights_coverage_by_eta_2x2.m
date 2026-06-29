function plot_weights_coverage_by_eta_2x2(T_values, eta_values, eta_file_tags, ...
    file_pattern, p_value, method_fields, save_picture, close_picture)
% PLOT_WEIGHTS_COVERAGE_BY_ETA_2X2  Paper-ready eta comparison figure.
%
% Layout:
%   row 1: mean TLP method7 weights by horizon, one line per eta
%   row 2: average coverage across horizons by eta, one line per method
%   columns: sample sizes, normally T = 200 and T = 800
%
% Example:
%   plot_weights_coverage_by_eta_2x2();

    if nargin < 8 || isempty(close_picture), close_picture = false; end
    if nargin < 7 || isempty(save_picture), save_picture = 1; end
    if nargin < 6 || isempty(method_fields)
        method_fields = {'LP.Studentized.method7', ...
                         'SLP.Studentized.method7', ...
                         'VAR.Studentized.method7', ...
                         'BLP', ...
                         'TLP.Studentized.method7'};
    end
    if nargin < 5 || isempty(p_value), p_value = 8; end
    if nargin < 4 || isempty(file_pattern)
        file_pattern = 'Final16June_T%d_sim1000_pL10_pV8_eta%s_DGP_OleaSW_P1.mat';
    end
    if nargin < 3 || isempty(eta_file_tags)
        eta_file_tags = {'1', '2', '4', '8', '16', '32', '64'};
    end
    if nargin < 2 || isempty(eta_values), eta_values = [1, 2, 4, 8, 16, 32, 64]; end
    if nargin < 1 || isempty(T_values), T_values = [200, 800]; end

    if numel(T_values) ~= 2
        error('T_values must contain exactly two sample sizes for the 2x2 plot.');
    end
    if numel(eta_values) ~= numel(eta_file_tags)
        error('eta_values and eta_file_tags must have the same length.');
    end
    if ischar(method_fields) || (isstring(method_fields) && isscalar(method_fields))
        method_fields = {char(method_fields)};
    elseif isstring(method_fields)
        method_fields = cellstr(method_fields);
    end

    data = cell(1, 2);
    for col = 1:2
        data{col} = load_eta_panel_data(T_values(col), eta_values, ...
            eta_file_tags, file_pattern, p_value, method_fields);
    end

    [eta_colors, eta_line_styles, eta_legend] = get_eta_styles(eta_values);
    [method_colors, method_line_styles, method_markers] = get_method_styles(method_fields);
    method_legend = method_labels_from_paths(method_fields);

    n_eta = numel(eta_values);
    figure;
    fig_width = min(0.72, 0.37 + 0.035 * max(0, n_eta - 6));
    fig_height = min(0.70, 0.52 + 0.025 * max(0, n_eta - 6));
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05 0.24 fig_width fig_height]);

    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    axis_font = 16;
    label_font = 18;
    title_font = 20;
    legend_font = max(10, 14 - max(0, n_eta - 6));
    line_width = 3.2;
    ref_line_width = 2;

    eta_handles = [];
    method_handles = [];

    for col = 1:2
        nexttile(col);
        hold on;
        for eta_idx = 1:numel(eta_values)
            h_line = plot(data{col}.horizons, data{col}.weight_by_eta(:, eta_idx), ...
                'Color', eta_colors(eta_idx, :), ...
                'LineStyle', eta_line_styles{eta_idx}, ...
                'LineWidth', line_width);
            if col == 1
                eta_handles = [eta_handles, h_line]; %#ok<AGROW>
            end
        end
        yline(0.5, 'Color', [0.45 0.45 0.45], ...
            'LineStyle', '--', 'LineWidth', ref_line_width);
        hold off;

        xlim([min(data{col}.horizons), max(data{col}.horizons)]);
        ylim([0, 1]);
        xticks(unique([0:4:max(data{col}.horizons), max(data{col}.horizons)]));
        yticks(0:0.25:1);
        grid on;
        box on;
        ax = gca;
        ax.YAxis.TickLabelFormat = '%.1f';
        set(gca, 'FontSize', axis_font, 'LineWidth', 1.2, 'Layer', 'top');

        title(sprintf('T = %d', round(T_values(col))), ...
            'FontSize', title_font, 'FontWeight', 'bold');
        xlabel('Horizon', 'FontSize', label_font);
        if col == 1
            ylabel({'Mean v(\lambda)', '(LP share)'}, 'FontSize', label_font);
        end
    end

    for col = 1:2
        nexttile(2 + col);
        hold on;
        for method_idx = 1:numel(method_fields)
            h_line = plot(eta_values, data{col}.coverage_by_eta(:, method_idx), ...
                'Color', method_colors(method_idx, :), ...
                'LineStyle', method_line_styles{method_idx}, ...
                'Marker', method_markers{method_idx}, ...
                'MarkerSize', 8, ...
                'MarkerFaceColor', method_colors(method_idx, :), ...
                'LineWidth', line_width);
            if col == 1
                method_handles = [method_handles, h_line]; %#ok<AGROW>
            end
        end
        yline(0.9, 'Color', [0.1 0.1 0.1], ...
            'LineStyle', '--', 'LineWidth', ref_line_width);
        hold off;

        set(gca, 'XScale', 'log', 'FontSize', axis_font, ...
            'LineWidth', 0.8, 'Layer', 'top');
        xlim([min(eta_values) * 0.9, max(eta_values) * 1.1]);
        ylim([0.5, 1.0]);
        xticks(eta_values);
        xticklabels(compose('%g', eta_values));
        if n_eta > 8
            xtickangle(35);
        end
        yticks(0.5:0.1:1.0);
        grid on;
        box on;
        ax = gca;
        ax.YAxis.TickLabelFormat = '%.1f';
        set(gca, 'LineWidth', 1.2, 'Layer', 'top');

        xlabel('\eta', 'FontSize', label_font);
        if col == 1
            ylabel('Average coverage', 'FontSize', label_font);
        end
    end

    if ~isempty(eta_handles)
        eta_legend_columns = min(4, numel(eta_values));
        lgd_eta = legend(eta_handles, eta_legend, ...
            'Orientation', 'horizontal', ...
            'NumColumns', eta_legend_columns, ...
            'FontSize', legend_font, ...
            'Box', 'on');
        lgd_eta.Layout.Tile = 'north';
        lgd_eta.ItemTokenSize = [46, 14];
    end

    if ~isempty(method_handles)
        lgd_method = legend(method_handles, method_legend, ...
            'Orientation', 'horizontal', ...
            'NumColumns', numel(method_fields), ...
            'FontSize', legend_font, ...
            'Box', 'on');
        lgd_method.Layout.Tile = 'south';
        lgd_method.ItemTokenSize = [46, 14];
    end

    set(gcf, 'Color', [1 1 1]);
    align_columns_like_vertical_dual_plot();

    if should_save_figure(save_picture)
        save_weights_coverage_figure(T_values, p_value, save_picture);
    end

    if close_picture
        close(gcf);
    end
end

function data = load_eta_panel_data(T_value, eta_values, eta_file_tags, ...
    file_pattern, p_value, method_fields)

    horizons = [];
    weight_by_eta = [];
    coverage_by_eta = nan(numel(eta_values), numel(method_fields));

    for eta_idx = 1:numel(eta_values)
        mat_file = sprintf(file_pattern, round(T_value), eta_file_tags{eta_idx});
        mat_file = resolve_mat_file(mat_file);
        S = load(mat_file, 'output', 'stats', 'str');

        if ~isfield(S, 'output') || ~isfield(S.output, 'TLP2') ...
                || ~isfield(S.output.TLP2, 'Studentized') ...
                || ~isfield(S.output.TLP2.Studentized, 'method7') ...
                || ~isfield(S.output.TLP2.Studentized.method7, 'v_lambda')
            error('Could not find output.TLP2.Studentized.method7.v_lambda in %s', mat_file);
        end

        v_raw = S.output.TLP2.Studentized.method7.v_lambda;
        p_idx = resolve_p_index(S, v_raw, p_value);
        if ndims(v_raw) >= 3
            v_eta = squeeze(v_raw(:, :, p_idx));
        else
            v_eta = squeeze(v_raw);
        end

        [h_eta, v_eta] = orient_horizon_matrix(S, v_eta);
        if isempty(horizons)
            horizons = h_eta;
            weight_by_eta = nan(numel(horizons), numel(eta_values));
        elseif numel(h_eta) ~= numel(horizons) || any(h_eta ~= horizons)
            error('Unexpected horizon grid in %s.', mat_file);
        end
        weight_by_eta(:, eta_idx) = mean_by_row_omitnan(v_eta);

        for method_idx = 1:numel(method_fields)
            method = resolve_field_path(S.stats, method_fields{method_idx});
            if isempty(method) || ~isfield(method, 'coverage')
                error('Could not find coverage for %s in %s.', ...
                    method_fields{method_idx}, mat_file);
            end
            coverage_h = extract_coverage_vector(method.coverage, S, p_value);
            coverage_by_eta(eta_idx, method_idx) = mean_omitnan(coverage_h(:));
        end
    end

    data = struct();
    data.horizons = horizons;
    data.weight_by_eta = weight_by_eta;
    data.coverage_by_eta = coverage_by_eta;
end

function [horizons, values] = orient_horizon_matrix(S, values)
    if isfield(S, 'str') && isfield(S.str, 'H_min') && isfield(S.str, 'H_max')
        expected_horizons = S.str.H_max - S.str.H_min + 1;
        h_min = S.str.H_min;
    else
        expected_horizons = 25;
        h_min = 0;
    end

    if size(values, 1) ~= expected_horizons && size(values, 2) == expected_horizons
        values = values';
    end

    horizons = h_min:(h_min + size(values, 1) - 1);
end

function coverage_h = extract_coverage_vector(coverage_raw, S, p_value)
    coverage_raw = squeeze(coverage_raw);
    p_idx = resolve_p_index(S, coverage_raw, p_value);

    if isvector(coverage_raw)
        coverage_h = coverage_raw(:);
    elseif size(coverage_raw, 2) >= p_idx
        coverage_h = coverage_raw(:, p_idx);
    elseif size(coverage_raw, 1) >= p_idx
        coverage_h = coverage_raw(p_idx, :)';
    else
        coverage_h = coverage_raw(:);
    end
end

function p_idx = resolve_p_index(S, values, p_value)
    p_idx = 1;
    if isfield(S, 'str') && isfield(S.str, 'P_VAR')
        p_values = S.str.P_VAR(:)';
        idx = find(abs(p_values - p_value) < 1e-10, 1);
        if ~isempty(idx)
            p_idx = idx;
            return;
        end
    end

    if ~isvector(values) && p_value >= 1 ...
            && abs(p_value - round(p_value)) < 1e-10 ...
            && round(p_value) <= max(size(values))
        p_idx = round(p_value);
    end
end

function method = resolve_field_path(s, path)
    if isstring(path), path = char(path); end
    parts = strsplit(path, '.');
    method = s;

    for i = 1:numel(parts)
        if isstruct(method) && isfield(method, parts{i})
            method = method.(parts{i});
        else
            method = [];
            return;
        end
    end
end

function values = mean_by_row_omitnan(x)
    values = nan(size(x, 1), 1);
    for row = 1:size(x, 1)
        values(row) = mean_omitnan(x(row, :));
    end
end

function value = mean_omitnan(x)
    x = x(~isnan(x));
    if isempty(x)
        value = NaN;
    else
        value = mean(x);
    end
end

function [eta_colors, eta_line_styles, eta_legend] = get_eta_styles(eta_values)
    n_eta = numel(eta_values);
    base_colors = [0.0000 0.4470 0.7410; ...
                   0.8500 0.3250 0.0980; ...
                   0.4660 0.6740 0.1880; ...
                   0.4940 0.1840 0.5560; ...
                   0.6350 0.0780 0.1840; ...
                   0.3010 0.7450 0.9330; ...
                   0.9290 0.6940 0.1250; ...
                   0.2500 0.2500 0.2500];

    if n_eta <= size(base_colors, 1)
        eta_colors = base_colors(1:n_eta, :);
    else
        eta_colors = lines(n_eta);
    end

    line_style_cycle = {'-', '--', '-.', ':'};
    eta_line_styles = cell(1, n_eta);
    for eta_idx = 1:n_eta
        eta_line_styles{eta_idx} = line_style_cycle{mod(eta_idx - 1, numel(line_style_cycle)) + 1};
    end

    eta_legend = cell(1, n_eta);
    for eta_idx = 1:n_eta
        eta_legend{eta_idx} = sprintf('\\eta = %g', eta_values(eta_idx));
    end
end

function [colors, line_styles, markers] = get_method_styles(method_fields)
    n_methods = numel(method_fields);
    colors = nan(n_methods, 3);
    line_styles = cell(1, n_methods);
    markers = cell(1, n_methods);

    for method_idx = 1:n_methods
        method_name = method_labels_from_paths(method_fields(method_idx));
        [colors(method_idx, :), line_styles{method_idx}, markers{method_idx}] = ...
            get_method_style(method_name{1});
    end
end

function [color, line_style, marker] = get_method_style(method_name)
    switch method_name
        case 'TLP'
            color = [0.9, 0.5, 0.1];
            line_style = '-';
            marker = 'o';
        case 'LP'
            color = [0.2, 0.6, 0.4];
            line_style = '--';
            marker = 's';
        case 'VAR'
            color = [0.8, 0.4, 0.8];
            line_style = ':';
            marker = '^';
        case 'SLP'
            color = [0.1, 0.6, 0.8];
            line_style = '-.';
            marker = 'd';
        case 'BLP'
            color = [0.2, 0.4, 0.6];
            line_style = '--';
            marker = 'v';
        otherwise
            color = [0.3, 0.3, 0.3];
            line_style = '-';
            marker = 'o';
    end
end

function labels = method_labels_from_paths(method_fields)
    labels = cell(1, numel(method_fields));
    for i = 1:numel(method_fields)
        if isstring(method_fields{i}), method_fields{i} = char(method_fields{i}); end
        parts = strsplit(method_fields{i}, '.');
        labels{i} = parts{1};
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

function resolved_file = resolve_mat_file(mat_file)
    if exist(mat_file, 'file') == 2
        resolved_file = mat_file;
        return;
    end

    function_dir = fileparts(mfilename('fullpath'));
    candidate = fullfile(function_dir, mat_file);
    if exist(candidate, 'file') == 2
        resolved_file = candidate;
        return;
    end

    error('Missing eta result file: %s', mat_file);
end

function value = should_save_figure(save_picture)
    value = false;
    if ischar(save_picture) || (isstring(save_picture) && isscalar(save_picture))
        value = true;
    elseif (isnumeric(save_picture) || islogical(save_picture)) && any(save_picture)
        value = true;
    end
end

function save_weights_coverage_figure(T_values, p_value, save_picture)
    output_dir = 'TablesAndPlots';
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    base_name = sprintf('weights_coverage_by_eta_2x2_T%d_T%d_pV%g_OleaSW_P1', ...
        round(T_values(1)), round(T_values(2)), p_value);
    if ischar(save_picture) || (isstring(save_picture) && isscalar(save_picture))
        base_name = sprintf('%s_%s', char(save_picture), base_name);
    end

    filename = fullfile(output_dir, [base_name '.pdf']);
    exportgraphics(gcf, filename, 'ContentType', 'vector');
    fprintf('Saved: %s\n', filename);
end
