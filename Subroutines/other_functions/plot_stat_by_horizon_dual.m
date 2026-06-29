function plot_stat_by_horizon_dual(stat_of_interest, mat_files, highlighted_field, other_fields, ...
    p_picture, save_picture, close_picture, plot_markers, legend_names)
% PLOT_STAT_BY_HORIZON_DUAL  Publication-ready k×2 comparison across two sample sizes
%
%==========================================================================
% PURPOSE:
%   Creates a k×2 tiled figure where each row is a metric and each column
%   corresponds to a different .mat file (typically T=200 vs T=800).
%   A single shared legend sits at the top, eliminating duplication.
%
%==========================================================================
% SYNTAX:
%   plot_stat_by_horizon_dual(stat_of_interest, mat_files, highlighted_field, ...
%       other_fields, p_picture, save_picture, close_picture, plot_markers, legend_names)
%
%==========================================================================
% INPUTS:
%   stat_of_interest  - Cell array of metric names (one per row)
%                       Supported: 'coverage', 'length', 'bias', 'bias_sq',
%                       'bias_abs', 'bias_sqrt', 'variance', 'std', 'mse',
%                       'rmse', 'below_ci', 'above_ci'
%                       Example: {'coverage', 'length', 'bias_sq', 'variance', 'MSE'}
%
%   mat_files         - Cell array of two .mat file paths {file_left, file_right}
%                       Each must contain 'stats' and 'str' variables
%                       Example: {'NewTLP_T200_...mat', 'NewTLP_T800_...mat'}
%
%   highlighted_field - String specifying path to highlighted method in stats
%                       Uses dot notation: 'TLP.Studentized.method7'
%                       Set to '' or [] to skip highlighting
%
%   other_fields      - Cell array of dot-notation strings for comparison methods
%                       Example: {'LP.Studentized.method7', 'VAR.Studentized.method7',
%                                 'BLP', 'SLP.Studentized.method7'}
%
%   p_picture         - Scalar index into str.P_VAR for which lag order to plot
%
%   save_picture      - Save control:
%                       0 or false: do not save
%                       1 or true: save with auto-generated filename
%                       'prefix': save with custom prefix before auto filename
%                       Example: 'myexp' saves as TablesAndPlots/myexp_dgp_3_DUAL_...pdf
%
%   close_picture     - Logical/numeric flag to close figure after saving
%
%   plot_markers      - Logical/numeric flag for markers on lines (default: 0)
%
%   legend_names      - Legend control, same syntax as plot_stat_by_horizon:
%                       0: no legend
%                       1: auto names from method_name fields
%                       {'TLP','LP','VAR','BLP','SLP'}: custom method names
%                       {{'TLP','LP',...}, 1}: custom names + size mode
%                       (reference lines like '90%' are added automatically)
%
%==========================================================================
% USAGE EXAMPLE:
%
%   p_picture = [1];
%   save_picture = 0;
%   close_picture = 0;
%   plot_markers = 0;
%   legend_names = {{"Targeted Local Projection","Local Projection",...
%       "Vector Autoregression","Bayesian Local Projection",...
%       "Smooth Local Projection"},1};
%   stat_of_interest = {'coverage', 'length', 'bias_sq', 'variance', 'MSE'};
%
%   plot_stat_by_horizon_dual(stat_of_interest, ...
%       {'NewTLP_T200_sim1000_pL10_pV8_DGP_OleaSW_P1.mat', ...
%        'NewTLP_T800_sim1000_pL10_pV8_DGP_OleaSW_P1.mat'}, ...
%       'TLP.Studentized.method7', ...
%       {'LP.Studentized.method7', 'VAR.Studentized.method7', ...
%        'BLP', 'SLP.Studentized.method7'}, ...
%       p_picture, save_picture, close_picture, plot_markers, legend_names);
%
%==========================================================================
% OUTPUT LAYOUT:
%           [  shared horizontal legend across top  ]
%              T = 200          T = 800
%           +--------------+--------------+
%           |  coverage    |  coverage    |
%           +--------------+--------------+
%           |  length      |  length      |
%           +--------------+--------------+
%           |  bias_sq     |  bias_sq     |
%           +--------------+--------------+
%           |  variance    |  variance    |
%           +--------------+--------------+
%           |  MSE         |  MSE         |
%           +--------------+--------------+
%
%   - Y-axis labels on left column only
%   - X-axis label ("Horizon") on bottom row only
%   - Column headers ("T = 200", "T = 800") as titles on top row
%   - Single shared legend at top of figure
%
%==========================================================================
% SAVING CONVENTION:
%   TablesAndPlots/dgp_{D}_DUAL_{METRICS}_t{T1}_t{T2}_p{P}.jpeg
%
%==========================================================================
% SEE ALSO:
%   plot_stat_by_horizon, plot_coverage_comparison_vertical
%==========================================================================

    %% Default arguments
    if nargin < 9 || isempty(legend_names)
        legend_names = 1;
    end
    if nargin < 8 || isempty(plot_markers)
        plot_markers = 0;
    end
    if nargin < 7
        close_picture = 0;
    end
    if nargin < 6
        save_picture = 0;
    end
    if nargin < 5
        p_picture = 1;
    end

    % Convert metrics to cell array if needed
    if ~iscell(stat_of_interest)
        stat_of_interest = {stat_of_interest};
    end
    other_fields = normalize_field_specs(other_fields);
    highlighted_field = normalize_highlighted_field(highlighted_field);

    %% Load mat files
    data1 = load_dual_plot_mat(mat_files{1});
    data2 = load_dual_plot_mat(mat_files{2});

    stats_all = {data1.stats, data2.stats};
    str_all   = {data1.str,   data2.str};

    %% Resolve methods and IV flags from both files
    highlighted = cell(1, 2);
    others      = cell(1, 2);
    iv_mode_by_col = false(1, 2);
    highlighted_names = cell(1, 2);
    other_names_by_col = cell(1, 2);

    default_bases = {'LP', 'VAR', 'BLP', 'SLP', 'Method5', 'Method6', 'Method7', 'Method8'};
    for col = 1:2
        iv_mode_by_col(col) = isfield(str_all{col}, 'is_iv_setup') && ...
            ~isempty(str_all{col}.is_iv_setup) && logical(str_all{col}.is_iv_setup);

        if is_resolvable_field_spec(highlighted_field)
            highlighted{col} = resolve_field_path(stats_all{col}, highlighted_field);
        else
            highlighted{col} = [];
        end

        others{col} = cell(1, length(other_fields));
        other_names_by_col{col} = cell(1, length(other_fields));
        for m = 1:length(other_fields)
            others{col}{m} = resolve_field_path(stats_all{col}, other_fields{m});
            if m <= length(default_bases)
                fallback = default_bases{m};
            else
                fallback = sprintf('Method%d', m);
            end
            other_names_by_col{col}{m} = get_method_name( ...
                others{col}{m}, extract_base_name(other_fields{m}, fallback), iv_mode_by_col(col));
        end

        if ~isempty(highlighted{col})
            highlighted_names{col} = get_method_name( ...
                highlighted{col}, extract_base_name(highlighted_field), iv_mode_by_col(col));
        else
            highlighted_names{col} = '';
        end
    end

    if iv_mode_by_col(1) ~= iv_mode_by_col(2)
        warning('plot_stat_by_horizon_dual:MixedIVModes', ...
            ['Left and right MAT files differ in IV mode. Shared legend entries ' ...
             'are inferred from the left file.']);
    end

    %% Collect shared legend/color names (from left file)
    highlighted_name = highlighted_names{1};
    other_names = other_names_by_col{1};

    %% Parse legend_names (same logic as plot_stat_by_horizon)
    use_custom_legend = false;
    show_legend = true;
    legend_size_mode = 1;
    custom_names = {};

    if iscell(legend_names)
        if length(legend_names) == 2 && (iscell(legend_names{1}) || isempty(legend_names{1})) && isnumeric(legend_names{2})
            legend_size_mode = legend_names{2};
            if ~isempty(legend_names{1})
                custom_names = legend_names{1};
                use_custom_legend = true;
            end
            if legend_size_mode == 0
                show_legend = false;
            end
        else
            custom_names = legend_names;
            use_custom_legend = true;
        end

        if use_custom_legend
            has_highlighted_field = is_resolvable_field_spec(highlighted_field);
            n_methods = length(other_fields) + has_highlighted_field;
            if length(custom_names) ~= n_methods
                warning('Custom legend has %d entries but %d methods. Using auto names.', ...
                    length(custom_names), n_methods);
                use_custom_legend = false;
            end
        end
    elseif isnumeric(legend_names)
        legend_size_mode = legend_names;
        if legend_names == 0
            show_legend = false;
        end
    end

    %% Colors and styles
    base_colors = containers.Map();
    base_colors('TLP') = [0.9, 0.5, 0.1];
    base_colors('LP')  = [0.2, 0.6, 0.4];
    base_colors('VAR') = [0.8, 0.4, 0.8];
    base_colors('SLP') = [0.1, 0.6, 0.8];
    base_colors('BLP') = [0.2, 0.4, 0.6];

    color_palette = [
        0.9, 0.5, 0.1;  0.2, 0.6, 0.4;  0.8, 0.4, 0.8;
        0.1, 0.6, 0.8;  0.2, 0.4, 0.6;  0.8, 0.2, 0.2;
        0.5, 0.5, 0.1;  0.1, 0.4, 0.3;  0.6, 0.2, 0.6;
        0.3, 0.7, 0.5;
    ];

    linestyle_cycle = {'-', '--', ':', '-.'};
    symbols = {'s', 'o', '^', 'v', 'd'};
    size_marker = 400;

    if legend_size_mode == 2
        line_width   = 6;
        marker_scale = 1.3;
    elseif legend_size_mode == 3 || legend_size_mode == -3
        line_width   = 2.5;
        marker_scale = 0.8;
    else
        line_width   = 4;
        marker_scale = 1.0;
    end
    line_width_highlight = line_width;
    line_width_other = line_width;

    % All method names for color assignment
    all_names = other_names;
    if ~isempty(highlighted_name)
        all_names = [{highlighted_name}, other_names];
    end

    % Pre-compute colors
    method_colors = cell(1, length(all_names));
    for i = 1:length(all_names)
        method_colors{i} = get_method_color(all_names{i}, base_colors, color_palette, i, all_names);
    end

    %% Metric labels
    metric_labels = containers.Map();
    metric_labels('coverage')  = 'Coverage';
    metric_labels('below_ci')  = 'Miss Rate (Below CI)';
    metric_labels('above_ci')  = 'Miss Rate (Above CI)';
    metric_labels('length')    = 'Length';
    metric_labels('rmse')      = 'RMSE';
    metric_labels('variance')  = 'Variance';
    metric_labels('std')       = 'Std Dev';
    metric_labels('bias_sq')   = 'Bias Squared';
    metric_labels('bias')      = 'Bias';
    metric_labels('bias_abs')  = 'Absolute Bias';
    metric_labels('bias_sqrt') = 'Sqrt Bias';
    metric_labels('mse')       = 'Mean Squared Error';

    %% Figure setup
    n_metrics = length(stat_of_interest);

    figure;
    fig_width = 0.37;  % 2x single (0.18) + small column gap
    fig_height = min(0.95, 0.14 * n_metrics + 0.40);  % Same as plot_stat_by_horizon publication mode
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05 0.025 fig_width fig_height]);

    tl = tiledlayout(n_metrics, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % We'll adjust horizontal gap after all tiles are created

    %% Font sizes (matched to plot_stat_by_horizon publication mode)
    axis_font = 18;
    label_font = 20;
    title_font = 22;
    legend_font = 14;

    %% Column headers
    T1_label = sprintf('T = %d', str_all{1}.T);
    T2_label = sprintf('T = %d', str_all{2}.T);

    %% Build legend from first tile only
    legend_handles = [];
    legend_entries = {};
    legend_built = false;

    first_p = p_picture(1);

    %% Plot loop: metrics x columns
    has_shared_highlight = ~isempty(highlighted_name);
    for metric_idx = 1:n_metrics
        metric = lower(stat_of_interest{metric_idx});

        if ~isKey(metric_labels, metric)
            error('Unknown metric: %s', metric);
        end
        ylabel_text = metric_labels(metric);

        % Reference line info
        [ref_val, ref_label, ref_style] = get_reference_line(metric);

        for col = 1:2
            nexttile;
            hold on;
            highlighted_name_col = highlighted_names{col};
            other_names_col = other_names_by_col{col};

            % Reference line (behind data)
            h_ref = [];
            if ~isempty(ref_val)
                if strcmp(ref_style, '--')
                    h_ref = yline(ref_val, 'LineWidth', line_width, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
                else
                    h_ref = yline(ref_val, 'LineWidth', line_width, 'Color', [0.5 0.5 0.5]);
                end
            end

            max_horizon = 0;
            method_legend_idx = 1;
            tile_handles = [];

            % --- Highlighted method ---
            if ~isempty(highlighted{col})
                [horizons, values] = extract_metric_data(highlighted{col}, metric, first_p);
                if ~isempty(values)
                    max_horizon = max(max_horizon, max(horizons));
                    hl_linestyle = get_linestyle_from_name(highlighted_name_col);
                    if plot_markers
                        h = plot(horizons, values, ...
                            'Color', method_colors{1}, 'LineWidth', line_width_highlight, ...
                            'LineStyle', hl_linestyle, ...
                            'Marker', symbols{1}, 'MarkerSize', sqrt(size_marker) * 1.3 * marker_scale);
                    else
                        h = plot(horizons, values, ...
                            'Color', method_colors{1}, 'LineWidth', line_width_highlight, ...
                            'LineStyle', hl_linestyle);
                    end
                    tile_handles = [tile_handles h];

                    if ~legend_built
                        if use_custom_legend
                            legend_entries{end+1} = custom_names{method_legend_idx};
                        else
                            legend_entries{end+1} = highlighted_name;
                        end
                        legend_handles = [legend_handles h];
                    end
                end
                method_legend_idx = method_legend_idx + 1;
            end

            % --- Other methods ---
            linestyle_idx = 1;
            for m = 1:length(others{col})
                method = others{col}{m};
                color_idx = m + has_shared_highlight;

                [horizons, values] = extract_metric_data(method, metric, first_p, max_horizon);
                if isempty(values) && first_p ~= 1
                    [horizons, values] = extract_metric_data(method, metric, 1, max_horizon);
                end
                if isempty(values)
                    method_legend_idx = method_legend_idx + 1;
                    linestyle_idx = linestyle_idx + 1;
                    continue;
                end

                max_horizon = max(max_horizon, max(horizons));
                name_lower = lower(char(other_names_col{m}));
                if contains(name_lower, 'clp')
                    ls = get_linestyle_from_name(other_names_col{m});
                else
                    ls = linestyle_cycle{mod(linestyle_idx - 1, length(linestyle_cycle)) + 1};
                end
                linestyle_idx = linestyle_idx + 1;

                if plot_markers
                    h = plot(horizons, values, ...
                        'Color', method_colors{color_idx}, 'LineWidth', line_width_other, ...
                        'LineStyle', ls, ...
                        'Marker', symbols{mod(m, length(symbols)) + 1}, ...
                        'MarkerSize', sqrt(size_marker) * marker_scale);
                else
                    h = plot(horizons, values, ...
                        'Color', method_colors{color_idx}, 'LineWidth', line_width_other, ...
                        'LineStyle', ls);
                end
                tile_handles = [tile_handles h];

                if ~legend_built
                    if use_custom_legend
                        legend_entries{end+1} = custom_names{method_legend_idx};
                    else
                        legend_entries{end+1} = other_names{m};
                    end
                    legend_handles = [legend_handles h];
                end

                method_legend_idx = method_legend_idx + 1;
            end

            % Reference line legend entry (first tile only)
            if ~legend_built && ~isempty(ref_label) && ~isempty(h_ref)
                legend_entries{end+1} = ref_label;
                legend_handles = [legend_handles h_ref];
            end

            if ~legend_built
                legend_built = true;
            end

            hold off;

            %% Tile formatting
            if max_horizon > 0
                xlim([0 max_horizon]);
                base_ticks = 0:5:max_horizon;
                if base_ticks(end) ~= max_horizon
                    base_ticks = [base_ticks, max_horizon];
                end
                xticks(base_ticks);
            end

            set_ylim(metric);

            set(gca, 'FontSize', axis_font);
            grid on;
            box on;
            format_axis_tick_labels(gca);

            % Column title on top row
            if metric_idx == 1
                if col == 1
                    title(T1_label, 'FontSize', title_font, 'FontWeight', 'bold');
                else
                    title(T2_label, 'FontSize', title_font, 'FontWeight', 'bold');
                end
            end

            % Y-label only on left column
            if col == 1
                ylabel(ylabel_text, 'FontSize', label_font);
            end

            % X-label only on bottom row
            if metric_idx == n_metrics
                xlabel('Horizon', 'FontSize', label_font);
            end

            % Highlighted method on top
            if ~isempty(highlighted{col}) && ~isempty(tile_handles)
                uistack(tile_handles(1), 'top');
            end
        end
    end

    %% Shared legend at top
    if show_legend && ~isempty(legend_handles)
        n_entries = length(legend_entries);
        n_cols_legend = ceil(n_entries / 2);  % 2 rows
        lgd = legend(legend_handles, legend_entries, ...
            'Orientation', 'horizontal', ...
            'FontSize', legend_font, ...
            'NumColumns', n_cols_legend);
        lgd.Layout.Tile = 'north';
        lgd.ItemTokenSize = [40, 14];
    end

    set(gcf, 'Color', [1 1 1]);
    
    % Force consistent column alignment across all rows
    % Collect all axes (excluding legend)
    all_axes = flipud(findobj(gcf, 'Type', 'axes'));  % flipud to get tile order
    % Filter out legend axes
    real_axes = [];
    for i = 1:length(all_axes)
        if isempty(all_axes(i).Tag) || ~strcmp(all_axes(i).Tag, 'legend')
            real_axes = [real_axes; all_axes(i)];
        end
    end
    
    % Separate into left and right columns based on position
    left_axes = [];
    right_axes = [];
    for i = 1:length(real_axes)
        pos = real_axes(i).Position;
        if pos(1) < 0.5
            left_axes = [left_axes; real_axes(i)];
        else
            right_axes = [right_axes; real_axes(i)];
        end
    end
    
    % Force all left-column axes to share the same x-position and width
    if ~isempty(left_axes)
        max_left = max(arrayfun(@(ax) ax.Position(1), left_axes));
        min_width_L = min(arrayfun(@(ax) ax.Position(3), left_axes));
        for i = 1:length(left_axes)
            pos = left_axes(i).Position;
            left_axes(i).Position = [max_left, pos(2), min_width_L, pos(4)];
        end
    end
    
    % Force all right-column axes to share the same x-position and width
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
        save_dual_figure(stat_of_interest, p_picture, str_all{1}, str_all{2}, save_picture);
    end

    if close_picture
        close(gcf);
    end
end

%% ========================================================================
%  HELPER FUNCTIONS
%% ========================================================================

function data = load_dual_plot_mat(mat_file)
% Load a plotting MAT file and validate the required stats/str variables.
    if isstring(mat_file)
        if ~isscalar(mat_file)
            error('load_dual_plot_mat:InvalidInput', ...
                'mat_file must be a string scalar or character vector.');
        end
        mat_file = char(mat_file);
    elseif ~ischar(mat_file)
        error('load_dual_plot_mat:InvalidInput', ...
            'mat_file must be a string scalar or character vector.');
    end

    mat_file = strtrim(mat_file);
    if isempty(mat_file)
        error('load_dual_plot_mat:EmptyPath', 'mat_file cannot be empty.');
    end

    if exist(mat_file, 'file') ~= 2
        [~, ~, ext] = fileparts(mat_file);
        if isempty(ext)
            mat_candidate = [mat_file '.mat'];
            if exist(mat_candidate, 'file') == 2
                mat_file = mat_candidate;
            else
                error('load_dual_plot_mat:FileNotFound', ...
                    'Could not find MAT file "%s" (also tried "%s").', mat_file, mat_candidate);
            end
        else
            error('load_dual_plot_mat:FileNotFound', ...
                'Could not find MAT file "%s".', mat_file);
        end
    end

    fprintf('Loading: %s\n', mat_file);
    data = load(mat_file);

    if ~isfield(data, 'stats') || ~isfield(data, 'str')
        error('load_dual_plot_mat:MissingFields', ...
            'MAT file "%s" must contain variables ''stats'' and ''str''.', mat_file);
    end
end

function result = resolve_field_path(s, path)
% Navigate nested struct: 'TLP.Studentized.method7' -> s.TLP.Studentized.method7
    if isempty(path)
        result = [];
        return;
    end

    if isstring(path)
        if isscalar(path)
            path = char(path);
        else
            warning('Field path string array has %d elements; using the first one.', numel(path));
            path = char(path(1));
        end
    end

    if ~ischar(path)
        warning('Field path must be text, got %s. Returning empty.', class(path));
        result = [];
        return;
    end

    candidates = build_compatible_paths(path);
    for k = 1:length(candidates)
        [ok, resolved] = try_resolve_exact(s, candidates{k});
        if ok
            result = resolved;
            return;
        end
    end

    warning('Field path "%s" not found under legacy/new Studentized naming variants.', path);
    result = [];
end

function paths = build_compatible_paths(path)
% Build compatible field-path alternatives for Studentized/Studentised variants.
    parts = strsplit(path, '.');
    expanded = {parts};

    for i = 1:length(parts)
        options = token_compatibility(parts{i});
        if length(options) == 1
            continue;
        end

        new_expanded = cell(1, 0);
        for j = 1:length(expanded)
            for o = 1:length(options)
                candidate = expanded{j};
                candidate{i} = options{o};
                new_expanded{end + 1} = candidate; %#ok<AGROW>
            end
        end
        expanded = new_expanded;
    end

    paths = cell(1, length(expanded));
    for i = 1:length(expanded)
        paths{i} = strjoin(expanded{i}, '.');
    end

    [~, keep_idx] = unique(paths, 'stable');
    paths = paths(sort(keep_idx));
end

function options = token_compatibility(token)
    switch token
        case 'Studentized_DB'
            options = {'Studentized_DB', 'Studentised_DB', 'Studentized', 'Studentised'};
        case 'Studentized_SB'
            options = {'Studentized_SB', 'Studentised_SB', 'Studentized', 'Studentised'};
        case 'Studentized'
            options = {'Studentized', 'Studentised', 'Studentized_DB', 'Studentized_SB', 'Studentised_DB', 'Studentised_SB'};
        case 'Studentised_DB'
            options = {'Studentised_DB', 'Studentized_DB', 'Studentised', 'Studentized'};
        case 'Studentised_SB'
            options = {'Studentised_SB', 'Studentized_SB', 'Studentised', 'Studentized'};
        case 'Studentised'
            options = {'Studentised', 'Studentized', 'Studentised_DB', 'Studentised_SB', 'Studentized_DB', 'Studentized_SB'};
        otherwise
            options = {token};
    end
end

function [ok, result] = try_resolve_exact(s, path)
    result = s;
    ok = true;
    parts = strsplit(path, '.');
    for i = 1:length(parts)
        if ~isstruct(result) || ~isfield(result, parts{i})
            ok = false;
            result = [];
            return;
        end
        result = result.(parts{i});
    end
end

function base = extract_base_name(field_path, fallback)
% 'TLP.Studentized.method7' -> 'TLP'
    if nargin < 2 || isempty(fallback)
        fallback = 'Method';
    end

    if isstring(field_path)
        if isscalar(field_path)
            field_path = char(field_path);
        else
            warning('Field path string array has %d elements; using the first one.', numel(field_path));
            field_path = char(field_path(1));
        end
    end

    if ischar(field_path) && ~isempty(field_path)
        parts = strsplit(field_path, '.');
        base = parts{1};
    else
        base = fallback;
    end
end

function name = get_method_name(method, fallback, iv_mode)
    if nargin < 3
        iv_mode = false;
    end

    name = fallback;
    if isempty(method) || ~isstruct(method)
        if isstring(method) && isscalar(method)
            name = char(method);
        elseif ischar(method)
            name = method;
        end
        name = normalize_method_display_name(name, iv_mode);
        return;
    end
    if isfield(method, 'B') && isfield(method.B, 'method_name')
        name = char(method.B.method_name);
    elseif isfield(method, 'method_name')
        name = char(method.method_name);
    elseif isfield(method, 'name')
        name = char(method.name);
    end
    name = normalize_method_display_name(name, iv_mode);
end

function color = get_method_color(method_name, base_colors, color_palette, palette_idx, all_names)
    base_types_list = {'TLP', 'LP', 'VAR', 'SLP', 'BLP'};
    this_base = infer_method_base(method_name);

    if nargin >= 5 && ~isempty(all_names)
        base_counts = zeros(1, length(base_types_list));
        for i = 1:length(base_types_list)
            base = base_types_list{i};
            for j = 1:length(all_names)
                if strcmp(infer_method_base(all_names{j}), base)
                    base_counts(i) = base_counts(i) + 1;
                end
            end
        end
        if any(base_counts > 1)
            color = color_palette(mod(palette_idx - 1, size(color_palette, 1)) + 1, :);
            return;
        end
    end

    if ismember(this_base, base_types_list)
        color = base_colors(this_base);
        return;
    end

    color = color_palette(mod(palette_idx - 1, size(color_palette, 1)) + 1, :);
end

function base = infer_method_base(method_name)
    name = lower(strtrim(char(method_name)));

    if startsWith(name, 'tlp')
        base = 'TLP'; return;
    end
    if startsWith(name, 'blp')
        base = 'BLP'; return;
    end
    if startsWith(name, 'slp')
        base = 'SLP'; return;
    end
    if contains(name, 'proxy var') || startsWith(name, 'var')
        base = 'VAR'; return;
    end
    if contains(name, 'lp iv') || startsWith(name, 'lp')
        base = 'LP'; return;
    end

    base = '';
end

function name = normalize_method_display_name(raw_name, iv_mode)
    raw = char(raw_name);
    normalized = strtrim(strrep(raw, '_', ' '));
    normalized = regexprep(normalized, '\s+', ' ');
    lower_name = lower(normalized);

    use_iv_label = iv_mode || contains(lower_name, ' iv') || contains(lower_name, 'proxy var');

    if startsWith(lower_name, 'lp iv') || strcmp(lower_name, 'lp') || startsWith(lower_name, 'lp ')
        if use_iv_label || startsWith(lower_name, 'lp iv')
            name = 'LP IV';
        else
            name = 'LP';
        end
        return;
    end

    if startsWith(lower_name, 'proxy var') || startsWith(lower_name, 'var')
        if use_iv_label || startsWith(lower_name, 'proxy var')
            name = 'Proxy VAR';
        else
            name = 'VAR';
        end
        return;
    end

    if startsWith(lower_name, 'tlp')
        name = 'TLP';
        return;
    end
    if startsWith(lower_name, 'blp')
        name = 'BLP';
        return;
    end
    if startsWith(lower_name, 'slp')
        name = 'SLP';
        return;
    end

    name = normalized;
end

function ls = get_linestyle_from_name(method_name)
% '--' for Ideal methods, '-' for everything else.
    if contains(lower(char(method_name)), 'ideal')
        ls = '--';
    else
        ls = '-';
    end
end

function normalized = normalize_field_specs(field_specs)
    normalized = {};

    if nargin < 1 || isempty(field_specs)
        return;
    end

    if ischar(field_specs) || (isstring(field_specs) && isscalar(field_specs))
        normalized = {field_specs};
    elseif isstruct(field_specs)
        error(['Invalid other_fields input type: struct. ', ...
            'For plot_stat_by_horizon_dual, pass method paths as text, e.g. ', ...
            '{''LP.Studentized.method7'',''VAR.Studentized.method7'',''BLP'',''SLP.Studentized.method7''}.']);
    elseif isstring(field_specs)
        normalized = cellstr(field_specs(:)).';
    elseif iscell(field_specs)
        for i = 1:numel(field_specs)
            item = field_specs{i};
            if isempty(item)
                continue;
            end

            if isstruct(item)
                error(['Invalid other_fields{%d}: struct provided. ', ...
                    'Use dot-path text for each method (e.g., ''LP.Studentized.method7'').'], i);
            end

            if isstring(item) && ~isscalar(item)
                expanded = cellstr(item(:));
                for j = 1:numel(expanded)
                    normalized{end+1} = expanded{j}; %#ok<AGROW>
                end
            else
                normalized{end+1} = item; %#ok<AGROW>
            end
        end
    else
        normalized = {field_specs};
    end

    for i = 1:numel(normalized)
        if isstring(normalized{i}) && isscalar(normalized{i})
            normalized{i} = char(normalized{i});
        end
    end
end

function highlighted = normalize_highlighted_field(highlighted_field)
    if nargin < 1 || isempty(highlighted_field)
        highlighted = '';
        return;
    end

    if iscell(highlighted_field)
        if isempty(highlighted_field)
            highlighted = '';
            return;
        end
        highlighted_field = highlighted_field{1};
    end

    if isstruct(highlighted_field)
        error(['Invalid highlighted_field input type: struct. ', ...
            'Use a dot-path text value like ''TLP.Studentized.method7'' or '''' to skip highlighting.']);
    end

    if isstring(highlighted_field)
        if isscalar(highlighted_field)
            highlighted = char(highlighted_field);
        else
            warning('highlighted_field string array has %d elements; using the first one.', numel(highlighted_field));
            highlighted = char(highlighted_field(1));
        end
        return;
    end

    highlighted = highlighted_field;
end

function tf = is_resolvable_field_spec(field_spec)
    tf = ~isempty(field_spec) && (ischar(field_spec) || (isstring(field_spec) && isscalar(field_spec)));
end

function [ref_val, ref_label, ref_style] = get_reference_line(metric)
    ref_val = [];
    ref_label = '';
    ref_style = '-';

    switch metric
        case 'coverage'
            ref_val = 0.9;
            ref_label = '90 %';
        case 'bias'
            ref_val = 0;
            ref_label = 'Zero';
        case {'below_ci', 'above_ci'}
            ref_val = 0.05;
            ref_label = '5 % (ideal)';
            ref_style = '--';
    end
end

function set_ylim(metric)
    switch metric
        case 'coverage'
            ylim([0.7 1]);
            yticks([0.5 0.6 0.7 0.8 0.9 1.0]);
        case {'below_ci', 'above_ci'}
            ylim([0 0.2]);
    end
end

function [horizons, values] = extract_metric_data(method, metric, p_picture, max_horizon)
    horizons = [];
    values = [];

    if nargin < 4
        max_horizon = [];
    end
    if isempty(method) || ~isstruct(method)
        return;
    end

    field_variants = {metric, upper(metric), lower(metric)};

    if strcmpi(metric, 'below_ci')
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
        if p_picture <= size(data_all, 3)
            data = squeeze(data_all(:, :, p_picture));
            if isrow(data)
                values = data;
            else
                values = data(:)';
            end
        end
    elseif ismatrix(data_all)
        [nrows, ncols] = size(data_all);
        if nrows == 1
            values = data_all;
        elseif ncols == 1
            values = data_all(:)';
        elseif ncols >= p_picture
            values = data_all(:, p_picture)';
        end
    end

    if isempty(values)
        return;
    end

    horizons = 0:(length(values) - 1);

    if length(values) == 1 && ~isempty(max_horizon) && max_horizon > 0
        values = repmat(values, 1, max_horizon + 1);
        horizons = 0:max_horizon;
    end
end

function format_axis_tick_labels(ax)
% Use compact numeric labels: 0 instead of 0.0, with trailing zeros trimmed.
    ticks = yticks(ax);
    labels = arrayfun(@compact_tick_label, ticks, 'UniformOutput', false);
    yticklabels(ax, labels);
end

function label = compact_tick_label(value)
    if abs(value) < 1e-12
        label = '0';
        return;
    end

    label = sprintf('%.6g', value);
    if strcmp(label, '-0')
        label = '0';
    end
end

function save_dual_figure(metrics, p_picture, str1, str2, save_picture)
    metric_str = upper(strjoin(metrics, '_'));

    T1 = str1.T;
    T2 = str2.T;

    if isfield(str1, 'OleaDGP') && str1.OleaDGP == 1
        dgp = 3;
    elseif isfield(str1, 'misspec_VARMA') && str1.misspec_VARMA == 1
        dgp = 2;
    else
        dgp = 1;
    end

    p_val = str1.P_VAR(p_picture(1));

    if ~exist('TablesAndPlots', 'dir')
        mkdir('TablesAndPlots');
    end

    % Build filename with optional prefix
    base_name = sprintf('dgp_%d_DUAL_%s_t%d_t%d_p%d', dgp, metric_str, T1, T2, p_val);
    if (ischar(save_picture) || isstring(save_picture))
        base_name = sprintf('%s_%s', char(save_picture), base_name);
    end
    
    filename = fullfile('TablesAndPlots', base_name);
    exportgraphics(gcf, [filename '.pdf'], 'ContentType', 'vector');
    fprintf('Saved: %s.pdf\n', filename);
end
