function table_out = BLP_application_significance_table(mat_file, table_vars, horizons, output_dir, latex_filename)
%BLP_APPLICATION_SIGNIFICANCE_TABLE Build BLP application significance table.
%
%   table_out = BLP_application_significance_table(mat_file)
%   loads the saved BLP application workspace, prints the default
%   significance table, saves the LaTeX file, returns the table, and assigns
%   BLP_significance_table in the base workspace.
%
%   table_vars can be a cell array of response field names, numeric indices
%   into irf_specs, or [] for the default first four responses.
%
%   horizons can be a scalar H, interpreted as 0:H, or a vector of explicit
%   horizons such as [0 4 8 12 16 20 24].

    if nargin < 1 || isempty(mat_file)
        mat_file = 'BLP_application_T259_pL10_pV4_FFRshock.mat';
    end
    if nargin < 2
        table_vars = [];
    end
    if nargin < 3 || isempty(horizons)
        horizons = [];
    end
    if nargin < 4 || isempty(output_dir)
        output_dir = 'TablesAndPlots';
    end
    if nargin < 5 || isempty(latex_filename)
        latex_filename = 'IRF_Significance_Table.tex';
    end

    required_fields = {'all_outputs', 'irf_specs', 'H_max'};
    S = load(mat_file, required_fields{:});
    for f = 1:numel(required_fields)
        if ~isfield(S, required_fields{f})
            error('BLP_application_significance_table:missingField', ...
                'Missing required field "%s" in %s.', required_fields{f}, mat_file);
        end
    end

    all_outputs = S.all_outputs;
    irf_specs = S.irf_specs;
    H_max = S.H_max;

    [table_vars, var_labels] = resolve_table_variables(table_vars, irf_specs);
    horizons = resolve_horizons(horizons, H_max);

    model_keys = {'LP2', 'TLP2', 'VAR2'};
    model_names_tab = {'LP', 'TLP', 'VAR'};

    n_table_vars = length(table_vars);
    n_models = length(model_keys);

    table_cell = cell(length(horizons) + 2, 1 + n_table_vars * n_models);
    table_cell{1, 1} = 'Horizon';
    table_cell{2, 1} = 'h';

    col_idx = 2;
    for v = 1:n_table_vars
        table_cell{1, col_idx} = var_labels{v};
        for m = 1:n_models
            table_cell{2, col_idx} = model_names_tab{m};
            col_idx = col_idx + 1;
        end
    end

    p_picture = 1;
    for h_pos = 1:length(horizons)
        h = horizons(h_pos);
        h_idx = h + 1;
        row_idx = h_pos + 2;
        table_cell{row_idx, 1} = sprintf('%d', h);

        col_idx = 2;
        for v = 1:n_table_vars
            output = all_outputs.(table_vars{v});
            for m = 1:n_models
                d = output.(model_keys{m}).Studentized.method7;

                est = d.irf(h_idx, 1, p_picture);
                ci_low = d.CI(h_idx, 1, 1, p_picture);
                ci_high = d.CI(h_idx, 2, 1, p_picture);

                is_sig = (ci_low > 0) || (ci_high < 0);
                if is_sig
                    table_cell{row_idx, col_idx} = sprintf('%.2f*', est);
                else
                    table_cell{row_idx, col_idx} = sprintf('%.2f', est);
                end
                col_idx = col_idx + 1;
            end
        end
    end

    disp(' ');
    disp('Significance Table (90% CI based on alpha=10):');
    print_table_cell(table_cell);

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    latex_path = fullfile(output_dir, latex_filename);
    write_latex_table(latex_path, table_cell, var_labels, model_names_tab);
    fprintf('\nLaTeX table saved to: %s\n', latex_path);

    display_table = cell2table(table_cell(3:end, :), ...
        'VariableNames', matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(table_cell(2, :))));

    table_out = struct();
    table_out.cell = table_cell;
    table_out.display_table = display_table;
    table_out.latex_file = latex_path;
    table_out.variables = table_vars;
    table_out.variable_labels = var_labels;
    table_out.horizons = horizons;

    assignin('base', 'BLP_significance_table', table_out);
end

function [table_vars, var_labels] = resolve_table_variables(table_vars, irf_specs)
    all_vars = irf_specs(:, 3)';
    default_labels = {'Real GDP', 'Real Consumption', 'Real Investment', ...
        'Hours Worked', 'Real Wage', 'GDP Deflator', 'Fed Funds Rate'};

    if isempty(table_vars)
        idx = 1:min(4, numel(all_vars));
        table_vars = all_vars(idx);
        var_labels = default_labels(idx);
        return;
    end

    if isnumeric(table_vars)
        idx = table_vars(:)';
        if any(idx < 1) || any(idx > numel(all_vars))
            error('BLP_application_significance_table:badVariableIndex', ...
                'Variable indices must be between 1 and %d.', numel(all_vars));
        end
        table_vars = all_vars(idx);
        var_labels = default_labels(idx);
        return;
    end

    if ischar(table_vars) || isstring(table_vars)
        table_vars = cellstr(table_vars);
    end

    var_labels = cell(size(table_vars));
    for v = 1:numel(table_vars)
        match_idx = find(strcmp(all_vars, table_vars{v}), 1);
        if isempty(match_idx)
            error('BLP_application_significance_table:badVariableName', ...
                'Unknown response "%s".', table_vars{v});
        end
        var_labels{v} = default_labels{match_idx};
    end
end

function horizons = resolve_horizons(horizons, H_max)
    if isempty(horizons)
        horizons = 0:H_max;
    elseif isscalar(horizons)
        horizons = 0:horizons;
    else
        horizons = horizons(:)';
    end

    if any(horizons < 0) || any(horizons > H_max) || any(fix(horizons) ~= horizons)
        error('BLP_application_significance_table:badHorizons', ...
            'Horizons must be integer values between 0 and %d.', H_max);
    end
end

function print_table_cell(table_cell)
    for r = 1:size(table_cell, 1)
        row_str = '';
        for c = 1:size(table_cell, 2)
            if isempty(table_cell{r,c})
                val = '';
            else
                val = table_cell{r,c};
            end
            row_str = [row_str, sprintf('%-12s', val)];
        end
        disp(row_str);
    end
end

function write_latex_table(latex_path, table_cell, var_labels, model_names_tab)
    n_table_vars = numel(var_labels);
    n_models = numel(model_names_tab);

    fid = fopen(latex_path, 'w');
    if fid < 0
        error('BLP_application_significance_table:cannotOpenFile', ...
            'Could not open %s for writing.', latex_path);
    end

    fprintf(fid, '\\begin{table}[H]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\resizebox{\\textwidth}{!}{\n');
    fprintf(fid, '\\begin{tabular}{l *{%d}{c}}\n', n_table_vars * n_models);
    fprintf(fid, '\\hline\\hline\n');

    fprintf(fid, 'Horizon');
    for v = 1:n_table_vars
        fprintf(fid, ' & \\multicolumn{%d}{c}{%s}', n_models, var_labels{v});
    end
    fprintf(fid, ' \\\\\n');

    fprintf(fid, '$h$');
    for v = 1:n_table_vars
        for m = 1:n_models
            fprintf(fid, ' & %s', model_names_tab{m});
        end
    end
    fprintf(fid, ' \\\\\n\\hline\n');

    for r = 3:size(table_cell, 1)
        fprintf(fid, '%s', table_cell{r, 1});
        for c = 2:size(table_cell, 2)
            val = table_cell{r, c};
            model_idx = mod(c - 2, n_models) + 1;
            variable_idx = floor((c - 2) / n_models) + 1;
            lp_col = 2 + (variable_idx - 1) * n_models;
            lp_val = table_cell{r, lp_col};
            highlight_tlp = strcmp(model_names_tab{model_idx}, 'TLP') ...
                && contains(val, '*') && ~contains(lp_val, '*');
            fprintf(fid, ' & %s', latex_table_value(val, highlight_tlp));
        end
        fprintf(fid, ' \\\\\n');
    end

    fprintf(fid, '\\hline\\hline\n');
    fprintf(fid, '\\end{tabular}}\n');
    fprintf(fid, '\\caption{Impulse Responses to 100 bp Federal Funds Rate Shock (bold blue indicates TLP(8,4) significance where LP(8) is insignificant; a star indicates 90\\%% significance). The table is truncated to include first 16 quarters.}\n');
    fprintf(fid, '\\label{tab:irf_sig}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);
end

function latex_value = latex_table_value(raw_value, highlight)
    if contains(raw_value, '*')
        numeric_value = strrep(raw_value, '*', '');
        if highlight
            latex_value = sprintf('$\\textcolor{blue}{\\boldsymbol{%s^{*}}}$', numeric_value);
        else
            latex_value = sprintf('$%s^{*}$', numeric_value);
        end
    else
        latex_value = sprintf('$%s$', raw_value);
    end
end
