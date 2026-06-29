%% BLP (Ferreira, Miranda-Agrippino, Ricco 2023) Application for TLP
% Replicates the BLP empirical setup: 7-variable quarterly US VAR
% with Cholesky-identified FFR shock, 1965Q1–2019Q4
% Variables: RGDP, RCON, RINV, HOURS, WAGE, DEFL, FFR
clc
clear all
close all

rng(1)

% horizon as a row, LP/VAR/TLP estimate, significance with a star? 

% --- paths (robust to current working directory) ---
script_dir  = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, 'Subroutines')));
results_dir = fullfile(script_dir, 'SimulationsResults');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
addpath(results_dir);

%% ========================================================================
%  PART 1: LOAD DATA FROM BLP REPLICATION EXCEL FILE
%  ========================================================================
fprintf('=== Part 1: Loading Data from BLP Excel File ===\n');

% Set time parameters
tStart = 1954.50;  % First non-NaN observation in BLP data (~1954Q3)
tEnd   = 2019.75;  % Stop before COVID-19
periods_per_year = 4;
P_VAR = 4;  % BLP uses 5 lags

% Load directly from the BLP replication data file
[tempData, tempText] = xlsread(fullfile(script_dir, 'DataNew2019.xlsx'), 'FRED_Data_Levels');

% Extract dates and convert Excel serial dates to decimal years
excel_dates = tempData(:, 1);
matlab_dates = x2mdate(excel_dates);
dv = datevec(matlab_dates);
tVec_all = dv(:,1) + (dv(:,2)-1)/4;

% Variable labels from the spreadsheet
varLabels = tempText(2, 2:end);  % {'RGDP','RCON','RINV','HOUR','WAGE','DEFL','FFR'}
varData   = tempData(:, 2:end);  % Already in 100*log levels (FFR in levels)

fprintf('Variables: '); fprintf('%s ', varLabels{:}); fprintf('\n');
fprintf('Full data: %.2f to %.2f (%d obs)\n', tVec_all(1), tVec_all(end), size(varData,1));

% ---- Select estimation sample ----
est_idx = (tVec_all >= tStart) & (tVec_all <= tEnd);
Z0 = varData(est_idx, :);
tVec_final = tVec_all(est_idx);

% Clean any NaN rows
valid_rows = all(isfinite(Z0), 2);
Z0 = Z0(valid_rows, :);
tVec_final = tVec_final(valid_rows);

fprintf('Estimation sample: %.2f to %.2f\n', tVec_final(1), tVec_final(end));

% ---- Detrend ----
names_all = {'Real GDP', 'Real Consumption', 'Real Investment', ...
             'Hours', 'Real Wage', 'GDP Deflator', 'Fed Funds Rate'};
period_sine = 0;
plot_detrend = 0;
trend_option = 1;  % linear trend

[Z1, trend_VAR] = detrend_function(Z0, names_all, tVec_final, period_sine, plot_detrend, trend_option);

T = size(Z1, 1);
k = size(Z1, 2);
fprintf('Data processed. T = %d, k = %d\n', T, k);
fprintf('Sample: %.2f to %.2f\n', tVec_final(1), tVec_final(end));

%% ========================================================================
%  PART 2: TLP ESTIMATION — ONE RESPONSE AT A TIME TO SAVE MEMORY
%  ========================================================================
fprintf('\n=== Part 2: Estimating IRFs ===\n');

H_min = 0;
H_max = 24;   % Match BLP horizon
P_VAR = 4;    % Match BLP lags
P_LP  = 10;    % Match BLP LP lags
one_matrix = 0;
quadratic_method = 2;
HAC_kernel = 2;
method = 0;
use_studentized_boot = 0;
use_TLP_boot = 1;
use_cov_from_var = 0;
use_double_bootstrap = 1;
use_iid_boot_TLP = 0;
use_second_level_tlp_variance = 1;
alpha = 10;          % 90% bands to match BLP
bootstrapN = 200;
bootstrap_seed_vector_B1 = randi([0 (2^32 - 1)], bootstrapN, 1);
block_length = H_max + 1 - H_min;

% FFR is variable 7 — shock variable
% IRF specs: shock to FFR (var 7), responses in each of the 7 variables
irf_specs = {
    1, 7, 'RGDP_to_FFR';
    2, 7, 'RCON_to_FFR';
    3, 7, 'RINV_to_FFR';
    4, 7, 'HOURS_to_FFR';
    5, 7, 'WAGE_to_FFR';
    6, 7, 'DEFL_to_FFR';
    7, 7, 'FFR_to_FFR';
};

all_outputs = struct();

for spec_idx = 1:size(irf_specs, 1)
    
    which_irf_y = irf_specs{spec_idx, 1};
    which_irf_x = irf_specs{spec_idx, 2};
    spec_name   = irf_specs{spec_idx, 3};
    
    fprintf('Estimating: %s (%d/%d)\n', spec_name, spec_idx, size(irf_specs,1));
    
    iii = [[1:k^2]' repelem([1:k]',k,1) repmat([1:k]',k,1)];
    which_irf_var = iii(iii(:,2)==which_irf_y & iii(:,3)==which_irf_x);
    
    P_LP2 = P_LP + 1;
    
    str.H_min = H_min;
    str.H_max = H_max;
    str.P_VAR = P_VAR;
    str.P_VARq = P_VAR;
    str.P_LP = P_LP2;
    str.one_matrix = one_matrix;
    str.BootN = bootstrapN;
    str.BootN_TLP = bootstrapN;
    str.alpha = alpha;
    str.which_irf_y = which_irf_y;
    str.which_irf_x = which_irf_x;
    str.which_irf_var = which_irf_var;
    str.k = k;
    str.T = T;
    str.HAC_kernel = HAC_kernel;
    str.block_length = block_length;
    str.method = method;
    str.use_double_bootstrap = use_double_bootstrap;
    str.use_second_level_tlp_variance = use_second_level_tlp_variance;
    str.bootstrap_seed_vector_B1 = bootstrap_seed_vector_B1;
    
    data = struct;
    data.y_t = Z1;
    
    % Run standard estimators
    temp2 = irf_SLP_function(data, str, quadratic_method, 1, 1, use_studentized_boot, 1);
    irf_LP  = temp2.LP;
    irf_SLP = temp2.SLP;
    clear temp2;  % free memory
    
    irf_VAR  = cell(1, 1);
    irf_BLP  = cell(1, 1);
    irf_TLP2 = cell(1, 1);
    
    str2 = str;
    str2.use_double_bootstrap = use_double_bootstrap;
    
    irf_VAR{1}  = irf_VAR_function(data, str2, str2.P_VAR);
    irf_BLP{1}  = irf_BLP_function(data, str2);
    compute_SLP = 1;
    irf_TLP2{1} = irf_TLP_function(data, str2, use_TLP_boot, use_iid_boot_TLP, use_cov_from_var, [], compute_SLP);
    
    output_s = struct();
    output_s.irf_LP   = irf_LP;
    output_s.irf_BLP  = irf_BLP;
    output_s.irf_SLP  = irf_SLP;
    output_s.irf_VAR  = irf_VAR;
    output_s.irf_TLP2 = irf_TLP2;
    
    % Free memory before unpack
    clear irf_LP irf_SLP irf_VAR irf_BLP irf_TLP2;
    
    output_cell = cell(1,1);
    output_cell{1} = output_s;
    clear output_s;  % free the struct
    
    output = unpack_cell(output_cell, str);
    clear output_cell;  % kill the cell immediately to free memory
    
    all_outputs.(spec_name) = output;
    all_outputs.([spec_name '_str']) = str;
    clear output;  % free this copy too
    
    fprintf('  Done. Memory freed.\n');
end
fprintf('All estimations complete.\n');

%% ========================================================================
%  PART 3: PLOTTING DATA + SELECTED FIGURE
%  ========================================================================
fprintf('\n=== Part 3: Preparing BLP application plots ===\n');

plot_blp_4x7 = 0;
plot_blp_overlay_4x2 = 1;
plot_blp_weights_4x2 = 0;

p_picture = 1;
H_truncate = H_max + 1;
h_start = H_min;
h_end   = min(H_max, H_min + H_truncate - 1);
x       = h_start:h_end;
n_h     = length(x);

k_models = 4;
model_names = {'TLP', 'VAR', 'LP', 'BLP'};

colors = [
    0.9  0.5  0.1;   % TLP     — orange
    0.8  0.4  0.8;   % VAR     — purple
    0.1  0.6  0.8;   % LP      — blue
    0.2  0.6  0.4;   % BLP     — green
];
linestyles = {'-', '-.', '-', '--'};
lineWidth  = 2.0;

n_resp = size(irf_specs, 1);
resp_names = irf_specs(:, 3)';
col_titles = names_all;

% Extract IRFs and CIs
IRFs = cell(k_models, n_resp);
CIs  = cell(k_models, n_resp);

for r = 1:n_resp
    output = all_outputs.(resp_names{r});
    
    % 1. TLP
    d = output.TLP2.Studentized.method7;
    IRFs{1,r} = d.irf(1:n_h, 1, p_picture);
    CIs{1,r}  = squeeze(d.CI(1:n_h, :, 1, p_picture));
    
    % 2. VAR
    d = output.VAR2.Studentized.method7;
    IRFs{2,r} = d.irf(1:n_h, 1, p_picture);
    CIs{2,r}  = squeeze(d.CI(1:n_h, :, 1, p_picture));
    
    % 3. LP
    d = output.LP2.Studentized.method7;
    IRFs{3,r} = d.irf(1:n_h, 1, p_picture);
    CIs{3,r}  = squeeze(d.CI(1:n_h, :, 1, p_picture));
    
    % 4. BLP
    d = output.BLP;
    IRFs{4,r} = d.irf(1:n_h, 1, p_picture);
    CIs{4,r}  = squeeze(d.CI(1:n_h, :, 1, p_picture));
end

% Column-specific y-limits
ylims_cols = cell(1, n_resp);
for r = 1:n_resp
    ylo = inf; yhi = -inf;
    for m = 1:k_models
        vals = [IRFs{m,r}; CIs{m,r}(:)];
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            ylo = min(ylo, min(vals));
            yhi = max(yhi, max(vals));
        end
    end
    pad = 0.08 * (yhi - ylo);
    if pad == 0, pad = 0.1; end
    ylims_cols{r} = [ylo - pad, yhi + pad];
end

if plot_blp_4x7
    fprintf('\n=== Part 3a: Plotting 4x7 BLP-style Figure ===\n');
    fig = figure('Units', 'Normalized', 'OuterPosition', [0.02 0.02 0.96 0.92]);
    tl  = tiledlayout(k_models, n_resp, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax_font    = 11;
    label_font = 12;
    title_font = 14;

    for m = 1:k_models
        for r = 1:n_resp
            nexttile; hold on;

            ci = CIs{m,r};
            valid = isfinite(ci(:,1));
            xv = x(valid);
            fill([xv fliplr(xv)], ...
                 [ci(valid,1)' fliplr(ci(valid,2)')], ...
                 colors(m,:), 'FaceAlpha', 0.20, 'EdgeColor', 'none');

            plot(x, IRFs{m,r}, 'Color', colors(m,:), ...
                 'LineWidth', lineWidth, 'LineStyle', linestyles{m});

            yline(0, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
            xlim([h_start h_end]);
            ylim(ylims_cols{r});

            set(gca, 'FontSize', ax_font);
            grid on; box on;

            base_ticks = 0:5:h_end;
            if base_ticks(end) ~= h_end, base_ticks = [base_ticks h_end]; end
            xticks(base_ticks);

            if m == 1
                title(col_titles{r}, 'FontSize', title_font, 'FontWeight', 'bold');
            end

            if r == 1
                ylabel(model_names{m}, 'FontSize', label_font, 'FontWeight', 'bold');
            end

            if m == k_models
                xlabel('Quarters', 'FontSize', label_font);
            end
        end
    end
    set(gcf, 'Color', 'w');
    sgtitle(sprintf('Response to 100 bp FFR Shock — BLP Application (%.0f–%.0f)', floor(tStart), floor(tEnd)), ...
        'FontSize', 18, 'FontWeight', 'bold');

    if ~exist('TablesAndPlots', 'dir'), mkdir('TablesAndPlots'); end
    fname = fullfile('TablesAndPlots', sprintf('BLP_application_4x%d_pL%d_pV%d', n_resp, P_LP, P_VAR));
    exportgraphics(gcf, [fname '.pdf'], 'ContentType', 'vector');
    fprintf('Saved: %s.pdf\n', fname);
end

%% ========================================================================
%  PART 3b: OVERLAY PLOT — 4×2 (7 responses + legend tile)
%  ========================================================================
if plot_blp_overlay_4x2
fprintf('\n=== Part 3b: Plotting 4x2 Overlay Figure ===\n');

n_rows_ov = 4;
n_cols_ov = 2;

colors_ov = [
    0.9  0.5  0.1;   % TLP     — orange
    0.8  0.4  0.8;   % VAR     — purple
    0.1  0.6  0.8;   % LP      — blue
];
ls_ov = {'-', '-.', '-'};
lw_ov = [3.0, 2.0, 2.0];
weight_color = [0.15 0.15 0.15];
weight_ls = '--';
weight_lw = 2.4;
model_names_ov = {'Targeted Local Projection', 'Vector Autoregression', 'Local Projection', 'Weight (LP share)'};
k_ov = 3;  % TLP=1, VAR=2, LP=3 (indices into IRFs/CIs from Part 3)

fig2 = figure('Units', 'Normalized', 'OuterPosition', [0.05 0.025 0.42 0.88]);
tl2  = tiledlayout(n_rows_ov, n_cols_ov, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- INCREASED FONT SIZES ---
ax_font2    = 16;
label_font2 = 18;
title_font2 = 20;

% Compute y-limits across TLP/VAR/LP per response
ylims_ov = cell(1, n_resp);
for r = 1:n_resp
    ylo = inf; yhi = -inf;
    for m = 1:k_ov
        vals = [IRFs{m,r}; CIs{m,r}(:)];
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            ylo = min(ylo, min(vals));
            yhi = max(yhi, max(vals));
        end
    end
    pad = 0.10 * (yhi - ylo);
    if pad == 0, pad = 0.1; end
    ylims_ov{r} = [ylo - pad, yhi + pad];
end

% Extract TLP LP-share weights (v_lambda) for dual-axis plotting
W_LP = nan(n_h, n_resp);
for r = 1:n_resp
    output = all_outputs.(resp_names{r});
    if isfield(output, 'TLP2') && isfield(output.TLP2, 'Studentized') ...
            && isfield(output.TLP2.Studentized, 'method7') ...
            && isfield(output.TLP2.Studentized.method7, 'v_lambda')
        v_tmp = squeeze(output.TLP2.Studentized.method7.v_lambda(1:n_h, 1, p_picture));
        W_LP(:, r) = v_tmp(:);
    end
end
W_LP = max(0, min(1, W_LP));

h_lines = gobjects(k_ov, 1);

% Plot order: bands back-to-front, then lines back-to-front
% so TLP (plotted last) is on top
plot_order = [3, 2, 1];  % LP, VAR, TLP
for r = 1:n_resp
    nexttile; hold on;
    yyaxis left;
    
    % CIs: plot in reverse order so TLP band is on top
    alpha_vals = [0.15, 0.15, 0.25];  % LP, VAR, TLP gets more opaque
    for idx = 1:k_ov
        m = plot_order(idx);
        ci = CIs{m,r};
        valid = isfinite(ci(:,1));
        xv = x(valid);
        fill([xv fliplr(xv)], ...
             [ci(valid,1)' fliplr(ci(valid,2)')], ...
             colors_ov(m,:), 'FaceAlpha', alpha_vals(idx), 'EdgeColor', colors_ov(m,:), ...
             'EdgeAlpha', 0.3, 'LineWidth', 0.5);
    end
    
    % Lines: plot in reverse order so TLP is on top
    for idx = 1:k_ov
        m = plot_order(idx);
        h = plot(x, IRFs{m,r}, 'Color', colors_ov(m,:), ...
             'LineWidth', lw_ov(m), 'LineStyle', ls_ov{m});
        if r == 1
            h_lines(m) = h;
        end
    end
    
    yline(0, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
    xlim([h_start h_end]);
    ylim(ylims_ov{r});
    
    yyaxis right;
    w_lp = W_LP(:, r);
    valid_w = isfinite(w_lp);
    plot(x(valid_w), w_lp(valid_w), 'Color', weight_color, ...
        'LineWidth', weight_lw, 'LineStyle', weight_ls);
    % If weights hit corners, draw near-corner dashed guides.
    hit_one = any((w_lp >= 1) & valid_w);
    hit_zero = any((w_lp <= 0) & valid_w);
    if hit_one
        yline(0.999, '--', 'Color', [0.35 0.35 0.35], ...
            'LineWidth', 1.0, 'HandleVisibility', 'off');
    end
    if hit_zero
        yline(0.001, '--', 'Color', [0.35 0.35 0.35], ...
            'LineWidth', 1.0, 'HandleVisibility', 'off');
    end
    ylim([0 1]);
    set(gca, 'YColor', [0.2 0.2 0.2]);

    set(gca, 'FontSize', ax_font2);
    grid on; box on;
    
    base_ticks = 0:5:h_end;
    if base_ticks(end) ~= h_end, base_ticks = [base_ticks h_end]; end
    xticks(base_ticks);
    
    title(col_titles{r}, 'FontSize', title_font2, 'FontWeight', 'bold');
    
    if mod(r-1, n_cols_ov) == 0
        yyaxis left;
        ylabel('Percentage points', 'FontSize', label_font2);
    end
    if mod(r-1, n_cols_ov) == 1
        yyaxis right;
        ylabel('Weight', 'FontSize', label_font2);
    end
    
    if r > n_resp - n_cols_ov
        xlabel('Quarters', 'FontSize', label_font2);
    end
end

% --- FIXED LEGEND TILE ---
% Last tile: legend only
nexttile;
axis off;

% Dummy plots for clean legend handles
hold on;
h_leg = gobjects(k_ov + 1, 1);
for m = 1:k_ov
    h_leg(m) = plot(NaN, NaN, 'Color', colors_ov(m,:), ...
        'LineWidth', lw_ov(m), 'LineStyle', ls_ov{m});
end
h_leg(k_ov + 1) = plot(NaN, NaN, 'Color', weight_color, ...
    'LineWidth', weight_lw, 'LineStyle', weight_ls);

% Let MATLAB handle the centering automatically to prevent spillover
lgd = legend(h_leg, model_names_ov, ...
    'FontSize', label_font2, ...
    'Box', 'off', ...
    'Location', 'east');

% Force consistent axis alignment across all tiles
drawnow;  % let MATLAB compute positions
all_ax = flipud(findobj(fig2, 'Type', 'axes'));
plot_axes = [];
for i = 1:length(all_ax)
    if ~isempty(get(all_ax(i), 'Children')) || strcmp(get(all_ax(i), 'Visible'), 'on')
        % skip the legend-only tile (last one, axis off)
        if strcmp(get(all_ax(i), 'Visible'), 'off')
            continue;
        end
        plot_axes = [plot_axes; all_ax(i)];
    end
end

% Separate left and right columns
left_ax = []; right_ax = [];
for i = 1:length(plot_axes)
    pos = plot_axes(i).Position;
    if pos(1) < 0.5
        left_ax = [left_ax; plot_axes(i)];
    else
        right_ax = [right_ax; plot_axes(i)];
    end
end

% Align left column
if ~isempty(left_ax)
    max_x = max(arrayfun(@(a) a.Position(1), left_ax));
    min_w = min(arrayfun(@(a) a.Position(3), left_ax));
    for i = 1:length(left_ax)
        p = left_ax(i).Position;
        left_ax(i).Position = [max_x, p(2), min_w, p(4)];
    end
end

% Align right column
if ~isempty(right_ax)
    min_x = min(arrayfun(@(a) a.Position(1), right_ax));
    min_w = min(arrayfun(@(a) a.Position(3), right_ax));
    for i = 1:length(right_ax)
        p = right_ax(i).Position;
        right_ax(i).Position = [min_x, p(2), min_w, p(4)];
    end
end

set(gcf, 'Color', 'w');

% --- INCREASED SUPER TITLE SIZE ---
sgtitle(sprintf('Response to 100 bp FFR Shock (%.0f–%.0f)', floor(tStart), floor(tEnd)), ...
    'FontSize', 24, 'FontWeight', 'bold');

fname2 = fullfile('TablesAndPlots', sprintf('BLP_application_overlay_4x2_pL%d_pV%d', P_LP, P_VAR));
if ~exist('TablesAndPlots', 'dir'), mkdir('TablesAndPlots'); end
exportgraphics(gcf, [fname2 '.pdf'], 'ContentType', 'vector');
fprintf('Saved: %s.pdf\n', fname2);
end

%% ========================================================================
%  PART 3c: TLP WEIGHTS PLOT - 4x2 (7 responses + legend tile)
%  ========================================================================
if plot_blp_weights_4x2
fprintf('\n=== Part 3c: Plotting TLP weights (4x2) ===\n');

n_rows_w = 4;
n_cols_w = 2;

color_w = [0.1  0.6  0.8];
ls_w = '-';
lw_w = 3.0;
weight_name = 'Weight (LP share)';

W_LP = nan(n_h, n_resp);
for r = 1:n_resp
    output = all_outputs.(resp_names{r});
    if isfield(output, 'TLP2') && isfield(output.TLP2, 'Studentized') ...
            && isfield(output.TLP2.Studentized, 'method7') ...
            && isfield(output.TLP2.Studentized.method7, 'v_lambda')
        v_tmp = squeeze(output.TLP2.Studentized.method7.v_lambda(1:n_h, 1, p_picture));
        W_LP(:, r) = v_tmp(:);
    end
end
W_LP = max(0, min(1, W_LP));

figw = figure('Units', 'Normalized', 'OuterPosition', [0.05 0.025 0.42 0.88]);
tlw  = tiledlayout(n_rows_w, n_cols_w, 'Padding', 'compact', 'TileSpacing', 'compact');

ax_font_w    = 16;
label_font_w = 18;
title_font_w = 20;

for r = 1:n_resp
    nexttile; hold on;
    
    w_lp = W_LP(:, r);
    valid = isfinite(w_lp);
    xv = x(valid);
    
    plot(xv, w_lp(valid), 'Color', color_w, ...
        'LineWidth', lw_w, 'LineStyle', ls_w);
    % If weights hit corners, draw near-corner dashed guides.
    hit_one = any((w_lp >= 1) & valid);
    hit_zero = any((w_lp <= 0) & valid);
    if hit_one
        yline(0.999, '--', 'Color', [0.35 0.35 0.35], ...
            'LineWidth', 1.0, 'HandleVisibility', 'off');
    end
    if hit_zero
        yline(0.001, '--', 'Color', [0.35 0.35 0.35], ...
            'LineWidth', 1.0, 'HandleVisibility', 'off');
    end
    
    yline(0.5, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
    xlim([h_start h_end]);
    ylim([0 1]);
    
    set(gca, 'FontSize', ax_font_w);
    grid on; box on;
    
    base_ticks = 0:5:h_end;
    if base_ticks(end) ~= h_end, base_ticks = [base_ticks h_end]; end
    xticks(base_ticks);
    
    title(col_titles{r}, 'FontSize', title_font_w, 'FontWeight', 'bold');
    
    if mod(r-1, n_cols_w) == 0
        ylabel('Weight', 'FontSize', label_font_w);
    end
    
    if r > n_resp - n_cols_w
        xlabel('Quarters', 'FontSize', label_font_w);
    end
end

% Legend tile
nexttile;
axis off;
hold on;
h_leg_w = plot(NaN, NaN, 'Color', color_w, 'LineWidth', lw_w, 'LineStyle', ls_w);
legend(h_leg_w, weight_name, 'FontSize', label_font_w, 'Box', 'off', 'Location', 'east');

set(gcf, 'Color', 'w');
sgtitle(sprintf('TLP weights after 100 bp FFR shock (%.0f-%.0f)', floor(tStart), floor(tEnd)), ...
    'FontSize', 24, 'FontWeight', 'bold');

fnamew = fullfile('TablesAndPlots', sprintf('BLP_application_weights_4x2_pL%d_pV%d', P_LP, P_VAR));
if ~exist('TablesAndPlots', 'dir'), mkdir('TablesAndPlots'); end
exportgraphics(gcf, [fnamew '.pdf'], 'ContentType', 'vector');
fprintf('Saved: %s.pdf\n', fnamew);
end

%% ========================================================================
%  PART 4: SAVE WORKSPACE
%  ========================================================================
fprintf('\n=== Saving Workspace ===\n');

mat_filename = fullfile(results_dir, sprintf('BLP_application_T%d_pL%d_pV%d_FFRshock.mat', T, P_LP, P_VAR));
clear fig fig2 figw tl tl2 tlw lgd h_leg h_leg_w h_lines plot_axes left_ax right_ax all_ax;
save(mat_filename);
fprintf('Workspace saved to: %s\n', mat_filename);


fprintf('Done.\n');


%% ========================================================================
%  PART 5: SIGNIFICANCE TABLE FOR FIRST 4 VARIABLES
%  ========================================================================
fprintf('\n=== Part 5: Generating Significance Table ===\n');

% Variables to include (first 4 from your irf_specs)
table_vars = {'RGDP_to_FFR', 'RCON_to_FFR', 'RINV_to_FFR', 'HOURS_to_FFR'};
var_labels = {'Real GDP', 'Real Consumption', 'Real Investment', 'Hours Worked'};
model_keys = {'LP2', 'TLP2', 'VAR2'}; % Keys in your output struct
model_names_tab = {'LP', 'TLP', 'VAR'};

n_table_vars = length(table_vars);
n_models = length(model_keys);
H_table = min(16, H_max);

% Initialize cell array for the table
% Rows: Header1 (Variables), Header2 (Models), Data (0 to 16)
% Cols: Horizon + (4 variables * 3 models) = 13
table_cell = cell(H_table + 3, 1 + n_table_vars * n_models);

% Fill Headers
table_cell{1, 1} = 'Horizon';
table_cell{2, 1} = 'h';

col_idx = 2;
for v = 1:n_table_vars
    % Center the variable name over its 3 columns in the cell array concept
    table_cell{1, col_idx} = var_labels{v}; 
    for m = 1:n_models
        table_cell{2, col_idx} = model_names_tab{m};
        col_idx = col_idx + 1;
    end
end

% Fill Data
p_picture = 1; 
for h_idx = 1:(H_table + 1)
    h = h_idx - 1; % Horizon 0 to 16
    row_idx = h_idx + 2;
    table_cell{row_idx, 1} = sprintf('%d', h); % Horizon column
    
    col_idx = 2;
    for v = 1:n_table_vars
        output = all_outputs.(table_vars{v});
        for m = 1:n_models
            % Extract data struct
            d = output.(model_keys{m}).Studentized.method7;
            
            % Get estimate and confidence bounds at horizon h_idx
            est = d.irf(h_idx, 1, p_picture);
            ci_low = d.CI(h_idx, 1, 1, p_picture);
            ci_high = d.CI(h_idx, 2, 1, p_picture);
            
            % Check significance: if 0 is NOT between ci_low and ci_high
            is_sig = (ci_low > 0) || (ci_high < 0);
            
            % Format string with or without star
            if is_sig
                table_cell{row_idx, col_idx} = sprintf('%.2f*', est);
            else
                table_cell{row_idx, col_idx} = sprintf('%.2f', est);
            end
            col_idx = col_idx + 1;
        end
    end
end

% --- 1. Display table in MATLAB Command Window ---
disp(' ');
disp('Significance Table (90% CI based on alpha=10):');
for r = 1:size(table_cell, 1)
    row_str = '';
    for c = 1:size(table_cell, 2)
        if isempty(table_cell{r,c})
            val = '';
        else
            val = table_cell{r,c};
        end
        % Pad strings to keep columns aligned in the console
        row_str = [row_str, sprintf('%-12s', val)];
    end
    disp(row_str);
end

% --- 2. Export table directly to a LaTeX file ---
if ~exist('TablesAndPlots', 'dir'), mkdir('TablesAndPlots'); end
latex_filename = fullfile('TablesAndPlots', 'IRF_Significance_Table.tex');
fid = fopen(latex_filename, 'w');

fprintf(fid, '\\begin{table}[H]\n\\centering\n');
fprintf(fid, '\\resizebox{\\textwidth}{!}{\n');
fprintf(fid, '\\begin{tabular}{l *{12}{c}}\n');
fprintf(fid, '\\hline\\hline\n');

% Top Header (Variable names spanning 3 columns each)
fprintf(fid, 'Horizon');
for v = 1:n_table_vars
    fprintf(fid, ' & \\multicolumn{3}{c}{%s}', var_labels{v});
end
fprintf(fid, ' \\\\\n');

% Second Header (Model names)
fprintf(fid, '$h$');
for v = 1:n_table_vars
    for m = 1:n_models
        fprintf(fid, ' & %s', model_names_tab{m});
    end
end
fprintf(fid, ' \\\\\n\\hline\n');

% Data Rows
for r = 3:size(table_cell, 1)
    fprintf(fid, '%s', table_cell{r, 1}); % Horizon
    for c = 2:size(table_cell, 2)
        val = table_cell{r, c};
        model_idx = mod(c - 2, n_models) + 1;
        variable_idx = floor((c - 2) / n_models) + 1;
        lp_col = 2 + (variable_idx - 1) * n_models;
        lp_val = table_cell{r, lp_col};
        highlight_tlp = strcmp(model_names_tab{model_idx}, 'TLP') ...
            && contains(val, '*') && ~contains(lp_val, '*');
        if contains(val, '*') && highlight_tlp
            numeric_val = strrep(val, '*', '');
            val = sprintf('$\\textcolor{blue}{\\boldsymbol{%s^{*}}}$', numeric_val);
        elseif contains(val, '*')
            numeric_val = strrep(val, '*', '');
            val = sprintf('$%s^{*}$', numeric_val);
        else
            val = sprintf('$%s$', val);
        end
        fprintf(fid, ' & %s', val);
    end
    fprintf(fid, ' \\\\\n');
end
fprintf(fid, '\\hline\\hline\n\\end{tabular}}\n');
fprintf(fid, '\\caption{Impulse Responses to 100 bp Federal Funds Rate Shock (bold blue indicates TLP(8,4) significance where LP(8) is insignificant; a star indicates 90\\%% significance). The table is truncated to include first 16 quarters.}\n');
fprintf(fid, '\\label{tab:irf_sig}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);

fprintf('\nLaTeX table saved to: %s\n', latex_filename);
