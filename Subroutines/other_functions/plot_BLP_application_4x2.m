function fig_handle = plot_BLP_application_4x2(mat_file, save_picture, close_picture, output_dir)
%PLOT_BLP_APPLICATION_4X2 Plot BLP application IRFs from saved results.
%
%   plot_BLP_application_4x2(mat_file, save_picture, close_picture)
%   loads the saved ApplicationQuarterly_BLP workspace and recreates the
%   4x2 overlay impulse-response figure without rerunning estimation.

    if nargin < 1 || isempty(mat_file)
        mat_file = 'BLP_application_T259_pL10_pV4_FFRshock.mat';
    end
    if nargin < 2 || isempty(save_picture)
        save_picture = 1;
    end
    if nargin < 3 || isempty(close_picture)
        close_picture = 0;
    end
    if nargin < 4 || isempty(output_dir)
        output_dir = 'TablesAndPlots';
    end

    required_fields = {'all_outputs', 'irf_specs', 'names_all', 'H_min', ...
        'H_max', 'P_LP', 'P_VAR', 'tStart', 'tEnd'};
    S = load(mat_file, required_fields{:});
    for f = 1:numel(required_fields)
        if ~isfield(S, required_fields{f})
            error('plot_BLP_application_4x2:missingField', ...
                'Missing required field "%s" in %s.', required_fields{f}, mat_file);
        end
    end

    all_outputs = S.all_outputs;
    irf_specs = S.irf_specs;
    names_all = S.names_all;
    H_min = S.H_min;
    H_max = S.H_max;
    P_LP = S.P_LP;
    P_VAR = S.P_VAR;
    tStart = S.tStart;
    tEnd = S.tEnd;

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    p_picture = 1;
    H_truncate = H_max + 1;
    h_start = H_min;
    h_end = min(H_max, H_min + H_truncate - 1);
    x = h_start:h_end;
    n_h = length(x);

    n_resp = size(irf_specs, 1);
    resp_names = irf_specs(:, 3)';
    col_titles = names_all;

    [IRFs, CIs] = collect_blp_overlay_data(all_outputs, resp_names, n_h, p_picture);

    n_rows_ov = 4;
    n_cols_ov = 2;

    colors_ov = [
        0.9  0.5  0.1;   % TLP
        0.8  0.4  0.8;   % VAR
        0.1  0.6  0.8;   % LP
    ];
    ls_ov = {'-', '-.', '-'};
    lw_ov = [3.0, 2.0, 2.0];
    weight_color = [0.15 0.15 0.15];
    weight_ls = '--';
    weight_lw = 2.4;
    model_names_ov = {'Targeted Local Projection', 'Vector Autoregression', ...
        'Local Projection', 'Weight (LP share)'};
    k_ov = 3;

    fig_handle = figure('Units', 'Normalized', 'OuterPosition', [0.05 0.025 0.42 0.88]);
    tiledlayout(n_rows_ov, n_cols_ov, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax_font2 = 16;
    label_font2 = 18;
    title_font2 = 20;

    ylims_ov = compute_overlay_ylims(IRFs, CIs, n_resp);
    W_LP = collect_blp_weights(all_outputs, resp_names, n_h, p_picture);

    plot_order = [3, 2, 1];  % LP, VAR, TLP
    for r = 1:n_resp
        nexttile; hold on;
        yyaxis left;

        alpha_vals = [0.15, 0.15, 0.25];
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

        for idx = 1:k_ov
            m = plot_order(idx);
            plot(x, IRFs{m,r}, 'Color', colors_ov(m,:), ...
                 'LineWidth', lw_ov(m), 'LineStyle', ls_ov{m});
        end

        yline(0, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
        xlim([h_start h_end]);
        ylim(ylims_ov{r});

        yyaxis right;
        w_lp = W_LP(:, r);
        valid_w = isfinite(w_lp);
        plot(x(valid_w), w_lp(valid_w), 'Color', weight_color, ...
            'LineWidth', weight_lw, 'LineStyle', weight_ls);
        draw_corner_weight_guides(w_lp, valid_w);
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

    nexttile;
    axis off;
    hold on;
    h_leg = gobjects(k_ov + 1, 1);
    for m = 1:k_ov
        h_leg(m) = plot(NaN, NaN, 'Color', colors_ov(m,:), ...
            'LineWidth', lw_ov(m), 'LineStyle', ls_ov{m});
    end
    h_leg(k_ov + 1) = plot(NaN, NaN, 'Color', weight_color, ...
        'LineWidth', weight_lw, 'LineStyle', weight_ls);

    legend(h_leg, model_names_ov, ...
        'FontSize', label_font2, ...
        'Box', 'off', ...
        'Location', 'east');

    drawnow;
    align_overlay_axes(fig_handle);

    set(gcf, 'Color', 'w');
    sgtitle(sprintf('Response to 100 bp FFR Shock (%.0f–%.0f)', floor(tStart), floor(tEnd)), ...
        'FontSize', 24, 'FontWeight', 'bold');

    if save_picture
        fname = fullfile(output_dir, sprintf('BLP_application_overlay_4x2_pL%d_pV%d', P_LP, P_VAR));
        exportgraphics(gcf, [fname '.pdf'], 'ContentType', 'vector');
        fprintf('Saved: %s.pdf\n', fname);
    end
    if close_picture
        close(fig_handle);
    end
end

function [IRFs, CIs] = collect_blp_overlay_data(all_outputs, resp_names, n_h, p_picture)
    k_ov = 3;
    n_resp = numel(resp_names);
    IRFs = cell(k_ov, n_resp);
    CIs = cell(k_ov, n_resp);

    for r = 1:n_resp
        output = all_outputs.(resp_names{r});

        d = output.TLP2.Studentized.method7;
        IRFs{1,r} = d.irf(1:n_h, 1, p_picture);
        CIs{1,r} = squeeze(d.CI(1:n_h, :, 1, p_picture));

        d = output.VAR2.Studentized.method7;
        IRFs{2,r} = d.irf(1:n_h, 1, p_picture);
        CIs{2,r} = squeeze(d.CI(1:n_h, :, 1, p_picture));

        d = output.LP2.Studentized.method7;
        IRFs{3,r} = d.irf(1:n_h, 1, p_picture);
        CIs{3,r} = squeeze(d.CI(1:n_h, :, 1, p_picture));
    end
end

function W_LP = collect_blp_weights(all_outputs, resp_names, n_h, p_picture)
    n_resp = numel(resp_names);
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
end

function ylims_ov = compute_overlay_ylims(IRFs, CIs, n_resp)
    k_ov = 3;
    ylims_ov = cell(1, n_resp);
    for r = 1:n_resp
        ylo = inf;
        yhi = -inf;
        for m = 1:k_ov
            vals = [IRFs{m,r}; CIs{m,r}(:)];
            vals = vals(isfinite(vals));
            if ~isempty(vals)
                ylo = min(ylo, min(vals));
                yhi = max(yhi, max(vals));
            end
        end
        pad = 0.10 * (yhi - ylo);
        if pad == 0 || ~isfinite(pad)
            pad = 0.1;
        end
        ylims_ov{r} = [ylo - pad, yhi + pad];
    end
end

function draw_corner_weight_guides(w_lp, valid)
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
end

function align_overlay_axes(fig_handle)
    all_ax = flipud(findobj(fig_handle, 'Type', 'axes'));
    plot_axes = [];
    for i = 1:length(all_ax)
        if ~isempty(get(all_ax(i), 'Children')) || strcmp(get(all_ax(i), 'Visible'), 'on')
            if strcmp(get(all_ax(i), 'Visible'), 'off')
                continue;
            end
            plot_axes = [plot_axes; all_ax(i)];
        end
    end

    left_ax = [];
    right_ax = [];
    for i = 1:length(plot_axes)
        pos = plot_axes(i).Position;
        if pos(1) < 0.5
            left_ax = [left_ax; plot_axes(i)];
        else
            right_ax = [right_ax; plot_axes(i)];
        end
    end

    if ~isempty(left_ax)
        max_x = max(arrayfun(@(a) a.Position(1), left_ax));
        min_w = min(arrayfun(@(a) a.Position(3), left_ax));
        for i = 1:length(left_ax)
            p = left_ax(i).Position;
            left_ax(i).Position = [max_x, p(2), min_w, p(4)];
        end
    end

    if ~isempty(right_ax)
        min_x = min(arrayfun(@(a) a.Position(1), right_ax));
        min_w = min(arrayfun(@(a) a.Position(3), right_ax));
        for i = 1:length(right_ax)
            p = right_ax(i).Position;
            right_ax(i).Position = [min_x, p(2), min_w, p(4)];
        end
    end
end
