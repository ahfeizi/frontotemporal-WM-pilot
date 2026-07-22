% s05_figures.m  --  Publication-grade figures (300 dpi PNG).
% Fig1a: frontal-midline THETA power, low vs high  (non-significant control)
% Fig1b: frontal-midline BETA power, low vs high   (headline local effect)
% Fig2 : frontotemporal dwPLI, theta & beta, low vs high  (headline connectivity)
% Each panel shows per-subject paired lines, group mean +- SEM, and permutation p.

clear; clc;
cfg = config();
load(fullfile(cfg.deriv_dir,'conn_results.mat'), 'R');
load(fullfile(cfg.deriv_dir,'stats_results.mat'), 'results');
outdir = fullfile(cfg.deriv_dir,'figures'); if ~exist(outdir,'dir'), mkdir(outdir); end
subs = unique(R.subject);

val = @(bn,mt,ld) arrayfun(@(i) getval(R, subs{i}, bn, mt, ld), 1:numel(subs));

% ---- Figure 1: frontal-midline power, both bands (theta = control, beta = headline) ----
f1 = figure('Color','w','Position',[100 100 760 380]);

% -- 1a: theta (non-significant control) --
subplot(1,2,1);
lo_t = val('theta','fm_power','low'); hi_t = val('theta','fm_power','high');
if cfg.log_power
    lo_t = log(lo_t); hi_t = log(hi_t); ylab_t = 'log frontal-midline \theta power';
else
    ylab_t = 'Frontal-midline \theta power (\muV^2)';
end
pairplot(gca, lo_t, hi_t, ylab_t, results.theta.fm_power.p, {'Low','High'});
title('\theta power (control, n.s.)');

% -- 1b: beta (headline local effect) --
subplot(1,2,2);
lo_b = val('beta','fm_power','low'); hi_b = val('beta','fm_power','high');
if cfg.log_power
    lo_b = log(lo_b); hi_b = log(hi_b); ylab_b = 'log frontal-midline \beta power';
else
    ylab_b = 'Frontal-midline \beta power (\muV^2)';
end
pairplot(gca, lo_b, hi_b, ylab_b, results.beta.fm_power.p, {'Low','High'});
title('\beta power (headline load effect)');

exportgraphics(f1, fullfile(outdir,'fig1_fm_power.png'), 'Resolution', cfg.fig.dpi);

% ---- Figure 2: headline frontotemporal dwPLI ---- (unchanged)
f2 = figure('Color','w','Position',[100 100 760 380]);
bn = {'theta','beta'};
for k = 1:2
    subplot(1,2,k);
    pairplot(gca, val(bn{k},'dwpli','low'), val(bn{k},'dwpli','high'), ...
        sprintf('Frontotemporal dwPLI'), results.(bn{k}).dwpli.p, {'Low','High'});
    title(sprintf('%s band', bn{k}));
end
exportgraphics(f2, fullfile(outdir,'fig2_ft_dwpli.png'), 'Resolution', cfg.fig.dpi);

fprintf('Figures saved to %s\n', outdir);

% ======================= local functions =======================
function v = getval(R, sub, bn, mt, ld)
    r = R(strcmp(R.subject,sub) & strcmp(R.band,bn) & strcmp(R.load,ld), :);
    if isempty(r), v = NaN; else, v = r.(mt)(1); end
end

function pairplot(ax, lo, hi, ylab, pval, xl)
    lo = lo(:); hi = hi(:); ok = ~isnan(lo) & ~isnan(hi);
    lo = lo(ok); hi = hi(ok); n = numel(lo);
    hold(ax,'on');
    for i = 1:n
        plot(ax, [1 2], [lo(i) hi(i)], '-', 'Color', [.7 .7 .7]);
    end
    plot(ax, ones(n,1), lo, 'o', 'MarkerFaceColor',[.3 .5 .8], 'MarkerEdgeColor','none');
    plot(ax, 2*ones(n,1), hi, 'o', 'MarkerFaceColor',[.85 .4 .3], 'MarkerEdgeColor','none');
    m = [mean(lo) mean(hi)]; e = [std(lo) std(hi)]/sqrt(n);
    errorbar(ax, [1 2], m, e, 'k', 'LineWidth', 2, 'CapSize', 12);
    xlim(ax,[.6 2.4]); set(ax,'XTick',[1 2],'XTickLabel',xl);
    ylabel(ax, ylab); box(ax,'off');
    yl = ylim(ax); text(ax, 1.5, yl(2), sprintf('p = %.3f', pval), ...
        'HorizontalAlignment','center','VerticalAlignment','top');
end