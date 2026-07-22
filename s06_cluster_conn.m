% s06_cluster_conn.m  --  FWER-controlled test across the full connectivity space.
% Max-statistic sign-flip permutation over all frontal x temporal pairs x bands
% (Maris-Oostenveld family; single-element clusters, no arbitrary pair adjacency).
% Answers honestly: is there ANY localised FT dwPLI load effect that the
% pair-averaged test in s04 may have diluted?  Controls family-wise error at 0.05.

clear; clc;
S = load(fullfile(config().deriv_dir,'conn_results.mat'));  % gives R and cfg
R = S.R; cfg = S.cfg;
rng(1);

bands = {'theta','beta'};
nF = numel(cfg.roi_frontal); nT = numel(cfg.roi_temporal); np = nF*nT;
subs = unique(R.subject,'stable');

% ----- build paired difference matrix D : [nsub x (np*nbands)] -----
labs = {}; D = [];
for bi = 1:numel(bands)
    Db = nan(numel(subs), np);
    for si = 1:numel(subs)
        rl = R(strcmp(R.subject,subs{si}) & strcmp(R.band,bands{bi}) & strcmp(R.load,'low'), :);
        rh = R(strcmp(R.subject,subs{si}) & strcmp(R.band,bands{bi}) & strcmp(R.load,'high'), :);
        if isempty(rl) || isempty(rh), continue; end
        lo = rl.dwpli_pairs{1}; if iscell(lo), lo = lo{1}; end
        hi = rh.dwpli_pairs{1}; if iscell(hi), hi = hi{1}; end
        Db(si,:) = hi - lo;
    end
    D = [D Db]; %#ok<AGROW>
    for i = 1:nF, for j = 1:nT
        labs{end+1} = sprintf('%-5s %s-%s', bands{bi}, cfg.roi_frontal{i}, cfg.roi_temporal{j}); %#ok<AGROW,SAGROW>
    end, end %#ok<ALIGN>
end

bad = any(isnan(D),2); nDrop = sum(bad); D = D(~bad,:); n = size(D,1); M = size(D,2);
fprintf('Subjects used: %d (dropped %d with incomplete data)\n', n, nDrop);

% ----- observed and permutation statistics -----
tobs = mean(D) ./ (std(D)./sqrt(n));
nperm = cfg.stat.n_perm;
tnull = zeros(nperm, M); maxnull = zeros(nperm,1);
for p = 1:nperm
    sgn = sign(rand(n,1)-0.5);
    Dp = D .* sgn;                          % same flip across all variables
    t = mean(Dp) ./ (std(Dp)./sqrt(n));
    tnull(p,:) = t; maxnull(p) = max(abs(t));
end
thr     = prctile(maxnull, 95);
p_unc   = mean(abs(tnull) >= abs(tobs), 1);          % uncorrected
p_fwer  = mean(maxnull   >= abs(tobs), 1);           % max-stat corrected
sig     = p_fwer < 0.05;

% ----- report -----
T = table(labs(:), tobs(:), p_unc(:), p_fwer(:), sig(:), ...
    'VariableNames', {'pair_band','t','p_uncorr','p_fwer','sig_fwer'});
T = sortrows(T, 't', 'descend', 'ComparisonMethod','abs');
fprintf('\nFWER threshold |t| (95%% of max-null) = %.2f\n', thr);
fprintf('FWER-significant pairs: %d of %d\n\n', sum(sig), M);
disp(T(1:min(12,height(T)), :));
writetable(T, fullfile(cfg.deriv_dir,'conn_cluster_results.csv'));

% ----- t-map figure (frontal x temporal), one panel per band -----
f = figure('Color','w','Position',[100 100 820 360]);
for bi = 1:numel(bands)
    cols = (bi-1)*np + (1:np);
    Tm = reshape(tobs(cols), nT, nF)';        % rows = frontal, cols = temporal
    Sm = reshape(sig(cols),  nT, nF)';
    subplot(1,2,bi);
    imagesc(Tm); axis image; colorbar; clim([-max(abs(tobs)) max(abs(tobs))]);
    colormap(parula);
    set(gca,'XTick',1:nT,'XTickLabel',cfg.roi_temporal, ...
            'YTick',1:nF,'YTickLabel',cfg.roi_frontal);
    title(sprintf('%s band: paired t (High-Low)', bands{bi}));
    [ri,ci] = find(Sm);                        % mark significant cells
    hold on; plot(ci, ri, 'k*', 'MarkerSize',10, 'LineWidth',1.5);
end
outdir = fullfile(cfg.deriv_dir,'figures'); if ~exist(outdir,'dir'), mkdir(outdir); end
exportgraphics(f, fullfile(outdir,'fig3_conn_tmap.png'), 'Resolution', cfg.fig.dpi);
fprintf('\nSaved conn_cluster_results.csv and fig3_conn_tmap.png\n');
