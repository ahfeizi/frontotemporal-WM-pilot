% s10_topography.m -- REVISED (v2): fixed for topoplot's 'blank' style API.
clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
T = readtable(fullfile(cfg.deriv_dir,'conn_cluster_results.csv'));
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

key = char(manifest.key(1));
EEG = pop_loadset('filename',[key '_clean.set'],'filepath',cfg.deriv_dir);
labels = {EEG.chanlocs.labels};
roi_chans = [cfg.roi_frontal, cfg.roi_temporal];

bands = {'theta','beta'};
f = figure('Color','w','Position',[100 100 800 380]);
for b = 1:numel(bands)
    bn = bands{b};
        tvals = nan(1, numel(labels));
    tcounts = zeros(1, numel(labels));
    rows = contains(T.pair_band, bn);
    Tb = T(rows,:);
    for r = 1:height(Tb)
        parts = strsplit(strtrim(Tb.pair_band{r}));
        pair = parts{2};
        chs = strsplit(pair,'-');
        for c = 1:2
            idx = find(strcmpi(labels, chs{c}));
            if ~isempty(idx) && any(strcmpi(labels{idx}, roi_chans))
                if isnan(tvals(idx)), tvals(idx) = 0; end
                tvals(idx) = tvals(idx) + Tb.t(r);
                tcounts(idx) = tcounts(idx) + 1;
            end
        end
    end
    tvals(tcounts>0) = tvals(tcounts>0) ./ tcounts(tcounts>0);

    subplot(1,2,b);
    plotIdx = find(~isnan(tvals));

    % FIX: for 'style','blank', pass the channel-number VECTOR itself as the
    % first argument (this is how topoplot marks specific channels in blank
    % mode) -- no 'emarker2' here.
    topoplot(plotIdx, EEG.chanlocs, 'style','blank', 'electrodes','off');
    hold on;

    % now manually overlay our own colored, t-value-scaled markers
    cmap = colormap('jet'); vmax = max(abs(tvals(plotIdx)));
    for ci = plotIdx
        [x,y] = pol2cart(EEG.chanlocs(ci).theta/180*pi, EEG.chanlocs(ci).radius);
        cval = tvals(ci);
        cidx = round(((cval+vmax)/(2*vmax))*(size(cmap,1)-1))+1;
        cidx = max(1,min(size(cmap,1),cidx));
        plot(-y, x, 'o', 'MarkerSize', 16, 'MarkerFaceColor', cmap(cidx,:), ...
            'MarkerEdgeColor','k');
        text(-y, x, labels{ci}, 'FontSize', 7, 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom');
    end
    caxis([-vmax vmax]); colorbar;
    title(sprintf('%s band: analyzed channels only (t, High-Low)', bn));
end
outdir = fullfile(cfg.deriv_dir,'figures'); if ~exist(outdir,'dir'), mkdir(outdir); end
exportgraphics(f, fullfile(outdir,'fig4_topography.png'), 'Resolution', cfg.fig.dpi);