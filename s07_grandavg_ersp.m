% s07_grandavg_ersp.m  --  Grand-average ERSP (feasibility figure).
% For each SUBJECT: pool good (correct, non-bad) maintenance trials across runs,
% compute the FCz ERSP locked to the 'WM' (dash) onset, baseline-corrected to the
% pre-dash period. Then average the per-subject ERSPs -> grand average over N.
% Load conditions are collapsed: this shows canonical WM dynamics vs baseline,
% which is exactly the feasibility claim.

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

chan      = 'FCz';                 % frontal-midline
epoch_s   = [-1.0 3.5];            % around WM onset (s)
baseline  = [-500 0];             % ms, pre-dash
freqrange = [2 35];
cycles    = [3 0.5];

subs = unique(manifest.subject,'stable');
GA = []; times = []; freqs = []; usedN = 0; subjList = {};

for su = 1:numel(subs)
    subj = subs(su);
    keys = manifest.key(manifest.subject == subj);

    % gather single-channel maintenance epochs across this subject's runs
    data = []; eplen = []; srate = [];
    for k = 1:numel(keys)
        key = char(keys(k));
        EEG = pop_loadset('filename',[key '_clean.set'],'filepath',cfg.deriv_dir);
        T = load(fullfile(cfg.deriv_dir,[key '_trials.mat'])); trials = T.trials;
        srate = EEG.srate;
        ep = round(epoch_s * srate); if isempty(eplen), eplen = ep(2)-ep(1)+1; end
        ch = find(strcmpi({EEG.chanlocs.labels}, chan), 1);
        sel = trials([trials.good]);          % collapse load conditions
        for t = 1:numel(sel)
            a = sel(t).lock_sample + ep(1); b = a + eplen - 1;
            if a < 1 || b > EEG.pnts, continue; end
            data = [data, EEG.data(ch, a:b).']; %#ok<AGROW>  % column per trial
        end
    end
    ntr = size(data,2);
    if ntr < 5, warning('%s: only %d trials, skipping', subj, ntr); continue; end

    [ersp,~,~,t_ms,f_hz] = newtimef(data, eplen, epoch_s*1000, srate, cycles, ...
        'baseline', baseline, 'freqs', freqrange, ...
        'plotersp','off','plotitc','off','verbose','off');

    if isempty(GA)
        GA = zeros(size(ersp,1), size(ersp,2), 0); times = t_ms; freqs = f_hz;
    end
    if ~isequal(size(ersp,1),numel(freqs)) || ~isequal(size(ersp,2),numel(times))
        warning('%s: ERSP grid mismatch, skipping', subj); continue;
    end
    GA(:,:,end+1) = ersp; %#ok<AGROW>
    usedN = usedN + 1; subjList{end+1} = char(subj); %#ok<AGROW>
    fprintf('ERSP done %s  (%d trials)\n', subj, ntr);
end

grand = mean(GA, 3);

% ---- figure ----
f = figure('Color','w','Position',[100 100 720 460]);
imagesc(times, freqs, grand); axis xy; hold on;
cmax = max(abs(grand(:))); clim([-cmax cmax]); colormap(jet);
cb = colorbar; ylabel(cb, 'ERSP (dB)');
plot([0 0], freqrange, 'k--', 'LineWidth', 1.5);             % maintenance onset
yline(4,'w:'); yline(7,'w:'); yline(13,'w:'); yline(30,'w:'); % band edges
xlabel('Time from maintenance onset (ms)'); ylabel('Frequency (Hz)');
title(sprintf('Grand-average ERSP at %s  (N = %d), locked to WM maintenance', chan, usedN));

outdir = fullfile(cfg.deriv_dir,'figures'); if ~exist(outdir,'dir'), mkdir(outdir); end
exportgraphics(f, fullfile(outdir,'fig_grandavg_ersp.png'), 'Resolution', cfg.fig.dpi);
save(fullfile(cfg.deriv_dir,'grandavg_ersp.mat'), 'grand','times','freqs','subjList','usedN');
fprintf('\nGrand-average ERSP over %d subjects saved.\n', usedN);
