% s03_connectivity.m  --  THE SCIENTIFIC CORE (subject-level).
% Aggregates each subject's datasets (runs) from the manifest. For every
% subject / band / load condition it POOLS good-trial maintenance segments
% ACROSS that subject's runs, then computes:
%   dwPLI : debiased weighted phase-lag index (Vinck et al., 2011)  [PRIMARY]
%   PLV   : phase-locking value                                     [SECONDARY]
% over all frontal x temporal channel pairs (across trials per time sample,
% averaged over the maintenance window). One row per subject x band x load.
% => N for statistics = number of SUBJECTS.

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

subjects = unique(manifest.subject, 'stable');
bands = fieldnames(cfg.bands);
conds = {'low','high'};
rows  = {};

for su = 1:numel(subjects)
    subj = subjects(su);
    keys = manifest.key(manifest.subject == subj);

    % accumulate pooled analytic segments per band/cond across this subject's runs
    acc = struct();
    for b = 1:numel(bands), for c = 1:2, acc.(bands{b}).(conds{c}) = []; end, end %#ok<ALIGN>
    refEEG = [];

    for k = 1:numel(keys)
        key = char(keys(k));
        EEG = pop_loadset('filename',[key '_clean.set'],'filepath',cfg.deriv_dir);
        T = load(fullfile(cfg.deriv_dir,[key '_trials.mat'])); trials = T.trials;
        if isempty(refEEG), refEEG = EEG; end
        srate = EEG.srate; win = round(cfg.maint_window*srate);
        wlen = win(2)-win(1)+1;

        for b = 1:numel(bands)
            Xf = util_bandpass(double(EEG.data), srate, cfg.bands.(bands{b}), cfg.filt_order_factor);
            Z  = hilbert(Xf.').';
            for c = 1:2
                sel = trials([trials.good] & strcmp({trials.load},conds{c}));
                for t = 1:numel(sel)
                    a = sel(t).lock_sample + win(1); z = a + wlen - 1;
                    if a < 1 || z > size(Z,2), continue; end
                    acc.(bands{b}).(conds{c}) = cat(3, acc.(bands{b}).(conds{c}), Z(:, a:z));
                end
            end
        end
    end
    if isempty(refEEG), continue; end

    fr = util_chanidx(refEEG, cfg.roi_frontal);
    te = util_chanidx(refEEG, cfg.roi_temporal);
    fmIdx = util_chanidx(refEEG, {'FCz','Fz'}); fmIdx = fmIdx(1);

    for b = 1:numel(bands)
        for c = 1:2
            seg = acc.(bands{b}).(conds{c}); ntr = size(seg,3);
            if ntr < 5
                warning('%s %s %s: %d pooled trials, skipping', subj, bands{b}, conds{c}, ntr);
                continue;
            end
            wlen = size(seg,2);
            np = numel(fr)*numel(te);
            dwpli_pairs = zeros(1,np); plv_pairs = zeros(1,np); p = 0;
            for i = 1:numel(fr)
                zi = reshape(seg(fr(i),:,:), wlen, ntr);
                for j = 1:numel(te)
                    p = p+1;
                    zj = reshape(seg(te(j),:,:), wlen, ntr);
                    X = zi .* conj(zj); im = imag(X);
                    s_im = sum(im,2); s_abs = sum(abs(im),2); s_sq = sum(im.^2,2);
                    dwpli_t = (s_im.^2 - s_sq) ./ (s_abs.^2 - s_sq);
                    plv_t   = abs(mean(X ./ abs(X), 2));
                    dwpli_pairs(p) = mean(dwpli_t,'omitnan');
                    plv_pairs(p)   = mean(plv_t,'omitnan');
                end
            end
            fm_power = mean(mean(abs(reshape(seg(fmIdx,:,:),wlen,ntr)).^2));
            rows(end+1,:) = {char(subj), conds{c}, bands{b}, ...
                mean(dwpli_pairs,'omitnan'), mean(plv_pairs,'omitnan'), ...
                fm_power, ntr, {dwpli_pairs}, {plv_pairs}}; %#ok<SAGROW>
        end
    end
    fprintf('done %s (%d run-files)\n', subj, numel(keys));
end

R = cell2table(rows, 'VariableNames', {'subject','load','band','dwpli','plv', ...
    'fm_power','ntrials','dwpli_pairs','plv_pairs'});
save(fullfile(cfg.deriv_dir,'conn_results.mat'),'R','cfg');
writetable(R(:,1:7), fullfile(cfg.deriv_dir,'conn_results.csv'));
fprintf('\nSaved conn_results: %d subjects.\n', numel(unique(R.subject)));
