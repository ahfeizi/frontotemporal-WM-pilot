% E_global_null.m
% Global null for the count of significant cells: 5,000 random sign-flip
% patterns (seed 1) applied jointly to all 48 cells (preserving their
% dependence); per pattern, count cells in which any pair passes the
% within-cell max-|t| FWER threshold at alpha = 0.05; report the
% distribution and P(count >= 2).
%
% This replaces the "expectation of 2.4" in §4.2.
%
% Per-subject per-pair contrast values are computed here directly from the
% preprocessed EEG files (no dependency on battery_persubject.csv, which
% accumulates data across multiple datasets).
%
% Strategy:
%  1. Load EEG and compute Hilbert + wavelet connectivity for all 48 cells.
%  2. Store [22 x 30 x 48] contrast matrix (high - low, ROI pairs).
%  3. For 5,000 random sign-flip patterns (seed 1, same flip for all cells):
%     a. Apply flip to per-participant contrasts.
%     b. For each cell x pair, compute t from flipped values.
%     c. For each cell, record max|t| over 30 pairs.
%  4. Derive cell-specific FWER threshold (95th percentile).
%  5. Count cells passing per pattern; report P(count >= 2).
%
% Reads:  <deriv_dir>/<key>_clean.set + <key>_trials.mat (per run)
% Writes: revision/E_global_null.csv
%         revision/E_global_null_dist.csv
%         revision/E_global_null.md

clear; clc;
repo = fileparts(mfilename('fullpath'));
addpath(fileparts(repo));
cfg  = config();

SEED   = 1;
N_PERM = 5000;
ALPHA  = 0.05;

addpath(cfg.eeglab_dir); eeglab nogui;
addpath(cfg.fieldtrip_dir); ft_defaults;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

n_cycles = cfg.battery.n_cycles;
REFS  = {'average','csd'};
CONDS = {'low','high'};
BANDS = fieldnames(cfg.battery.bands);
N_BANDS = numel(BANDS);
MEAS  = {'PLV','dwPLI','ImC'};
DECOMPS = {'hilbert','wavelet'};

% All wavelet centre frequencies (union across bands)
all_freqs = [];
for b = 1:N_BANDS, all_freqs = [all_freqs cfg.battery.wavelet_freqs.(BANDS{b})]; end %#ok<AGROW>
all_freqs = unique(all_freqs);

% Cell ordering: 48 cells = 4 bands x 3 measures x 2 decomps x 2 refs
cell_defs = {};
for b = 1:N_BANDS
    for d = 1:2, for r = 1:2, for m = 1:3
        cell_defs{end+1} = struct('band',BANDS{b},'measure',MEAS{m}, ...
            'decomposition',DECOMPS{d},'reference',REFS{r}); %#ok<SAGROW>
    end, end, end
end
n_cells = numel(cell_defs);  % 48

subjects = unique(manifest.subject, 'stable');
% Restrict to the N=22 analysed set (manifest.mat has 23; sub-005 was
% excluded from the original battery -- see A_trial_equalisation.m).
excl_T = readtable(fullfile(fileparts(repo), 'results', 'battery_persubject.csv'));
included_subj = unique(string(excl_T.subject));
subjects = subjects(ismember(subjects, included_subj));
N = numel(subjects);
fprintf('N = %d subjects, %d cells\n', N, n_cells);

% Will determine n_pairs from first subject
n_pairs = [];
pairLabs = {};

% X: [N x n_pairs x n_cells], contrast (high - low) per subject per pair per cell
X = [];

for su = 1:N
    subj = char(subjects(su));
    t0 = tic;
    keys = manifest.key(manifest.subject == subjects(su));

    % Load all runs once
    runs = struct('srate',{},'sel_lo',{},'sel_hi',{},'average',{},'csd',{});
    for k = 1:numel(keys)
        key = char(keys(k));
        EEG = pop_loadset('filename',[key '_clean.set'],'filepath',cfg.deriv_dir);
        T_t = load(fullfile(cfg.deriv_dir,[key '_trials.mat'])); trials_k = T_t.trials;
        if k == 1 && su == 1
            have = {EEG.chanlocs.labels};
            labF = have(cellfun(@(l) any(strcmpi(have,l)), cfg.roi_frontal));
            labT = have(cellfun(@(l) any(strcmpi(have,l)), cfg.roi_temporal));
            roiLabs = [labF labT];
            nF = numel(labF); nT = numel(labT);
            fr = 1:nF; te = nF + (1:nT);
            n_pairs = nF * nT;
            X = nan(N, n_pairs, n_cells);
            for i = 1:nF
                for j = 1:nT
                    pairLabs{end+1} = sprintf('%s-%s', labF{i}, labT{j}); %#ok<SAGROW>
                end
            end
        end
        [csd_d, csd_lbs] = local_csd(EEG, key, cfg.fieldtrip_dir);
        iAvg = util_chanidx(EEG, roiLabs);
        iCsd = util_chanidx(local_labels2eeg(csd_lbs), roiLabs);
        runs(k).srate   = EEG.srate;
        [runs(k).sel_lo, runs(k).sel_hi] = local_sel_cond(trials_k, key);
        runs(k).average = double(EEG.data(iAvg,:));
        runs(k).csd     = double(csd_d(iCsd,:));
        clear EEG csd_d
    end
    if isempty(runs), continue; end

    for r = 1:numel(REFS)
        rf = REFS{r};

        % Hilbert: one set of segments per band
        for b = 1:N_BANDS
            band = cfg.battery.bands.(BANDS{b});
            s_lo = local_hil_segs(runs, rf, cfg.maint_window, band, cfg.filt_order_factor, 'lo');
            s_hi = local_hil_segs(runs, rf, cfg.maint_window, band, cfg.filt_order_factor, 'hi');
            if size(s_lo,3) < 5 || size(s_hi,3) < 5, continue; end
            [plv_lo,dw_lo,imc_lo,~] = util_conn_pairs(s_lo, fr, te);
            [plv_hi,dw_hi,imc_hi,~] = util_conn_pairs(s_hi, fr, te);
            for m = 1:3
                cid = local_cid(cell_defs, BANDS{b}, MEAS{m}, 'hilbert', rf);
                if cid == 0, continue; end
                switch m
                    case 1, X(su,:,cid) = (plv_hi - plv_lo)';
                    case 2, X(su,:,cid) = (dw_hi  - dw_lo)';
                    case 3, X(su,:,cid) = (imc_hi - imc_lo)';
                end
            end
        end

        % Wavelet: per-frequency then average within band
        % Precompute per-freq connectivity
        nf = numel(all_freqs);
        plv_lo_f = nan(nf, n_pairs); plv_hi_f = nan(nf, n_pairs);
        dw_lo_f  = nan(nf, n_pairs); dw_hi_f  = nan(nf, n_pairs);
        imc_lo_f = nan(nf, n_pairs); imc_hi_f = nan(nf, n_pairs);
        ok_lo = false(1,nf); ok_hi = false(1,nf);
        for fi = 1:nf
            fc = all_freqs(fi);
            wlo = local_wav_segs(runs, rf, cfg.maint_window, fc, n_cycles, 'lo');
            whi = local_wav_segs(runs, rf, cfg.maint_window, fc, n_cycles, 'hi');
            if size(wlo,3) >= 5
                [plv_lo_f(fi,:),dw_lo_f(fi,:),imc_lo_f(fi,:),~] = util_conn_pairs(wlo, fr, te);
                ok_lo(fi) = true;
            end
            if size(whi,3) >= 5
                [plv_hi_f(fi,:),dw_hi_f(fi,:),imc_hi_f(fi,:),~] = util_conn_pairs(whi, fr, te);
                ok_hi(fi) = true;
            end
        end
        for b = 1:N_BANDS
            band_freqs = cfg.battery.wavelet_freqs.(BANDS{b});
            fb = find(ismember(all_freqs, band_freqs));
            used = fb(ok_lo(fb) & ok_hi(fb));
            if isempty(used), continue; end
            nu = numel(used);
            for m = 1:3
                cid = local_cid(cell_defs, BANDS{b}, MEAS{m}, 'wavelet', rf);
                if cid == 0, continue; end
                switch m
                    case 1, X(su,:,cid) = (mean(plv_hi_f(used,:),1) - mean(plv_lo_f(used,:),1))';
                    case 2, X(su,:,cid) = (mean(dw_hi_f(used,:),1)  - mean(dw_lo_f(used,:),1))';
                    case 3, X(su,:,cid) = (mean(imc_hi_f(used,:),1) - mean(imc_lo_f(used,:),1))';
                end
            end
        end
    end
    fprintf('  %s done (%.1f min)\n', subj, toc(t0)/60);
end
fprintf('Contrast matrix: %d x %d x %d\n', N, n_pairs, n_cells);

% ── Observed t-statistics ─────────────────────────────────────────────────────
t_obs = local_tstat(X);        % [n_pairs x n_cells]
maxt_obs = max(abs(t_obs), [], 1);  % [1 x n_cells]

% ── Sign-flip permutations ────────────────────────────────────────────────────
rng(SEED);
maxt_perm = nan(N_PERM, n_cells);

fprintf('Running %d sign-flip permutations (seed %d) ...\n', N_PERM, SEED);
for b = 1:N_PERM
    signs = sign(rand(N,1) - 0.5);   % [N x 1], ±1
    Xb = X .* signs;                  % broadcast over pairs and cells
    tb = local_tstat(Xb);
    maxt_perm(b,:) = max(abs(tb), [], 1);
    if mod(b,500) == 0, fprintf('  %d / %d done\n', b, N_PERM); end
end

% Cell-specific FWER threshold (95th percentile under H0)
thresh = prctile(maxt_perm, 100*(1-ALPHA), 1);  % [1 x n_cells]

% Count of cells passing per permuted pattern
pass_perm = sum(maxt_perm > thresh, 2);   % [N_PERM x 1]
pass_obs  = sum(maxt_obs  > thresh);
p_count_ge2 = mean(pass_perm >= 2);

fprintf('\nObserved: %d cells pass within-cell FWER.\n', pass_obs);
for k = 0:min(6, max(pass_perm))
    fprintf('  P(count = %d) = %.4f\n', k, mean(pass_perm == k));
end
fprintf('P(count >= 2) = %.4f\n', p_count_ge2);

% ── Write outputs ─────────────────────────────────────────────────────────────
cell_band  = cellfun(@(c) c.band,          cell_defs, 'UniformOutput', false)';
cell_meas  = cellfun(@(c) c.measure,       cell_defs, 'UniformOutput', false)';
cell_dec   = cellfun(@(c) c.decomposition, cell_defs, 'UniformOutput', false)';
cell_ref   = cellfun(@(c) c.reference,     cell_defs, 'UniformOutput', false)';
out_cells  = table(cell_band, cell_meas, cell_dec, cell_ref, ...
    maxt_obs', thresh', (maxt_obs > thresh)', ...
    'VariableNames', {'band','measure','decomposition','reference', ...
                      'maxt_obs','fwer_thresh','obs_pass'});
writetable(out_cells, fullfile(repo, 'E_global_null.csv'));

k_vals = (0:max(pass_perm))';
dist_T = table(k_vals, ...
    arrayfun(@(k) sum(pass_perm==k), k_vals), ...
    arrayfun(@(k) mean(pass_perm==k), k_vals), ...
    'VariableNames', {'n_cells_pass','count','prob'});
writetable(dist_T, fullfile(repo, 'E_global_null_dist.csv'));

% Markdown
md = {};
md{end+1} = '# E. Global null: count of cells with within-cell FWER hits';
md{end+1} = '';
md{end+1} = sprintf('N = %d subjects; %d cells x %d pairs; seed = %d; N_perm = %d.', ...
    N, n_cells, n_pairs, SEED, N_PERM);
md{end+1} = 'Per-subject contrasts computed from EEG files directly (not battery_persubject.csv).';
md{end+1} = sprintf('Within-cell threshold: 95th pct of max|t| over %d pairs under H0', n_pairs);
md{end+1} = sprintf('(from same %d sign-flip patterns; compare exact ~3.67 in §4.2).', N_PERM);
md{end+1} = '';
md{end+1} = sprintf('Observed: **%d** cells pass FWER (cf. manuscript: 2 cells).', pass_obs);
md{end+1} = sprintf('Cell-specific thresholds: min = %.3f, max = %.3f, median = %.3f.', ...
    min(thresh), max(thresh), median(thresh));
md{end+1} = '';
md{end+1} = '## Distribution of count under joint null';
for k = 0:min(5, max(pass_perm))
    md{end+1} = sprintf('  P(count = %d) = %.4f', k, mean(pass_perm == k));
end
md{end+1} = sprintf('**P(count >= 2) = %.4f** (replaces "expectation of 2.4" in §4.2)', p_count_ge2);
md{end+1} = '';
md{end+1} = 'Outputs: E_global_null.csv (per-cell), E_global_null_dist.csv (distribution).';
fid = fopen(fullfile(repo, 'E_global_null.md'), 'w');
fprintf(fid, '%s\n', md{:});
fclose(fid);
fprintf('Done.\n');

% ===================== local functions =====================
function t = local_tstat(X)
    n = size(X,1);
    mn = squeeze(mean(X, 1, 'omitnan'));
    sd = squeeze(std(X, 0, 1));
    t  = mn ./ (sd / sqrt(n));
end

function cid = local_cid(cell_defs, band, measure, decomp, ref)
    cid = 0;
    for i = 1:numel(cell_defs)
        c = cell_defs{i};
        if strcmp(c.band,band) && strcmp(c.measure,measure) && ...
           strcmp(c.decomposition,decomp) && strcmp(c.reference,ref)
            cid = i; return;
        end
    end
end

function seg = local_hil_segs(runs, rf, maint_window, band, ord_fac, cond_str)
    seg = [];
    for k = 1:numel(runs)
        srate = runs(k).srate; win = round(maint_window*srate); wlen = win(2)-win(1)+1;
        Xf = util_bandpass(runs(k).(rf), srate, band, ord_fac);
        Z  = hilbert(Xf.').';
        sel = runs(k).(['sel_' cond_str]);
        tmp = complex(zeros(size(Z,1), wlen, numel(sel))); kept = 0;
        for t = 1:numel(sel)
            a = sel(t).lock_sample + win(1); z = a + wlen - 1;
            if a < 1 || z > size(Z,2), continue; end
            kept = kept + 1; tmp(:,:,kept) = Z(:,a:z);
        end
        seg = cat(3, seg, tmp(:,:,1:kept));
    end
end

function seg = local_wav_segs(runs, rf, maint_window, fc, nc, cond_str)
    sigma_t = nc / (2*pi*fc);
    seg = [];
    for k = 1:numel(runs)
        srate = runs(k).srate; win = round(maint_window*srate); wlen = win(2)-win(1)+1;
        t_wav = -3*sigma_t : 1/srate : 3*sigma_t;
        wav = exp(2i*pi*fc.*t_wav) .* exp(-t_wav.^2/(2*sigma_t^2));
        wav = wav / sqrt(sum(abs(wav).^2));
        Z = zeros(size(runs(k).(rf)));
        for ch = 1:size(Z,1), Z(ch,:) = conv(runs(k).(rf)(ch,:), wav, 'same'); end
        sel = runs(k).(['sel_' cond_str]);
        tmp = complex(zeros(size(Z,1), wlen, numel(sel))); kept = 0;
        for t = 1:numel(sel)
            a = sel(t).lock_sample + win(1); z = a + wlen - 1;
            if a < 1 || z > size(Z,2), continue; end
            kept = kept + 1; tmp(:,:,kept) = Z(:,a:z);
        end
        seg = cat(3, seg, tmp(:,:,1:kept));
    end
end

function [csd_data, csd_labels] = local_csd(EEG, key, ft_dir)
    addpath(ft_dir);
    noPos = arrayfun(@(c) isempty(c.X)||isempty(c.Y)||isempty(c.Z)||...
        any(~isfinite(double([c.X c.Y c.Z]))), EEG.chanlocs);
    if any(noPos)
        error('E_global_null:no_positions','%s: CSD requires 3-D positions', key);
    end
    ft_raw = eeglab2fieldtrip(EEG,'raw','none');
    cfg_c = []; cfg_c.method = 'spline';
    ft_csd = ft_scalpcurrentdensity(cfg_c, ft_raw);
    csd_data = ft_csd.trial{1}; csd_labels = ft_csd.label;
end

function E = local_labels2eeg(labels)
    E.chanlocs = struct('labels', reshape(labels,1,[]));
end

function [sel_lo, sel_hi] = local_sel_cond(trials, key)
    if isfield(trials, 'load')
        sel_lo = trials([trials.good] & strcmp({trials.load}, 'low'));
        sel_hi = trials([trials.good] & strcmp({trials.load}, 'high'));
    elseif isfield(trials, 'condition')
        sel_lo = trials([trials.good] & strcmp({trials.condition}, 'low'));
        sel_hi = trials([trials.good] & strcmp({trials.condition}, 'high'));
    elseif isfield(trials, 'wm_load')
        sel_lo = trials([trials.good] & strcmp({trials.wm_load}, 'low'));
        sel_hi = trials([trials.good] & strcmp({trials.wm_load}, 'high'));
    elseif isfield(trials, 'cue') && isnumeric(trials(1).cue)
        % ds004117: cue=1 = low load (set size 3), cue=2 = high load (set size 7)
        cue_vals = [trials.cue];
        sel_lo = trials([trials.good] & cue_vals == 1);
        sel_hi = trials([trials.good] & cue_vals == 2);
    else
        fn = fieldnames(trials);
        error('local_sel_cond: no condition field. key=%s; fields: %s', key, strjoin(fn, ', '));
    end
end
