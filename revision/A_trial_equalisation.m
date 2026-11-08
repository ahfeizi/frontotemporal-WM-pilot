% A_trial_equalisation.m
% Trial-count equalisation robustness check.
% Trial counts per participant range from 9-44 (low) and 11-39 (high).
%
% For 200 subsampling repetitions (seeds 1-200):
%   Each participant: draw min(n_low, n_high) trials from each condition.
%   All 12 pipelines x 4 bands: compute ROI-average connectivity contrast.
%   Seed per (rep, subject) = rep * 1e6 + subject_index (reproducible).
%
% Summary across reps: median dz, 2.5th/97.5th percentile, median BF01.
%
% Also compute PPC (pairwise phase consistency; trial-count-unbiased PLV^2
% estimator) on FULL (unequalized) trials: PPC = (n_tr*PLV^2 - 1)/(n_tr-1).
%
% Reads:  <deriv_dir>/<key>_clean.set + <key>_trials.mat (per run)
% Writes: revision/A_trial_equalisation.csv  (one row per cell, 48 rows)
%         revision/A_ppc_full_trials.csv      (PPC, PLV/Hilbert/average|csd)
%         revision/A_trial_equalisation.md

clear; clc;
repo = fileparts(mfilename('fullpath'));
addpath(fileparts(repo));
cfg  = config();

addpath(cfg.eeglab_dir); eeglab nogui;
addpath(cfg.fieldtrip_dir); ft_defaults;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

N_REPS    = 200;
N_SUBJ    = 22;
n_cycles  = cfg.battery.n_cycles;   % 5
REFS      = {'average','csd'};
CONDS     = {'low','high'};
BANDS     = fieldnames(cfg.battery.bands);
N_BANDS   = numel(BANDS);
MEAS      = {'PLV','dwPLI','ImC'};
DECOMPS   = {'hilbert','wavelet'};

subjects = unique(manifest.subject, 'stable');
% Restrict to the N=22 analysed set (manifest.mat has 23; sub-005 was
% excluded from the original battery for below-chance accuracy -- see
% s17_v2_checks.py -- and is absent from battery_persubject.csv).
excl_T = readtable(fullfile(fileparts(repo), 'results', 'battery_persubject.csv'));
included_subj = unique(string(excl_T.subject));
subjects = subjects(ismember(subjects, included_subj));
N = numel(subjects);
assert(N == N_SUBJ, 'Expected N=%d subjects, got %d', N_SUBJ, N);

% Cell order (48 = 4 bands x 2 refs x 2 decomps x 3 measures)
cell_defs = {};
for b = 1:N_BANDS
    for d = 1:2, for r = 1:2, for m = 1:3
        cell_defs{end+1} = struct('band',BANDS{b},'measure',MEAS{m}, ...
            'decomposition',DECOMPS{d},'reference',REFS{r}); %#ok<SAGROW>
    end, end, end
end
n_cells = numel(cell_defs);  % should be 48

% Wavelet centre frequencies across all bands (union)
all_freqs = [];
for b = 1:N_BANDS, all_freqs = [all_freqs cfg.battery.wavelet_freqs.(BANDS{b})]; end %#ok<AGROW>
all_freqs = unique(all_freqs);

% Storage: X_sub(su, rep, cell) = ROI-average contrast (high - low)
X_sub = nan(N, N_REPS, n_cells);
ppc_full = nan(N, n_cells);  % PPC on full trials (PLV/hilbert only)

fprintf('Processing %d subjects, %d reps ...\n', N, N_REPS);
for su = 1:N
    subj = char(subjects(su));
    t0 = tic;
    keys = manifest.key(manifest.subject == subjects(su));

    % ---- Load all runs once ------------------------------------------------
    runs = struct('srate',{},'sel_lo',{},'sel_hi',{},'average',{},'csd',{});
    for k = 1:numel(keys)
        key = char(keys(k));
        EEG = pop_loadset('filename',[key '_clean.set'],'filepath',cfg.deriv_dir);
        T_t = load(fullfile(cfg.deriv_dir,[key '_trials.mat'])); trials_k = T_t.trials;
        if k == 1
            have = {EEG.chanlocs.labels};
            labF = have(cellfun(@(l) any(strcmpi(have,l)), cfg.roi_frontal));
            labT = have(cellfun(@(l) any(strcmpi(have,l)), cfg.roi_temporal));
            roiLabs = [labF labT];
            nF = numel(labF); nT = numel(labT);
            fr = 1:nF; te = nF + (1:nT);
        end
        [csd_d, csd_lbs] = local_csd(EEG, key, cfg.fieldtrip_dir);
        runs(k).srate   = EEG.srate;
        [runs(k).sel_lo, runs(k).sel_hi] = local_sel_cond(trials_k, key);
        iAvg = util_chanidx(EEG, roiLabs);
        iCsd = util_chanidx(local_labels2eeg(csd_lbs), roiLabs);
        runs(k).average = double(EEG.data(iAvg,:));
        runs(k).csd     = double(csd_d(iCsd,:));
        clear EEG csd_d
    end
    if isempty(runs), continue; end

    % ---- Precompute per-trial analytic signal segments for all refs/bands ---
    % segs{ref_idx, band_idx, cond_idx}  = [nROI x wlen x n_trials] (Hilbert)
    % wsegs{ref_idx, freq_idx, cond_idx} = [nROI x wlen x n_trials] (wavelet)
    n_freqs = numel(all_freqs);
    segs  = cell(2, N_BANDS, 2);
    wsegs = cell(2, n_freqs,  2);
    n_lo_by_band = zeros(N_BANDS,1);
    n_hi_by_band = zeros(N_BANDS,1);

    for r = 1:numel(REFS)
        rf = REFS{r};
        % Hilbert: one set of segments per band
        for b = 1:N_BANDS
            band = cfg.battery.bands.(BANDS{b});
            [s_lo, s_hi] = local_hil_segs(runs, rf, cfg.maint_window, band, cfg.filt_order_factor);
            segs{r,b,1} = s_lo; segs{r,b,2} = s_hi;
            if r==1, n_lo_by_band(b) = size(s_lo,3); n_hi_by_band(b) = size(s_hi,3); end
        end
        % Wavelet: one set of segments per frequency
        for fi = 1:n_freqs
            fc = all_freqs(fi);
            [ws_lo, ws_hi] = local_wav_segs(runs, rf, cfg.maint_window, fc, n_cycles);
            wsegs{r,fi,1} = ws_lo; wsegs{r,fi,2} = ws_hi;
        end
    end

    % Min trial count across bands (use global min for consistency across cells)
    n_lo_total = size(segs{1,1,1},3);   % total pooled low trials (same for all bands)
    n_hi_total = size(segs{1,1,2},3);
    n_min = min(n_lo_total, n_hi_total);
    if n_min < 5
        warning('%s: min trials=%d < 5, skipping', subj, n_min); continue;
    end

    % ---- PPC on full trials (Hilbert / average ref only) -------------------
    for b = 1:N_BANDS
        for r = 1:numel(REFS)
            rf = REFS{r};
            s_lo = segs{r,b,1}; s_hi = segs{r,b,2};
            n_lo_b = size(s_lo,3); n_hi_b = size(s_hi,3);
            if n_lo_b < 5 || n_hi_b < 5, continue; end
            [plv_lo,~,~,~] = util_conn_pairs(s_lo, fr, te);
            [plv_hi,~,~,~] = util_conn_pairs(s_hi, fr, te);
            ppc_lo = (n_lo_b .* plv_lo.^2 - 1) ./ (n_lo_b - 1);
            ppc_hi = (n_hi_b .* plv_hi.^2 - 1) ./ (n_hi_b - 1);
            cid = local_cell_idx(cell_defs, BANDS{b}, 'PLV', 'hilbert', rf);
            if cid > 0, ppc_full(su, cid) = mean(ppc_hi - ppc_lo); end
        end
    end

    % ---- Subsampled reps ---------------------------------------------------
    for rep = 1:N_REPS
        % Unique reproducible seed per (rep, subject)
        rng(rep * 1000000 + su);
        idx_lo = randperm(n_lo_total, n_min);
        idx_hi = randperm(n_hi_total, n_min);

        for r = 1:numel(REFS)
            rf = REFS{r};

            % Hilbert: subsample from precomputed segments
            for b = 1:N_BANDS
                s_lo = segs{r,b,1}; s_hi = segs{r,b,2};
                nl = size(s_lo,3); nh = size(s_hi,3);
                if nl < n_min || nh < n_min, continue; end
                s_lo_s = s_lo(:,:, idx_lo(1:min(n_min,nl)));
                s_hi_s = s_hi(:,:, idx_hi(1:min(n_min,nh)));
                [plv_lo,dw_lo,imc_lo,~] = util_conn_pairs(s_lo_s, fr, te);
                [plv_hi,dw_hi,imc_hi,~] = util_conn_pairs(s_hi_s, fr, te);
                cid = local_cell_idx(cell_defs, BANDS{b},'PLV','hilbert',rf);
                if cid>0, X_sub(su,rep,cid) = mean(plv_hi - plv_lo); end
                cid = local_cell_idx(cell_defs, BANDS{b},'dwPLI','hilbert',rf);
                if cid>0, X_sub(su,rep,cid) = mean(dw_hi - dw_lo); end
                cid = local_cell_idx(cell_defs, BANDS{b},'ImC','hilbert',rf);
                if cid>0, X_sub(su,rep,cid) = mean(imc_hi - imc_lo); end
            end

            % Wavelet: subsample per-frequency, average within band
            for b = 1:N_BANDS
                band_freqs = cfg.battery.wavelet_freqs.(BANDS{b});
                plv_lo_b = zeros(1, nF*nT); dw_lo_b = zeros(1,nF*nT); imc_lo_b = zeros(1,nF*nT);
                plv_hi_b = zeros(1, nF*nT); dw_hi_b = zeros(1,nF*nT); imc_hi_b = zeros(1,nF*nT);
                n_used = 0;
                for fi = find(ismember(all_freqs, band_freqs))
                    ws_lo = wsegs{r,fi,1}; ws_hi = wsegs{r,fi,2};
                    nl = size(ws_lo,3); nh = size(ws_hi,3);
                    if nl < n_min || nh < n_min, continue; end
                    ws_lo_s = ws_lo(:,:, idx_lo(1:min(n_min,nl)));
                    ws_hi_s = ws_hi(:,:, idx_hi(1:min(n_min,nh)));
                    [p,d,im,~] = util_conn_pairs(ws_lo_s, fr, te);
                    plv_lo_b = plv_lo_b + p; dw_lo_b = dw_lo_b + d; imc_lo_b = imc_lo_b + im;
                    [p,d,im,~] = util_conn_pairs(ws_hi_s, fr, te);
                    plv_hi_b = plv_hi_b + p; dw_hi_b = dw_hi_b + d; imc_hi_b = imc_hi_b + im;
                    n_used = n_used + 1;
                end
                if n_used == 0, continue; end
                cid = local_cell_idx(cell_defs, BANDS{b},'PLV','wavelet',rf);
                if cid>0, X_sub(su,rep,cid) = mean((plv_hi_b - plv_lo_b)/n_used); end
                cid = local_cell_idx(cell_defs, BANDS{b},'dwPLI','wavelet',rf);
                if cid>0, X_sub(su,rep,cid) = mean((dw_hi_b - dw_lo_b)/n_used); end
                cid = local_cell_idx(cell_defs, BANDS{b},'ImC','wavelet',rf);
                if cid>0, X_sub(su,rep,cid) = mean((imc_hi_b - imc_lo_b)/n_used); end
            end
        end
    end
    fprintf('  %s: n_lo=%d, n_hi=%d, n_min=%d (%.1f min)\n', subj, n_lo_total, n_hi_total, n_min, toc(t0)/60);
end

% ── Compute statistics via Python ─────────────────────────────────────────────
tmp_in  = fullfile(repo, 'tmp_A_contrasts.csv');
tmp_out = fullfile(repo, 'tmp_A_stats.csv');
fid = fopen(tmp_in, 'w');
fprintf(fid, 'cell_key,n');
for s = 1:N, fprintf(fid, ',x%d', s); end
fprintf(fid, '\n');
for c = 1:n_cells
    cd_c = cell_defs{c};
    key_str = sprintf('%s__%s__%s__%s', cd_c.band, cd_c.measure, cd_c.decomposition, cd_c.reference);
    for rep = 1:N_REPS
        x = squeeze(X_sub(:, rep, c))';
        fvals = arrayfun(@(v) sprintf('%.10g', v), x, 'UniformOutput', false);
        fprintf(fid, '%s__rep%03d,%d,%s\n', key_str, rep, N, strjoin(fvals, ','));
    end
end
fclose(fid);

py_script = fullfile(repo, 'stats_from_contrasts.py');
cmd = sprintf('python "%s" "%s" "%s"', py_script, tmp_in, tmp_out);
fprintf('Calling Python for stats ...\n');
[status, cmdout] = system(cmd);
if status ~= 0
    warning('Python call failed:\n%s', cmdout);
    fprintf('Writing partial output (no BF01) ...\n');
    % Fall back to dz-only output
    agg_T = local_agg_dz_only(X_sub, N, N_REPS, n_cells, cell_defs);
    writetable(agg_T, fullfile(repo, 'A_trial_equalisation.csv'));
else
    disp(cmdout);
    stats_T = readtable(tmp_out, 'TextType', 'char');
    agg_rows = {};
    for c = 1:n_cells
        cd_c = cell_defs{c};
        key_str = sprintf('%s__%s__%s__%s', cd_c.band, cd_c.measure, cd_c.decomposition, cd_c.reference);
        idx_c = startsWith(string(stats_T.cell_key), string(key_str) + "__rep");
        sub_T = stats_T(idx_c, :);
        if isempty(sub_T), continue; end
        ok = isfinite(sub_T.dz);
        dz_all  = sub_T.dz(ok);
        bf01_all = sub_T.bf01(ok);
        n_rep = sum(ok);
        dz_med  = median(dz_all);  dz_lo = prctile(dz_all,2.5);  dz_hi = prctile(dz_all,97.5);
        bf01_med = median(bf01_all);
        agg_rows{end+1} = {cd_c.band, cd_c.measure, cd_c.decomposition, cd_c.reference, ...
            N, n_rep, dz_med, dz_lo, dz_hi, bf01_med, N_REPS, 1, N_REPS}; %#ok<SAGROW>
    end
    agg_T = cell2table(vertcat(agg_rows{:}), 'VariableNames', ...
        {'band','measure','decomposition','reference','n','n_reps_ok', ...
         'median_dz','ci_lo_2p5','ci_hi_97p5','median_bf01','n_reps_total','seed_start','seed_end'});
    writetable(agg_T, fullfile(repo, 'A_trial_equalisation.csv'));
    if isfile(tmp_in),  delete(tmp_in);  end
    if isfile(tmp_out), delete(tmp_out); end
    fprintf('Wrote A_trial_equalisation.csv (%d rows)\n', height(agg_T));
end

% ── PPC output ────────────────────────────────────────────────────────────────
ppc_rows = {};
for c = 1:n_cells
    cd_c = cell_defs{c};
    if ~strcmp(cd_c.measure,'PLV') || ~strcmp(cd_c.decomposition,'hilbert'), continue; end
    x = ppc_full(:, c);
    x = x(isfinite(x));
    if numel(x) < 3, continue; end
    mn = mean(x); sd_x = std(x); n_pp = numel(x);
    dz_pp = mn / sd_x; t_pp = mn / (sd_x / sqrt(n_pp));
    ppc_rows{end+1} = {cd_c.band, cd_c.reference, n_pp, dz_pp, t_pp}; %#ok<SAGROW>
end
if ~isempty(ppc_rows)
    ppc_T = cell2table(vertcat(ppc_rows{:}), 'VariableNames', ...
        {'band','reference','n','dz_ppc_contrast','t_ppc_contrast'});
    writetable(ppc_T, fullfile(repo, 'A_ppc_full_trials.csv'));
    fprintf('Wrote A_ppc_full_trials.csv (%d rows)\n', height(ppc_T));
end

% ── Markdown ─────────────────────────────────────────────────────────────────
md = {};
md{end+1} = '# A. Trial-count equalisation';
md{end+1} = '';
md{end+1} = sprintf('N = %d subjects; trial counts: low 9–44, high 11–39.', N);
md{end+1} = sprintf('Subsample to min(n_low, n_high) per participant, %d reps (seeds 1–%d).', N_REPS, N_REPS);
md{end+1} = 'Seed per (rep, subject) = rep × 10^6 + subject_index.';
md{end+1} = 'Wavelet: subsampled at trial level per frequency bin, averaged within band.';
md{end+1} = 'PPC: trial-count-unbiased PLV^2 = (n*PLV^2 - 1)/(n-1), full trials, PLV/Hilbert.';
md{end+1} = '';
md{end+1} = 'Outputs:';
md{end+1} = '  A_trial_equalisation.csv: 48 rows; median dz, 2.5/97.5 pct, median BF01.';
md{end+1} = '  A_ppc_full_trials.csv:    PPC x {average, CSD} x 4 bands = 8 rows.';
fid = fopen(fullfile(repo, 'A_trial_equalisation.md'), 'w');
fprintf(fid, '%s\n', md{:});
fclose(fid);
fprintf('Done.\n');

% ===================== local functions =====================
function [s_lo, s_hi] = local_hil_segs(runs, rf, maint_window, band, ord_fac)
    s_lo = []; s_hi = [];
    for k = 1:numel(runs)
        srate = runs(k).srate; win = round(maint_window*srate); wlen = win(2)-win(1)+1;
        Xf = util_bandpass(runs(k).(rf), srate, band, ord_fac);
        Z = hilbert(Xf.').';
        for c = 1:2
            cond_names = {'lo', 'hi'};
            field_name = ['sel_' cond_names{c}];
            sel = runs(k).(field_name);
            tmp = complex(zeros(size(Z,1), wlen, numel(sel))); kept = 0;
            for t = 1:numel(sel)
                a = sel(t).lock_sample + win(1); z = a + wlen - 1;
                if a < 1 || z > size(Z,2), continue; end
                kept = kept + 1; tmp(:,:,kept) = Z(:,a:z);
            end
            if c == 1, s_lo = cat(3, s_lo, tmp(:,:,1:kept));
            else,       s_hi = cat(3, s_hi, tmp(:,:,1:kept)); end
        end
    end
end

function [s_lo, s_hi] = local_wav_segs(runs, rf, maint_window, fc, nc)
    sigma_t = nc / (2*pi*fc);
    s_lo = []; s_hi = [];
    for k = 1:numel(runs)
        srate = runs(k).srate; win = round(maint_window*srate); wlen = win(2)-win(1)+1;
        t_wav = -3*sigma_t : 1/srate : 3*sigma_t;
        wav = exp(2i*pi*fc.*t_wav) .* exp(-t_wav.^2/(2*sigma_t^2));
        wav = wav / sqrt(sum(abs(wav).^2));
        Z = zeros(size(runs(k).(rf)));
        for ch = 1:size(Z,1), Z(ch,:) = conv(runs(k).(rf)(ch,:), wav, 'same'); end
        for c = 1:2
            cond_names = {'lo', 'hi'};
            field_name = ['sel_' cond_names{c}];
            sel = runs(k).(field_name);
            tmp = complex(zeros(size(Z,1), wlen, numel(sel))); kept = 0;
            for t = 1:numel(sel)
                a = sel(t).lock_sample + win(1); z = a + wlen - 1;
                if a < 1 || z > size(Z,2), continue; end
                kept = kept + 1; tmp(:,:,kept) = Z(:,a:z);
            end
            if c == 1, s_lo = cat(3, s_lo, tmp(:,:,1:kept));
            else,       s_hi = cat(3, s_hi, tmp(:,:,1:kept)); end
        end
    end
end

function cid = local_cell_idx(cell_defs, band, measure, decomp, ref)
    cid = 0;
    for i = 1:numel(cell_defs)
        cd = cell_defs{i};
        if strcmp(cd.band,band) && strcmp(cd.measure,measure) && ...
           strcmp(cd.decomposition,decomp) && strcmp(cd.reference,ref)
            cid = i; return;
        end
    end
end

function [csd_data, csd_labels] = local_csd(EEG, key, ft_dir)
    addpath(ft_dir);
    noPos = arrayfun(@(c) isempty(c.X)||isempty(c.Y)||isempty(c.Z)||...
        any(~isfinite(double([c.X c.Y c.Z]))), EEG.chanlocs);
    if any(noPos)
        error('A_trial_equalisation:no_positions','%s: CSD requires 3-D positions',key);
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

function T = local_agg_dz_only(X_sub, N, N_REPS, n_cells, cell_defs)
    rows = {};
    for c = 1:n_cells
        cd_c = cell_defs{c};
        dz_reps = nan(N_REPS,1);
        for rep = 1:N_REPS
            x = squeeze(X_sub(:,rep,c)); x = x(isfinite(x));
            if numel(x)>=3 && std(x)>0, dz_reps(rep) = mean(x)/std(x); end
        end
        ok = isfinite(dz_reps);
        rows{end+1} = {cd_c.band, cd_c.measure, cd_c.decomposition, cd_c.reference, ...
            N, sum(ok), median(dz_reps(ok)), prctile(dz_reps(ok),2.5), ...
            prctile(dz_reps(ok),97.5), NaN, N_REPS, 1, N_REPS}; %#ok<SAGROW>
    end
    T = cell2table(vertcat(rows{:}), 'VariableNames', ...
        {'band','measure','decomposition','reference','n','n_reps_ok', ...
         'median_dz','ci_lo_2p5','ci_hi_97p5','median_bf01','n_reps_total','seed_start','seed_end'});
end
