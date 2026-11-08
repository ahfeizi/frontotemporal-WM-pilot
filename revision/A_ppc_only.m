% A_ppc_only.m
% Runs only the per-subject load + PPC loop from A_trial_equalisation.m.
% Skips the 200-rep subsampling. Saves ppc_full to A_workspace.mat, then
% writes A_ppc_full_trials.csv and A_trial_equalisation.md.
%
% Expected run time: ~5-10 min (47 subjects, Hilbert only, no subsampling).
%
% Reads:  <deriv_dir>/<key>_clean.set + <key>_trials.mat (per run)
% Writes: revision/A_workspace.mat         (ppc_full, cell_defs for recovery)
%         revision/A_ppc_full_trials.csv   (PPC, PLV/Hilbert/average|csd)
%         revision/A_trial_equalisation.md (static summary)

clear; clc;
repo = fileparts(mfilename('fullpath'));
addpath(fileparts(repo));
cfg  = config();

addpath(cfg.eeglab_dir); eeglab nogui;
addpath(cfg.fieldtrip_dir); ft_defaults;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

N_REPS    = 200;   % kept for the .md output only
N_SUBJ    = 22;
REFS      = {'average','csd'};
BANDS     = fieldnames(cfg.battery.bands);
N_BANDS   = numel(BANDS);
MEAS      = {'PLV','dwPLI','ImC'};
DECOMPS   = {'hilbert','wavelet'};

subjects = unique(manifest.subject, 'stable');
% Restrict to the N=22 analysed set (manifest.mat has 23; sub-005 was
% excluded from the original battery -- see A_trial_equalisation.m).
excl_T = readtable(fullfile(fileparts(repo), 'results', 'battery_persubject.csv'));
included_subj = unique(string(excl_T.subject));
subjects = subjects(ismember(subjects, included_subj));
N = numel(subjects);
assert(N == N_SUBJ, 'Expected N=%d subjects, got %d', N_SUBJ, N);

% Cell order identical to A_trial_equalisation.m (48 = 4 bands x 2 refs x 2 decomps x 3 measures)
cell_defs = {};
for b = 1:N_BANDS
    for d = 1:2, for r = 1:2, for m = 1:3
        cell_defs{end+1} = struct('band',BANDS{b},'measure',MEAS{m}, ...
            'decomposition',DECOMPS{d},'reference',REFS{r}); %#ok<SAGROW>
    end, end, end
end
n_cells = numel(cell_defs);  % 48

ppc_full = nan(N, n_cells);  % PPC on full trials (PLV/hilbert only)

fprintf('Processing %d subjects (PPC only, no subsampling) ...\n', N);
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

    % ---- Precompute Hilbert segments for all refs/bands --------------------
    segs = cell(2, N_BANDS, 2);
    for r = 1:numel(REFS)
        rf = REFS{r};
        for b = 1:N_BANDS
            band = cfg.battery.bands.(BANDS{b});
            [s_lo, s_hi] = local_hil_segs(runs, rf, cfg.maint_window, band, cfg.filt_order_factor);
            segs{r,b,1} = s_lo; segs{r,b,2} = s_hi;
        end
    end

    % Sanity check: skip if global trial count too low
    n_lo_total = size(segs{1,1,1},3);
    n_hi_total = size(segs{1,1,2},3);
    n_min = min(n_lo_total, n_hi_total);
    if n_min < 5
        warning('%s: min trials=%d < 5, skipping', subj, n_min); continue;
    end

    % ---- PPC on full trials (Hilbert / average and csd ref) ----------------
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
    fprintf('  %s: n_lo=%d, n_hi=%d (%.1f min)\n', subj, n_lo_total, n_hi_total, toc(t0)/60);
end

% Save workspace so A_recover_missing_outputs.m can load ppc_full
save(fullfile(repo,'A_workspace.mat'), 'ppc_full', 'cell_defs');
fprintf('Saved A_workspace.mat\n');

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
        cue_vals = [trials.cue];
        sel_lo = trials([trials.good] & cue_vals == 1);
        sel_hi = trials([trials.good] & cue_vals == 2);
    else
        fn = fieldnames(trials);
        error('local_sel_cond: no condition field. key=%s; fields: %s', key, strjoin(fn, ', '));
    end
end
