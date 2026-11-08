% H_wavelet_cycles.m
% Wavelet cycles sensitivity: rerun the wavelet portion of the battery with
% n_cycles = 4, 5 (existing), and 7, for theta and beta bands only.
% Hilbert results are unchanged; CSD and average reference are both included.
%
% Reports per-cell (ROI-average) dz, 95% CI, t, parametric p, BF01 for each
% cycle count.  Full stats via Python stats_from_contrasts.py.
%
% Reads:  <deriv_dir>/<key>_clean.set + <key>_trials.mat
% Writes: revision/H_wavelet_cycles.csv  (one row per [cell x n_cycles]: 24 rows)
%         revision/H_wavelet_cycles.md

clear; clc;
repo = fileparts(mfilename('fullpath'));
addpath(fileparts(repo));
cfg  = config();

addpath(cfg.eeglab_dir); eeglab nogui;
addpath(cfg.fieldtrip_dir); ft_defaults;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

CYCLES_LIST = [4, 5, 7];   % cfg.battery.n_cycles is 5; test 4 and 7 as well
BANDS = {'theta','beta'};   % theta 4-7 Hz, beta 13-30 Hz
REFS  = {'average','csd'};
CONDS = {'low','high'};

subjects = unique(manifest.subject, 'stable');
% Restrict to the N=22 analysed set (manifest.mat has 23; sub-005 was
% excluded from the original battery -- see A_trial_equalisation.m).
excl_T = readtable(fullfile(fileparts(repo), 'results', 'battery_persubject.csv'));
included_subj = unique(string(excl_T.subject));
subjects = subjects(ismember(subjects, included_subj));
N = numel(subjects);

% Contrast arrays: [N x n_cycles] per cell (band x ref)
n_cells_hw = numel(BANDS) * numel(REFS) * 3;  % x3 for PLV/dwPLI/ImC

% Storage: band x measure x ref x n_cycles -> N per-subj contrasts
MEAS = {'PLV','dwPLI','ImC'};
X_all = nan(N, numel(CYCLES_LIST), numel(BANDS), numel(MEAS), numel(REFS));

for su = 1:N
    subj = char(subjects(su));
    t0 = tic;
    keys = manifest.key(manifest.subject == subjects(su));
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
        [csd_data, csd_labels] = local_csd(EEG, key, cfg.fieldtrip_dir);
        iAvg = util_chanidx(EEG, roiLabs);
        iCsd = util_chanidx(local_labels2eeg(csd_labels), roiLabs);
        runs(k).srate   = EEG.srate;
        [runs(k).sel_lo, runs(k).sel_hi] = local_sel_cond(trials_k, key);
        runs(k).average = double(EEG.data(iAvg,:));
        runs(k).csd     = double(csd_data(iCsd,:));
        clear EEG csd_data
    end
    if isempty(runs), continue; end

    for r = 1:numel(REFS)
        rf = REFS{r};
        for b = 1:numel(BANDS)
            bn = BANDS{b};
            band = cfg.battery.bands.(bn);
            wavfreqs = cfg.battery.wavelet_freqs.(bn);
            for nc_idx = 1:numel(CYCLES_LIST)
                nc = CYCLES_LIST(nc_idx);
                % Compute per-frequency wavelet connectivity then average within band
                nf = numel(wavfreqs);
                plv_lo  = nan(nf, nF*nT); plv_hi  = nan(nf, nF*nT);
                dw_lo   = nan(nf, nF*nT); dw_hi   = nan(nf, nF*nT);
                imc_lo  = nan(nf, nF*nT); imc_hi  = nan(nf, nF*nT);
                ok_lo   = false(1,nf); ok_hi = false(1,nf);
                for fi = 1:nf
                    fc = wavfreqs(fi);
                    seg_lo = local_pool_runs(runs, rf, cfg.maint_window, ...
                        @(X,sr) local_morlet(X, sr, fc, nc), 'lo');
                    seg_hi = local_pool_runs(runs, rf, cfg.maint_window, ...
                        @(X,sr) local_morlet(X, sr, fc, nc), 'hi');
                    if size(seg_lo,3) >= 5
                        [plv_lo(fi,:), dw_lo(fi,:), imc_lo(fi,:), ~] = util_conn_pairs(seg_lo, fr, te);
                        ok_lo(fi) = true;
                    end
                    if size(seg_hi,3) >= 5
                        [plv_hi(fi,:), dw_hi(fi,:), imc_hi(fi,:), ~] = util_conn_pairs(seg_hi, fr, te);
                        ok_hi(fi) = true;
                    end
                end
                used = find(ok_lo & ok_hi);
                nu = numel(used);
                if nu == 0, continue; end
                plv_c   = mean(plv_hi(used,:),1)  - mean(plv_lo(used,:),1);
                dw_c    = mean(dw_hi(used,:),1)   - mean(dw_lo(used,:),1);
                imc_c   = mean(imc_hi(used,:),1)  - mean(imc_lo(used,:),1);
                X_all(su, nc_idx, b, 1, r) = mean(plv_c);
                X_all(su, nc_idx, b, 2, r) = mean(dw_c);
                X_all(su, nc_idx, b, 3, r) = mean(imc_c);
            end
        end
    end
    fprintf('  %s done (%.1f min)\n', subj, toc(t0)/60);
end

% Write contrasts to temp CSV and call Python for stats
tmp_in  = fullfile(repo, 'tmp_H_contrasts.csv');
tmp_out = fullfile(repo, 'tmp_H_stats.csv');
fid = fopen(tmp_in, 'w');
fprintf(fid, 'cell_key,n');
for s = 1:N, fprintf(fid, ',x%d', s); end
fprintf(fid, '\n');
for b = 1:numel(BANDS)
    for m = 1:numel(MEAS)
        for r = 1:numel(REFS)
            for nc_idx = 1:numel(CYCLES_LIST)
                x = X_all(:, nc_idx, b, m, r);
                key_str = sprintf('%s__%s__wavelet__%s__ncyc%d', ...
                    BANDS{b}, MEAS{m}, REFS{r}, CYCLES_LIST(nc_idx));
                fvals = arrayfun(@(v) sprintf('%.10g', v), x', 'UniformOutput', false);
                fprintf(fid, '%s,%d,%s\n', key_str, N, strjoin(fvals, ','));
            end
        end
    end
end
fclose(fid);

py_script = fullfile(repo, 'stats_from_contrasts.py');
cmd = sprintf('python "%s" "%s" "%s"', py_script, tmp_in, tmp_out);
fprintf('Calling Python: %s\n', cmd);
[status, cmdout] = system(cmd);
if status ~= 0
    warning('Python call failed:\n%s', cmdout);
else
    disp(cmdout);
    stats_T = readtable(tmp_out, 'TextType', 'char');
    % Parse key to extract metadata
    parts_all = cellfun(@(k) strsplit(k,'__'), stats_T.cell_key, 'UniformOutput',false);
    n_rows = height(stats_T);
    out_band  = cell(n_rows,1); out_meas = cell(n_rows,1);
    out_ref   = cell(n_rows,1); out_ncyc = nan(n_rows,1);
    for i = 1:n_rows
        p = parts_all{i};
        out_band{i}  = p{1}; out_meas{i} = p{2};
        out_ref{i}   = p{4};
        out_ncyc(i)  = str2double(strrep(p{5},'ncyc',''));
    end
    stats_T.band          = out_band;
    stats_T.measure       = out_meas;
    stats_T.decomposition = repmat({'wavelet'}, n_rows, 1);
    stats_T.reference     = out_ref;
    stats_T.n_cycles      = out_ncyc;
    stats_T = removevars(stats_T, 'cell_key');
    writetable(stats_T, fullfile(repo, 'H_wavelet_cycles.csv'));
    delete(tmp_in); delete(tmp_out);
    fprintf('Wrote H_wavelet_cycles.csv (%d rows)\n', height(stats_T));
end

% Markdown
md = {};
md{end+1} = '# H. Wavelet cycles sensitivity (theta and beta)';
md{end+1} = '';
md{end+1} = sprintf('N = %d; wavelet cycles: 4, 5 (default), 7.', N);
md{end+1} = 'Hilbert results unchanged. Bands: theta (4-7 Hz), beta (13-30 Hz).';
md{end+1} = sprintf('Cells: 2 bands x 3 measures x 2 refs x 3 cycle counts = %d rows.', ...
    numel(BANDS)*numel(MEAS)*numel(REFS)*numel(CYCLES_LIST));
md{end+1} = '';
md{end+1} = 'Output: H_wavelet_cycles.csv (dz, 95% CI, t, p, BF01 per row).';
fid = fopen(fullfile(repo, 'H_wavelet_cycles.md'), 'w');
fprintf(fid, '%s\n', md{:});
fclose(fid);
fprintf('Done.\n');

% ===================== local functions =====================
function [csd_data, csd_labels] = local_csd(EEG, key, ft_dir)
    addpath(ft_dir);
    noPos = arrayfun(@(c) isempty(c.X)||isempty(c.Y)||isempty(c.Z)||...
        any(~isfinite(double([c.X c.Y c.Z]))), EEG.chanlocs);
    if any(noPos)
        error('H_wavelet_cycles:no_positions','%s: %d channel(s) lack positions',key,sum(noPos));
    end
    ft_raw = eeglab2fieldtrip(EEG,'raw','none');
    cfg_csd = []; cfg_csd.method = 'spline';
    ft_csd = ft_scalpcurrentdensity(cfg_csd, ft_raw);
    csd_data = ft_csd.trial{1}; csd_labels = ft_csd.label;
end

function Z = local_morlet(data, srate, fc, n_cycles)
    sigma_t = n_cycles / (2*pi*fc);
    t = -3*sigma_t : 1/srate : 3*sigma_t;
    wavelet = exp(2i*pi*fc.*t) .* exp(-t.^2/(2*sigma_t^2));
    wavelet = wavelet / sqrt(sum(abs(wavelet).^2));
    Z = zeros(size(data));
    for ch = 1:size(data,1), Z(ch,:) = conv(data(ch,:), wavelet,'same'); end
end

function seg = local_pool_runs(runs, rf, maint_window, transform, cond_str)
    seg = [];
    for k = 1:numel(runs)
        srate = runs(k).srate; win = round(maint_window*srate); wlen = win(2)-win(1)+1;
        Z = transform(runs(k).(rf), srate);
        sel = runs(k).(['sel_' cond_str]);
        tmp = complex(zeros(size(Z,1),wlen,numel(sel))); kept = 0;
        for t = 1:numel(sel)
            a = sel(t).lock_sample + win(1); z = a + wlen - 1;
            if a < 1 || z > size(Z,2), continue; end
            kept = kept + 1; tmp(:,:,kept) = Z(:,a:z);
        end
        seg = cat(3, seg, tmp(:,:,1:kept));
    end
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
