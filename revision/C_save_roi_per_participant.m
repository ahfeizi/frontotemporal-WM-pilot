% C_save_roi_per_participant.m
% Compute the ROI-average (high - low) connectivity contrast per participant
% for all 48 cells (4 bands x 3 measures x 2 decompositions x 2 references).
%
% Per-subject contrasts are computed from the preprocessed EEG files
% directly, not from battery_persubject.csv (which may include subjects
% from other datasets and may not contain all 48 cells).
%
% Writes: revision/roi_per_participant.csv  (N x 48 cells -> one row per
%         participant x cell; used by C_multiverse_summary.py bootstrap)
%         revision/C_save_roi_per_participant.md
%
% Reads: <deriv_dir>/<key>_clean.set + <key>_trials.mat (per run)

clear; clc;
repo = fileparts(mfilename('fullpath'));  % revision/
addpath(fileparts(repo));                % repo root contains config.m
cfg  = config();

addpath(cfg.eeglab_dir); eeglab nogui;
addpath(cfg.fieldtrip_dir); ft_defaults;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

n_cycles = cfg.battery.n_cycles;
REFS   = {'average','csd'};
BANDS  = fieldnames(cfg.battery.bands);
N_BANDS = numel(BANDS);
MEAS   = {'PLV','dwPLI','ImC'};
DECOMPS = {'hilbert','wavelet'};

all_freqs = [];
for b = 1:N_BANDS, all_freqs = [all_freqs cfg.battery.wavelet_freqs.(BANDS{b})]; end %#ok<AGROW>
all_freqs = unique(all_freqs);

% Cell ordering: 48 cells
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

% ROI average per subject x cell
roi_contrast = nan(N, n_cells);

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

        % Hilbert
        for b = 1:N_BANDS
            band = cfg.battery.bands.(BANDS{b});
            s_lo = local_hil_segs(runs, rf, cfg.maint_window, band, cfg.filt_order_factor, 'lo');
            s_hi = local_hil_segs(runs, rf, cfg.maint_window, band, cfg.filt_order_factor, 'hi');
            if size(s_lo,3) < 5 || size(s_hi,3) < 5, continue; end
            [plv_lo,dw_lo,imc_lo,~] = util_conn_pairs(s_lo, fr, te);
            [plv_hi,dw_hi,imc_hi,~] = util_conn_pairs(s_hi, fr, te);
            plv_c = mean(plv_hi - plv_lo); dw_c = mean(dw_hi - dw_lo); imc_c = mean(imc_hi - imc_lo);
            for m = 1:3
                cid = local_cid(cell_defs, BANDS{b}, MEAS{m}, 'hilbert', rf);
                if cid == 0, continue; end
                switch m
                    case 1, roi_contrast(su, cid) = plv_c;
                    case 2, roi_contrast(su, cid) = dw_c;
                    case 3, roi_contrast(su, cid) = imc_c;
                end
            end
        end

        % Wavelet: per-frequency then average within band
        nf = numel(all_freqs);
        plv_lo_f = nan(nf, nF*nT); plv_hi_f = nan(nf, nF*nT);
        dw_lo_f  = nan(nf, nF*nT); dw_hi_f  = nan(nf, nF*nT);
        imc_lo_f = nan(nf, nF*nT); imc_hi_f = nan(nf, nF*nT);
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
                    case 1, roi_contrast(su,cid) = mean(mean(plv_hi_f(used,:),1) - mean(plv_lo_f(used,:),1));
                    case 2, roi_contrast(su,cid) = mean(mean(dw_hi_f(used,:),1)  - mean(dw_lo_f(used,:),1));
                    case 3, roi_contrast(su,cid) = mean(mean(imc_hi_f(used,:),1) - mean(imc_lo_f(used,:),1));
                end
            end
        end
    end
    fprintf('  %s done (%.1f min)\n', subj, toc(t0)/60);
end

% Build long-format output table
out_rows = {};
for su = 1:N
    subj = char(subjects(su));
    for c = 1:n_cells
        cd = cell_defs{c};
        out_rows{end+1} = {subj, cd.band, cd.measure, cd.decomposition, cd.reference, roi_contrast(su,c)}; %#ok<SAGROW>
    end
end
T = cell2table(vertcat(out_rows{:}), ...
    'VariableNames', {'participant','band','measure','decomposition','reference','roi_value'});

fprintf('Wrote %d rows (%d subjects x %d cells).\n', height(T), N, n_cells);
writetable(T, fullfile(repo, 'roi_per_participant.csv'));

% Markdown
md_lines = {
    '# C_save_roi_per_participant'
    ''
    sprintf('Subjects: %d (from manifest.mat)', N)
    sprintf('Cells: %d', n_cells)
    sprintf('Output rows: %d', height(T))
    ''
    'Each row: participant, band, measure, decomposition, reference, roi_value (high - low, ROI average).'
    'roi_value = mean over ROI pairs of (high_condition_value - low_condition_value).'
    'Per-subject contrasts computed from EEG files (not battery_persubject.csv).'
    ''
    'Output: roi_per_participant.csv'
};
fid = fopen(fullfile(repo, 'C_save_roi_per_participant.md'), 'w');
fprintf(fid, '%s\n', md_lines{:});
fclose(fid);
fprintf('Done.\n');

% ===================== local functions =====================
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

function [csd_data, csd_labels] = local_csd(EEG, key, ft_dir)
    addpath(ft_dir);
    noPos = arrayfun(@(c) isempty(c.X)||isempty(c.Y)||isempty(c.Z)||...
        any(~isfinite(double([c.X c.Y c.Z]))), EEG.chanlocs);
    if any(noPos)
        error('C_save_roi:no_positions','%s: CSD requires 3-D positions', key);
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
