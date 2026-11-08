% G_cell_discrepancy.m
% Reviewer cites t = -2.23, exact p = 0.017 vs parametric p = 0.037 (§4.1):
% FCz beta-band power at maintenance onset, high vs low load.
%
% This script:
%  (1) Searches unified_contrasts_ds004117.csv for ROI cells where
%      |exact_p - parametric_p| is largest (connectivity, all 48 cells).
%  (2) Extracts per-participant FCz beta power contrast (high - low) from
%      the preprocessed EEG files (same window as main analysis).
%  (3) Reports: skewness, Shapiro-Wilk W and p, leave-one-out dz vector,
%      and the participant (sub-024) leave-one-out effect.
%
% Reads:  results/unified_contrasts_ds004117.csv (ROI rows, connectivity)
%         <deriv_dir>/<key>_clean.set + <key>_trials.mat (power extraction)
% Writes: revision/G_cell_discrepancy.csv  (per-participant FCz beta contrast)
%         revision/G_cell_discrepancy.md

clear; clc;
repo = fileparts(mfilename('fullpath'));
addpath(fileparts(repo));
cfg  = config();
rng(1);

% ── (1) Connectivity divergence check ────────────────────────────────────────
csv_in = fullfile(fileparts(repo), 'results', 'unified_contrasts_ds004117.csv');
opts = detectImportOptions(csv_in, 'Delimiter', ',');
opts = setvartype(opts, {'dataset','contrast','level','band','measure','decomposition', ...
    'reference','pair','perm_method'}, 'char');
opts = setvartype(opts, {'n','dz','dz_ci_lo','dz_ci_hi','t','df','p','p_perm', ...
    'p_fwer','bf10','bf01','p_tost'}, 'double');
T = readtable(csv_in, opts);
roi_T = T(strcmp(T.level, 'roi'), :);
roi_T.p_exact = roi_T.p_perm;
roi_T.abs_divergence = abs(roi_T.p - roi_T.p_exact);
[~, idx] = sort(roi_T.abs_divergence, 'descend');
top5 = roi_T(idx(1:min(5, height(roi_T))), ...
    {'band','measure','decomposition','reference','t','p','p_exact','abs_divergence','bf01'});
fprintf('Top 5 connectivity cells by |exact_p - parametric_p|:\n');
disp(top5);

% ── (2) FCz beta power per participant ───────────────────────────────────────
addpath(cfg.eeglab_dir); eeglab nogui;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;
subjects = unique(manifest.subject, 'stable');
% Restrict to the N=22 analysed set (manifest.mat has 23; sub-005 was
% excluded from the original battery -- see A_trial_equalisation.m).
excl_T = readtable(fullfile(fileparts(repo), 'results', 'battery_persubject.csv'));
included_subj = unique(string(excl_T.subject));
subjects = subjects(ismember(subjects, included_subj));
N = numel(subjects);

band_beta = [13 30];
fcz_label = 'FCz';
win_s = cfg.maint_window;     % [0.5 2.5] s
conds = {'low','high'};

pow_subj = nan(N, 2);  % N x 2 (low, high)
sub_labels = cell(N, 1);

for s = 1:N
    subj = char(subjects(s));
    sub_labels{s} = subj;
    keys = manifest.key(manifest.subject == subjects(s));
    all_lo = []; all_hi = [];
    for k = 1:numel(keys)
        key = char(keys(k));
        EEG = pop_loadset('filename', [key '_clean.set'], 'filepath', cfg.deriv_dir);
        T_t = load(fullfile(cfg.deriv_dir, [key '_trials.mat'])); trials = T_t.trials;
        iFCz = util_chanidx(EEG, {fcz_label});
        if isempty(iFCz)
            warning('%s: FCz not found, skipping run', key);
            continue;
        end
        srate = EEG.srate;
        win_samp = round(win_s * srate);   % [125 625] at 250 Hz
        wlen = win_samp(2) - win_samp(1) + 1;
        % Band-pass filter
        raw_fcz = double(EEG.data(iFCz, :));
        Xf = util_bandpass(raw_fcz, srate, band_beta, cfg.filt_order_factor);
        [sel_lo_t, sel_hi_t] = local_sel_cond(trials, key);
        for c = 1:2
            cn = conds{c};
            sel = sel_lo_t; if c == 2, sel = sel_hi_t; end
            for tr = 1:numel(sel)
                a = sel(tr).lock_sample + win_samp(1);
                z = a + wlen - 1;
                if a < 1 || z > numel(Xf), continue; end
                seg = Xf(a:z);
                % Log power: mean(log(bandpass_envelope^2)) as used elsewhere
                % envelope via Hilbert; consistent with cfg.log_power = true
                env2 = abs(hilbert(seg)).^2;
                if strcmpi(cn, 'low')
                    all_lo(end+1) = mean(log(env2 + eps)); %#ok<AGROW>
                else
                    all_hi(end+1) = mean(log(env2 + eps)); %#ok<AGROW>
                end
            end
        end
        clear EEG
    end
    pow_subj(s, 1) = mean(all_lo, 'omitnan');
    pow_subj(s, 2) = mean(all_hi, 'omitnan');
end

diff_pow = pow_subj(:, 2) - pow_subj(:, 1);  % high - low
ok = isfinite(diff_pow);
x = diff_pow(ok);
n_ok = sum(ok);

% Basic stats (same formulas as util_contrast_stats.py)
mn = mean(x); sd_x = std(x); t_stat = mn / (sd_x / sqrt(n_ok));
p_param = 2 * tcdf(-abs(t_stat), n_ok - 1);
dz = mn / sd_x;

% Skewness and Shapiro-Wilk
sk = skewness(x, 1);
[sw_h, sw_p, sw_w] = swtest(x);   % MATLAB Statistics toolbox

% Leave-one-out dz
loo_dz = nan(n_ok, 1);
x_idx = find(ok);
for i = 1:n_ok
    xi = x; xi(i) = [];
    loo_dz(i) = mean(xi) / std(xi);
end

% Role of sub-024
sub024_idx = find(cellfun(@(s) contains(s,'024'), sub_labels(ok)));
if ~isempty(sub024_idx)
    x_no024 = x; x_no024(sub024_idx) = [];
    t_no024 = mean(x_no024) / (std(x_no024) / sqrt(numel(x_no024)));
    dz_no024 = mean(x_no024) / std(x_no024);
else
    t_no024 = NaN; dz_no024 = NaN;
    fprintf('sub-024 not found among included subjects.\n');
end

fprintf('\nFCz beta power contrast (high - low), N = %d:\n', n_ok);
fprintf('  dz = %.4f,  t = %.4f,  p_param = %.4f\n', dz, t_stat, p_param);
fprintf('  skewness = %.4f,  Shapiro-Wilk W = %.4f,  p_SW = %.4f\n', sk, sw_w, sw_p);
fprintf('  dz without sub-024 = %.4f,  t = %.4f\n', dz_no024, t_no024);

% Write per-participant CSV
out_T = table(sub_labels(ok), x, loo_dz, ...
    'VariableNames', {'participant', 'contrast_high_minus_low', 'loo_dz'});
out_csv = fullfile(repo, 'G_cell_discrepancy.csv');
writetable(out_T, out_csv);
fprintf('Wrote %s\n', out_csv);

% Markdown
md = {};
md{end+1} = '# G. FCz beta power: exact vs parametric p discrepancy';
md{end+1} = '';
md{end+1} = '## Connectivity cells with largest |exact p − parametric p| (top 5)';
for i = 1:height(top5)
    md{end+1} = sprintf('- %s %s %s %s: t=%.4f, p_param=%.4f, p_exact=%.4f, divergence=%.4f', ...
        top5.band{i}, top5.measure{i}, top5.decomposition{i}, top5.reference{i}, ...
        top5.t(i), top5.p(i), top5.p_exact(i), top5.abs_divergence(i));
end
md{end+1} = '';
md{end+1} = '## FCz beta power (manuscript §4.1: t=-2.23, exact p=0.017, param p=0.037)';
md{end+1} = sprintf('N included = %d', n_ok);
md{end+1} = sprintf('Computed: dz = %.4f, t = %.4f, p_param = %.4f', dz, t_stat, p_param);
md{end+1} = sprintf('Skewness = %.4f', sk);
md{end+1} = sprintf('Shapiro-Wilk W = %.4f, p = %.4f (reject normality: %d)', sw_w, sw_p, sw_h);
md{end+1} = sprintf('LOO dz range = [%.4f, %.4f]', min(loo_dz), max(loo_dz));
md{end+1} = sprintf('Without sub-024: dz = %.4f, t = %.4f', dz_no024, t_no024);
md{end+1} = '';
md{end+1} = sprintf('Output: G_cell_discrepancy.csv (%d rows: participant, contrast, loo_dz)', n_ok);
fid = fopen(fullfile(repo, 'G_cell_discrepancy.md'), 'w');
fprintf(fid, '%s\n', md{:});
fclose(fid);
fprintf('Done.\n');

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
