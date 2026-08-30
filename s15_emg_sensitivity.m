% s15_emg_sensitivity.m -- Reviewer-requested check (B2), PARTIAL:
% Quantifies residual broadband power at temporal vs. frontal ROI electrodes
% and correlates it with frontotemporal beta dwPLI, as an index of whether
% unequal myogenic (muscle) noise exposure could be driving the
% power-connectivity dissociation, rather than genuine physiology.
%
% ============================ IMPORTANT CAVEATS ============================
% 1) RANGE: the reviewer requests residual power in 30-100 Hz. The
%    preprocessing pipeline (s01_preprocess.m) applies a 45 Hz low-pass
%    filter (cfg.lp_hz = 45), so frequencies above 45 Hz do NOT exist in the
%    _clean.set files. This script therefore uses 30-45 Hz as a TRUNCATED
%    proxy for broadband EMG contamination. This is a genuine limitation of
%    the check, not a coding choice -- state it explicitly if these numbers
%    are reported, since most EMG power lies above 45 Hz and this proxy will
%    understate true myogenic contamination.
% 2) NOT COVERED HERE: "mean number of muscle ICs removed per subject."
%    pop_subcomp() in s01_preprocess.m discards which components were
%    rejected once it removes them, so this information no longer exists in
%    the saved _clean.set files. Recovering it requires RE-RUNNING the full
%    preprocessing pipeline (ICA on every dataset, likely hours) with one
%    added logging line before pop_subcomp is called. Not done in this
%    script -- say the word if you want that version of s01 written; it's a
%    rerun, not a quick add-on.
% =============================================================================

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;
subjects = unique(manifest.subject, 'stable');
emg_band = [30 45];   % truncated proxy; see caveat (1) above

rows = {};
for su = 1:numel(subjects)
    subj = subjects(su);
    keys = manifest.key(manifest.subject == subj);
    fr_pow = []; te_pow = [];
    refEEG = [];

    for k = 1:numel(keys)
        key = char(keys(k));
        EEG = pop_loadset('filename',[key '_clean.set'],'filepath',cfg.deriv_dir);
        if isempty(refEEG), refEEG = EEG; end
        Xf = util_bandpass(double(EEG.data), EEG.srate, emg_band, cfg.filt_order_factor);
        fr = util_chanidx(EEG, cfg.roi_frontal);
        te = util_chanidx(EEG, cfg.roi_temporal);
        fr_pow = [fr_pow, mean(Xf(fr,:).^2, 2)']; %#ok<AGROW>
        te_pow = [te_pow, mean(Xf(te,:).^2, 2)']; %#ok<AGROW>
    end
    if isempty(refEEG), continue; end
    rows(end+1,:) = {char(subj), mean(fr_pow), mean(te_pow), mean(te_pow)/mean(fr_pow)}; %#ok<SAGROW>
    fprintf('done %s\n', subj);
end

R = cell2table(rows, 'VariableNames', {'subject','frontal_3045_power','temporal_3045_power','temporal_frontal_ratio'});
writetable(R, fullfile(cfg.deriv_dir,'emg_band_power.csv'));

% ---- correlate temporal EMG-band power with frontotemporal dwPLI (beta, high-low diff) ----
Cn = load(fullfile(cfg.deriv_dir,'conn_results.mat'), 'R'); Rconn = Cn.R;
subs = R.subject;
d_cn = nan(numel(subs),1);
for i = 1:numel(subs)
    rl = Rconn(strcmp(Rconn.subject,subs{i}) & strcmp(Rconn.band,'beta') & strcmp(Rconn.load,'low'), :);
    rh = Rconn(strcmp(Rconn.subject,subs{i}) & strcmp(Rconn.band,'beta') & strcmp(Rconn.load,'high'), :);
    if ~isempty(rl) && ~isempty(rh), d_cn(i) = rh.dwpli(1) - rl.dwpli(1); end
end
ok = ~isnan(d_cn);
[rho, p_rho] = corr(R.temporal_3045_power(ok), d_cn(ok), 'Type','Pearson');

fprintf('\nTemporal EMG-band (30-45Hz) power vs beta dwPLI load-diff, N=%d\n', sum(ok));
fprintf('  Pearson r = %.3f, p = %.4f\n', rho, p_rho);
fprintf('  Mean temporal/frontal 30-45Hz power ratio = %.2f\n', mean(R.temporal_frontal_ratio));

save(fullfile(cfg.deriv_dir,'emg_sensitivity_results.mat'), 'R','rho','p_rho');
