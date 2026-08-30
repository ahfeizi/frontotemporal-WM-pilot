% s14_amplitude_phase_confound.m -- Reviewer-requested checks (C2):
% Part 1: per-subject correlation between the beta power load-difference and
%         the beta dwPLI load-difference. A significant POSITIVE correlation
%         would support the account that amplitude-driven phase-estimation
%         noise (not a genuine physiological null) is behind the beta
%         connectivity equivalence result.
% Part 2: sensitivity check -- recompute beta dwPLI/PLV with a per-trial
%         amplitude normalization (z-score the bandpassed real signal within
%         each trial window) applied BEFORE the Hilbert transform, then
%         compare against the primary (unnormalized) result.
%
% CAVEAT on Part 2: normalization is applied per trial WINDOW rather than to
% the continuous signal, so the Hilbert transform here runs on each short
% windowed segment individually. This is more vulnerable to edge effects
% than the continuous-signal approach in s03_connectivity.m. Acceptable for
% a sensitivity check, but flag this if the numbers are reported.

clear; clc;
cfg = config();

% ---------------- Part 1: power-change vs connectivity-change correlation ----------------
load(fullfile(cfg.deriv_dir,'conn_results.mat'), 'R');
band = 'beta';
subs = unique(R.subject);
pw_lo = nan(numel(subs),1); pw_hi = nan(numel(subs),1);
cn_lo = nan(numel(subs),1); cn_hi = nan(numel(subs),1);
for i = 1:numel(subs)
    rl = R(strcmp(R.subject,subs{i}) & strcmp(R.band,band) & strcmp(R.load,'low'), :);
    rh = R(strcmp(R.subject,subs{i}) & strcmp(R.band,band) & strcmp(R.load,'high'), :);
    if ~isempty(rl), pw_lo(i) = rl.fm_power(1); cn_lo(i) = rl.dwpli(1); end
    if ~isempty(rh), pw_hi(i) = rh.fm_power(1); cn_hi(i) = rh.dwpli(1); end
end
if cfg.log_power, pw_lo = log(pw_lo); pw_hi = log(pw_hi); end
ok = ~isnan(pw_lo) & ~isnan(pw_hi) & ~isnan(cn_lo) & ~isnan(cn_hi);
d_pw = pw_hi(ok) - pw_lo(ok);
d_cn = cn_hi(ok) - cn_lo(ok);
n = sum(ok);

[rho, p_rho] = corr(d_pw, d_cn, 'Type','Pearson');
fprintf('\nPart 1: per-subject beta power-change vs dwPLI-change correlation (N=%d)\n', n);
fprintf('  Pearson r = %.3f, p = %.4f\n', rho, p_rho);
fprintf('  (A significant POSITIVE r would support the amplitude-noise-driven-null account.)\n');

save(fullfile(cfg.deriv_dir,'amplitude_phase_correlation.mat'), 'd_pw','d_cn','rho','p_rho','n');

% ---------------- Part 2: amplitude-normalized dwPLI/PLV sensitivity check ----------------
addpath(cfg.eeglab_dir); eeglab nogui;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;
subjects = unique(manifest.subject, 'stable');
conds = {'low','high'};
rows = {};

for su = 1:numel(subjects)
    subj = subjects(su);
    keys = manifest.key(manifest.subject == subj);
    acc = struct('low',[],'high',[]);
    refEEG = [];

    for k = 1:numel(keys)
        key = char(keys(k));
        EEG = pop_loadset('filename',[key '_clean.set'],'filepath',cfg.deriv_dir);
        T = load(fullfile(cfg.deriv_dir,[key '_trials.mat'])); trials = T.trials;
        if isempty(refEEG), refEEG = EEG; end
        srate = EEG.srate; win = round(cfg.maint_window*srate); wlen = win(2)-win(1)+1;
        Xf = util_bandpass(double(EEG.data), srate, cfg.bands.beta, cfg.filt_order_factor);

        for c = 1:2
            sel = trials([trials.good] & strcmp({trials.load},conds{c}));
            for t = 1:numel(sel)
                a = sel(t).lock_sample + win(1); z = a + wlen - 1;
                if a < 1 || z > size(Xf,2), continue; end
                seg = Xf(:, a:z);
                % --- per-trial amplitude normalization (z-score), per channel ---
                seg = (seg - mean(seg,2)) ./ std(seg,0,2);
                Zseg = hilbert(seg.').';   % Hilbert on the normalized windowed segment
                acc.(conds{c}) = cat(3, acc.(conds{c}), Zseg);
            end
        end
    end
    if isempty(refEEG), continue; end

    fr = util_chanidx(refEEG, cfg.roi_frontal);
    te = util_chanidx(refEEG, cfg.roi_temporal);

    for c = 1:2
        seg = acc.(conds{c}); ntr = size(seg,3);
        if ntr < 5, continue; end
        wlen2 = size(seg,2);
        np_ = numel(fr)*numel(te);
        dwpli_pairs = zeros(1,np_); plv_pairs = zeros(1,np_); p = 0;
        for i = 1:numel(fr)
            zi = reshape(seg(fr(i),:,:), wlen2, ntr);
            for j = 1:numel(te)
                p = p+1;
                zj = reshape(seg(te(j),:,:), wlen2, ntr);
                X = zi .* conj(zj); im = imag(X);
                s_im = sum(im,2); s_abs = sum(abs(im),2); s_sq = sum(im.^2,2);
                dwpli_t = (s_im.^2 - s_sq) ./ (s_abs.^2 - s_sq);
                plv_t   = abs(mean(X ./ abs(X), 2));
                dwpli_pairs(p) = mean(dwpli_t,'omitnan');
                plv_pairs(p)   = mean(plv_t,'omitnan');
            end
        end
        rows(end+1,:) = {char(subj), conds{c}, mean(dwpli_pairs,'omitnan'), ...
            mean(plv_pairs,'omitnan'), ntr}; %#ok<SAGROW>
    end
    fprintf('amplitude-normalized done %s\n', subj);
end

Rn = cell2table(rows, 'VariableNames', {'subject','load','dwpli','plv','ntrials'});
writetable(Rn, fullfile(cfg.deriv_dir,'conn_results_ampnorm.csv'));

% ---- compare vs primary (unnormalized) beta result, same style as s04 ----
rng(1); subs2 = unique(Rn.subject);
lo = nan(numel(subs2),1); hi = nan(numel(subs2),1);
for i = 1:numel(subs2)
    rl = Rn(strcmp(Rn.subject,subs2{i}) & strcmp(Rn.load,'low'), :);
    rh = Rn(strcmp(Rn.subject,subs2{i}) & strcmp(Rn.load,'high'), :);
    if ~isempty(rl), lo(i) = rl.dwpli(1); end
    if ~isempty(rh), hi(i) = rh.dwpli(1); end
end
ok2 = ~isnan(lo) & ~isnan(hi);
d2 = hi(ok2) - lo(ok2); n2 = sum(ok2);
tobs2 = mean(d2)/(std(d2)/sqrt(n2));
np_ = cfg.stat.n_perm; tnull2 = zeros(np_,1);
for pp = 1:np_
    dp = d2 .* sign(rand(n2,1)-0.5);
    tnull2(pp) = mean(dp)/(std(dp)/sqrt(n2));
end
pval2 = mean(abs(tnull2) >= abs(tobs2));
dz2 = mean(d2)/std(d2);

fprintf('\nPart 2: amplitude-normalized beta dwPLI, High-Low (N=%d)\n', n2);
fprintf('  dz = %.3f, t = %.3f, p = %.4f  (compare to primary dz=0.02, p=0.94)\n', dz2, tobs2, pval2);

save(fullfile(cfg.deriv_dir,'ampnorm_stats.mat'), 'Rn','dz2','tobs2','pval2','n2');
