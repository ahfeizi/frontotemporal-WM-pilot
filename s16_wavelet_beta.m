% s16_wavelet_beta.m -- Reviewer-requested reanalysis (B1): replace the
% Hilbert-transform-based beta connectivity estimates with a Morlet-wavelet
% decomposition at multiple centre frequencies spanning 13-30 Hz, which
% satisfies the narrowband condition at each frequency and is the
% field-standard approach for wideband phase analysis (Bruns, 2004),
% avoiding the Bedrosian's-theorem concern raised for broadband Hilbert
% phase estimation on a ~17 Hz-wide band (fractional bandwidth ~0.79).
%
% For each centre frequency f (13:2:29 Hz, i.e. 2 Hz spacing, 9 frequencies;
% adjust cf_freqs below for a different grid), a complex Morlet wavelet
% (5 cycles) is convolved with each subject's continuous, already-
% preprocessed signal (same _clean.set files used throughout the pilot).
% dwPLI/PLV are computed per frequency bin exactly as in s03_connectivity.m,
% then averaged across frequency bins within a band:
%   beta      -> bins 13-29 Hz  (replaces the pooled 13-30 Hz Hilbert result)
%   beta_low  -> bins 13-19 Hz  (replaces the low-beta Hilbert result)
%   beta_high -> bins 21-29 Hz  (replaces the high-beta Hilbert result)
% giving one connectivity value per subject/band/load, directly comparable
% to (and intended to replace, pending comparison) the Hilbert-based values
% currently in Tables 1, 2, and 4.
%
% Theta is NOT re-analyzed here: its narrower fractional bandwidth (~0.55)
% was already judged less exposed to the Bedrosian concern (Methods 2.4),
% and the reviewer's framework explicitly allows retaining Hilbert for theta.
%
% RUNTIME NOTE: this does 9 frequencies x time-domain wavelet convolution
% over the full continuous recording, per subject/run. Expect this to take
% noticeably longer than s03 (minutes, not the hours a full ICA rerun would
% cost), but budget accordingly for the full N=22 sample.

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

% ---------------- wavelet frequency grid ----------------
cf_freqs = 13:2:29;      % centre frequencies, 2 Hz spacing, 13-29 Hz
n_cycles = 5;            % fixed cycle count (5-7 suggested by the reviewer;
                          % 5 chosen for adequate time resolution given the
                          % 2 s maintenance window)

band_membership = struct( ...
    'beta',      cf_freqs >= 13 & cf_freqs <= 29, ...
    'beta_low',  cf_freqs >= 13 & cf_freqs <= 19, ...
    'beta_high', cf_freqs >= 21 & cf_freqs <= 29);
band_names = fieldnames(band_membership);

subjects = unique(manifest.subject, 'stable');
conds = {'low','high'};
rows = {};

for su = 1:numel(subjects)
    subj = subjects(su);
    keys = manifest.key(manifest.subject == subj);

    % accumulate, PER FREQUENCY BIN, complex analytic segments per load cond
    acc = cell(1, numel(cf_freqs));
    for f = 1:numel(cf_freqs), acc{f} = struct('low',[],'high',[]); end
    refEEG = [];

    for k = 1:numel(keys)
        key = char(keys(k));
        EEG = pop_loadset('filename',[key '_clean.set'],'filepath',cfg.deriv_dir);
        T = load(fullfile(cfg.deriv_dir,[key '_trials.mat'])); trials = T.trials;
        if isempty(refEEG), refEEG = EEG; end
        srate = EEG.srate; win = round(cfg.maint_window*srate); wlen = win(2)-win(1)+1;
        data = double(EEG.data);   % nchan x ntime

        for f = 1:numel(cf_freqs)
            fc = cf_freqs(f);
            sigma_t = n_cycles / (2*pi*fc);
            t = -3*sigma_t : 1/srate : 3*sigma_t;
            wavelet = exp(2i*pi*fc.*t) .* exp(-t.^2/(2*sigma_t^2));
            wavelet = wavelet / sqrt(sum(abs(wavelet).^2));   % unit energy

            Z = zeros(size(data));
            for ch = 1:size(data,1)
                Z(ch,:) = conv(data(ch,:), wavelet, 'same');
            end

            for c = 1:2
                sel = trials([trials.good] & strcmp({trials.load},conds{c}));
                for t_i = 1:numel(sel)
                    a = sel(t_i).lock_sample + win(1); z = a + wlen - 1;
                    if a < 1 || z > size(Z,2), continue; end
                    acc{f}.(conds{c}) = cat(3, acc{f}.(conds{c}), Z(:, a:z));
                end
            end
        end
    end
    if isempty(refEEG), continue; end

    fr = util_chanidx(refEEG, cfg.roi_frontal);
    te = util_chanidx(refEEG, cfg.roi_temporal);

    % ---- per-frequency-bin dwPLI/PLV, then average within each band ----
    np_pairs = numel(fr)*numel(te);
    for bn_i = 1:numel(band_names)
        bn = band_names{bn_i};
        fbins = find(band_membership.(bn));
        for c = 1:2
            dwpli_freqavg = zeros(1,np_pairs); plv_freqavg = zeros(1,np_pairs);
            nf_used = 0;
            for f = fbins
                seg = acc{f}.(conds{c}); ntr = size(seg,3);
                if ntr < 5, continue; end
                wlen2 = size(seg,2);
                dwpli_pairs = zeros(1,np_pairs); plv_pairs = zeros(1,np_pairs); p = 0;
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
                dwpli_freqavg = dwpli_freqavg + dwpli_pairs;
                plv_freqavg   = plv_freqavg   + plv_pairs;
                nf_used = nf_used + 1;
            end
            if nf_used == 0, continue; end
            dwpli_freqavg = dwpli_freqavg / nf_used;
            plv_freqavg   = plv_freqavg   / nf_used;
            rows(end+1,:) = {char(subj), conds{c}, bn, ...
                mean(dwpli_freqavg,'omitnan'), mean(plv_freqavg,'omitnan'), ...
                nf_used, {dwpli_freqavg}, {plv_freqavg}}; %#ok<SAGROW>
        end
    end
    fprintf('wavelet done %s\n', subj);
end

R = cell2table(rows, 'VariableNames', {'subject','load','band','dwpli','plv', ...
    'n_freq_bins','dwpli_pairs','plv_pairs'});
save(fullfile(cfg.deriv_dir,'conn_results_wavelet.mat'),'R','cfg','cf_freqs','n_cycles');
writetable(R(:,1:6), fullfile(cfg.deriv_dir,'conn_results_wavelet.csv'));

% ---- paired permutation stats, same style as s04/s09 ----
rng(1); metrics = {'dwpli','plv'};
fprintf('\n band        metric   n   low      high     dz     t       p\n');
fprintf('-------------------------------------------------------------------\n');
results = struct();
for bn_i = 1:numel(band_names)
    bn = band_names{bn_i};
    for m = 1:numel(metrics)
        mt = metrics{m};
        subs = unique(R.subject);
        lo = nan(numel(subs),1); hi = nan(numel(subs),1);
        for i = 1:numel(subs)
            rl = R(strcmp(R.subject,subs{i}) & strcmp(R.band,bn) & strcmp(R.load,'low'), :);
            rh = R(strcmp(R.subject,subs{i}) & strcmp(R.band,bn) & strcmp(R.load,'high'), :);
            if ~isempty(rl), lo(i) = rl.(mt)(1); end
            if ~isempty(rh), hi(i) = rh.(mt)(1); end
        end
        ok = ~isnan(lo) & ~isnan(hi);
        d = hi(ok) - lo(ok); n = sum(ok);
        if n < 3, continue; end
        tobs = mean(d)/(std(d)/sqrt(n));
        np_ = cfg.stat.n_perm; tnull = zeros(np_,1);
        for pp = 1:np_
            dp = d .* sign(rand(n,1)-0.5);
            tnull(pp) = mean(dp)/(std(dp)/sqrt(n));
        end
        pval = mean(abs(tnull) >= abs(tobs));
        dz = mean(d)/std(d);
        results.(bn).(mt) = struct('n',n,'dz',dz,'t',tobs,'p',pval);
        fprintf('%-11s %-8s %3d  %7.4f  %7.4f  %5.2f  %6.2f  %.4f\n', ...
            bn, mt, n, mean(lo(ok)), mean(hi(ok)), dz, tobs, pval);
    end
end
save(fullfile(cfg.deriv_dir,'stats_results_wavelet.mat'), 'results');

% ---- BF + TOST, reusing the corrected JZS formula verified in s12 ----
Delta = 0.4; r_prior = 0.707;
rows_bf = {};
for bn_i = 1:numel(band_names)
    bn = band_names{bn_i};
    for m = 1:numel(metrics)
        mt = metrics{m};
        if ~isfield(results,bn) || ~isfield(results.(bn),mt), continue; end
        r_ = results.(bn).(mt);
        n = r_.n; tobs = r_.t; v = n-1;

        num = (1 + tobs^2/v)^(-(v+1)/2);
        integrand = @(g) (1+n.*g).^(-0.5) .* ...
            (1 + tobs^2 ./ ((1+n.*g).*v)).^(-(v+1)/2) .* ...
            (r_prior/sqrt(2*pi)) .* g.^(-1.5) .* exp(-(r_prior^2)./(2.*g));
        denom = integral(integrand, 0, Inf);
        bf01 = num/denom; bf10 = 1/bf01;

        % recompute the paired differences directly for TOST (need raw SD)
        subs = unique(R.subject);
        lo = nan(numel(subs),1); hi = nan(numel(subs),1);
        for i = 1:numel(subs)
            rl = R(strcmp(R.subject,subs{i}) & strcmp(R.band,bn) & strcmp(R.load,'low'), :);
            rh = R(strcmp(R.subject,subs{i}) & strcmp(R.band,bn) & strcmp(R.load,'high'), :);
            if ~isempty(rl), lo(i) = rl.(mt)(1); end
            if ~isempty(rh), hi(i) = rh.(mt)(1); end
        end
        ok = ~isnan(lo) & ~isnan(hi);
        d = hi(ok) - lo(ok);
        se = std(d)/sqrt(n);
        bound_raw = Delta*std(d);
        t_lo = (mean(d) - (-bound_raw))/se; t_hi = (mean(d) - bound_raw)/se;
        p_lo = 1 - tcdf(t_lo, v); p_hi = tcdf(t_hi, v);
        p_tost = max(p_lo, p_hi); equiv = p_tost < 0.05;

        rows_bf(end+1,:) = {bn, mt, n, r_.dz, r_.t, bf10, bf01, p_tost, equiv}; %#ok<SAGROW>
    end
end
Tbf = cell2table(rows_bf, 'VariableNames', ...
    {'band','metric','n','dz','t','BF10','BF01','p_TOST','equivalent_at_0p4'});
disp(Tbf);
writetable(Tbf, fullfile(cfg.deriv_dir,'bayes_tost_wavelet.csv'));
fprintf('\nSaved conn_results_wavelet.csv, stats_results_wavelet.mat, bayes_tost_wavelet.csv\n');
fprintf('Compare directly against the Hilbert-based values in Table 1/2/4 and bayes_tost_expanded.csv.\n');
