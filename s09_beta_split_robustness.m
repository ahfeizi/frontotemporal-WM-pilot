% s09_beta_split_robustness.m -- Reviewer-requested check: split the 13-30 Hz
% beta band into low-beta (13-20 Hz) and high-beta (21-30 Hz) and re-run the
% FULL analysis pipeline (power, dwPLI/PLV, paired permutation stats, and
% FWER-controlled pair-level test) separately for each sub-band.
% This does NOT modify the original theta/beta results already reported;
% it is an additional, self-contained robustness analysis.

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

subbands.beta_low  = [13 20];
subbands.beta_high = [21 30];
subnames = fieldnames(subbands);

subjects = unique(manifest.subject, 'stable');
conds = {'low','high'};
rows = {};

% ---------- 1) power + dwPLI/PLV extraction (same logic as s03) ----------
for su = 1:numel(subjects)
    subj = subjects(su);
    keys = manifest.key(manifest.subject == subj);
    acc = struct();
    for b = 1:numel(subnames), for c = 1:2, acc.(subnames{b}).(conds{c}) = []; end, end %#ok<ALIGN>
    refEEG = [];

    for k = 1:numel(keys)
        key = char(keys(k));
        EEG = pop_loadset('filename',[key '_clean.set'],'filepath',cfg.deriv_dir);
        T = load(fullfile(cfg.deriv_dir,[key '_trials.mat'])); trials = T.trials;
        if isempty(refEEG), refEEG = EEG; end
        srate = EEG.srate; win = round(cfg.maint_window*srate); wlen = win(2)-win(1)+1;

        for b = 1:numel(subnames)
            Xf = util_bandpass(double(EEG.data), srate, subbands.(subnames{b}), cfg.filt_order_factor);
            Z  = hilbert(Xf.').';
            for c = 1:2
                sel = trials([trials.good] & strcmp({trials.load},conds{c}));
                for t = 1:numel(sel)
                    a = sel(t).lock_sample + win(1); z = a + wlen - 1;
                    if a < 1 || z > size(Z,2), continue; end
                    acc.(subnames{b}).(conds{c}) = cat(3, acc.(subnames{b}).(conds{c}), Z(:, a:z));
                end
            end
        end
    end
    if isempty(refEEG), continue; end

    fr = util_chanidx(refEEG, cfg.roi_frontal);
    te = util_chanidx(refEEG, cfg.roi_temporal);
    fmIdx = util_chanidx(refEEG, {'FCz','Fz'}); fmIdx = fmIdx(1);

    for b = 1:numel(subnames)
        for c = 1:2
            seg = acc.(subnames{b}).(conds{c}); ntr = size(seg,3);
            if ntr < 5, continue; end
            wlen2 = size(seg,2);
            np = numel(fr)*numel(te);
            dwpli_pairs = zeros(1,np); plv_pairs = zeros(1,np); p = 0;
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
            fm_power = mean(mean(abs(reshape(seg(fmIdx,:,:),wlen2,ntr)).^2));
            rows(end+1,:) = {char(subj), conds{c}, subnames{b}, ...
                mean(dwpli_pairs,'omitnan'), mean(plv_pairs,'omitnan'), ...
                fm_power, ntr, {dwpli_pairs}, {plv_pairs}}; %#ok<SAGROW>
        end
    end
    fprintf('done %s\n', subj);
end

R = cell2table(rows, 'VariableNames', {'subject','load','band','dwpli','plv', ...
    'fm_power','ntrials','dwpli_pairs','plv_pairs'});
save(fullfile(cfg.deriv_dir,'conn_results_betasplit.mat'),'R','cfg');
writetable(R(:,1:7), fullfile(cfg.deriv_dir,'conn_results_betasplit.csv'));

% ---------- 2) paired permutation stats (same logic as s04) ----------
rng(1); metrics = {'dwpli','plv','fm_power'};
results = struct();
fprintf('\n band       metric     n   low      high     dz     t       p\n');
fprintf('-------------------------------------------------------------------\n');
for b = 1:numel(subnames)
  bn = subnames{b};
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
    if strcmp(mt,'fm_power') && cfg.log_power, lo = log(lo); hi = log(hi); end
    ok = ~isnan(lo) & ~isnan(hi);
    d = hi(ok)-lo(ok); n = sum(ok);
    if n < 3, continue; end
    tobs = mean(d)/(std(d)/sqrt(n));
    np_ = cfg.stat.n_perm; tnull = zeros(np_,1);
    for pp = 1:np_
        dp = d .* sign(rand(n,1)-0.5);
        tnull(pp) = mean(dp)/(std(dp)/sqrt(n));
    end
    pval = mean(abs(tnull) >= abs(tobs));
    results.(bn).(mt) = struct('n',n,'dz',mean(d)/std(d),'t',tobs,'p',pval);
    fprintf('%-10s %-9s %3d  %7.4f  %7.4f  %5.2f  %6.2f  %.4f\n', ...
        bn, mt, n, mean(lo(ok)), mean(hi(ok)), mean(d)/std(d), tobs, pval);
  end
end
save(fullfile(cfg.deriv_dir,'stats_results_betasplit.mat'),'results');

% ---------- 3) FWER-controlled pair-level test (same logic as s06) ----------
nF = numel(cfg.roi_frontal); nT = numel(cfg.roi_temporal); npair = nF*nT;
labs = {}; D = [];
for b = 1:numel(subnames)
    bn = subnames{b};
    Db = nan(numel(subjects), npair);
    for si = 1:numel(subjects)
        rl = R(strcmp(R.subject,subjects{si}) & strcmp(R.band,bn) & strcmp(R.load,'low'), :);
        rh = R(strcmp(R.subject,subjects{si}) & strcmp(R.band,bn) & strcmp(R.load,'high'), :);
        if isempty(rl) || isempty(rh), continue; end
        lo = rl.dwpli_pairs{1}; if iscell(lo), lo = lo{1}; end
        hi = rh.dwpli_pairs{1}; if iscell(hi), hi = hi{1}; end
        Db(si,:) = hi - lo;
    end
    D = [D Db]; %#ok<AGROW>
    for i = 1:nF, for j = 1:nT
        labs{end+1} = sprintf('%-9s %s-%s', bn, cfg.roi_frontal{i}, cfg.roi_temporal{j}); %#ok<AGROW,SAGROW>
    end, end
end
bad = any(isnan(D),2); D = D(~bad,:); n = size(D,1); Mcol = size(D,2);
tobs = mean(D)./(std(D)./sqrt(n));
np_ = cfg.stat.n_perm; maxnull = zeros(np_,1);
for p = 1:np_
    sgn = sign(rand(n,1)-0.5); Dp = D.*sgn;
    t = mean(Dp)./(std(Dp)./sqrt(n));
    maxnull(p) = max(abs(t));
end
thr = prctile(maxnull, 95);
p_fwer = arrayfun(@(x) mean(maxnull >= abs(x)), tobs);
sig = p_fwer < 0.05;

T2 = table(labs(:), tobs(:), p_fwer(:), sig(:), 'VariableNames', {'pair_band','t','p_fwer','sig_fwer'});
T2 = sortrows(T2,'t','descend','ComparisonMethod','abs');
fprintf('\nBeta-split FWER threshold |t| = %.2f | significant: %d of %d\n', thr, sum(sig), Mcol);
writetable(T2, fullfile(cfg.deriv_dir,'conn_cluster_results_betasplit.csv'));