% s08_site_robustness.m -- Robustness check: does the beta-power load
% effect replicate at Fz and Cz, not just FCz?
% Re-uses the exact same accumulation/extraction logic as s03_connectivity.m
% and the exact same sign-flip permutation test as s04_stats.m.

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

sites = {'FCz','Fz','Cz'};
band  = cfg.bands.beta;   % [13 30]
subjects = unique(manifest.subject, 'stable');
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

        Xf = util_bandpass(double(EEG.data), srate, band, cfg.filt_order_factor);
        Z  = hilbert(Xf.').';

        for c = {'low','high'}
            cond = c{1};
            sel = trials([trials.good] & strcmp({trials.load},cond));
            for t = 1:numel(sel)
                a = sel(t).lock_sample + win(1); z = a + wlen - 1;
                if a < 1 || z > size(Z,2), continue; end
                acc.(cond) = cat(3, acc.(cond), Z(:, a:z));
            end
        end
    end
    if isempty(refEEG), continue; end

    for s = 1:numel(sites)
        ch = util_chanidx(refEEG, sites(s));
        if isempty(ch), continue; end
        for c = {'low','high'}
            cond = c{1};
            seg = acc.(cond); ntr = size(seg,3);
            if ntr < 5, continue; end
            wlen2 = size(seg,2);
            pw = mean(mean(abs(reshape(seg(ch,:,:),wlen2,ntr)).^2));
            rows(end+1,:) = {char(subj), sites{s}, cond, pw};
        end
    end
end

R = cell2table(rows, 'VariableNames', {'subject','site','load','beta_power'});

% ---- sign-flip permutation test per site (identical logic to s04_stats.m) ----
rng(1); np = cfg.stat.n_perm;
results = table();
for s = 1:numel(sites)
    site = sites{s};
    subs = unique(R.subject);
    lo = nan(numel(subs),1); hi = nan(numel(subs),1);
    for i = 1:numel(subs)
        rl = R(strcmp(R.subject,subs{i}) & strcmp(R.site,site) & strcmp(R.load,'low'), :);
        rh = R(strcmp(R.subject,subs{i}) & strcmp(R.site,site) & strcmp(R.load,'high'), :);
        if ~isempty(rl), lo(i) = rl.beta_power(1); end
        if ~isempty(rh), hi(i) = rh.beta_power(1); end
    end
    if cfg.log_power, lo = log(lo); hi = log(hi); end
    ok = ~isnan(lo) & ~isnan(hi);
    d = hi(ok) - lo(ok); n = sum(ok);

    tobs = mean(d)/(std(d)/sqrt(n));
    tnull = zeros(np,1);
    for pp = 1:np
        dp = d .* sign(rand(n,1)-0.5);
        tnull(pp) = mean(dp)/(std(dp)/sqrt(n));
    end
    pval = mean(abs(tnull) >= abs(tobs));
    dz = mean(d)/std(d);

    results = [results; {site, n, dz, tobs, pval}];
end
results.Properties.VariableNames = {'Site','N','dz','t','p'};
disp(results)
writetable(results, fullfile(cfg.deriv_dir,'fcz_robustness_check.csv'));
fprintf('\nSaved fcz_robustness_check.csv\n');