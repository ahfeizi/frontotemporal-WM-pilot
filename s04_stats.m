% s04_stats.m  --  Group-level within-subject contrast: HIGH vs LOW load.
% Non-parametric sign-flip permutation test on the paired subject-level
% differences, for each band and each metric (dwpli, plv, fm_power).
% Reports mean per condition, paired diff, Cohen's dz, t, and permutation p.

clear; clc;
cfg = config();
load(fullfile(cfg.deriv_dir,'conn_results.mat'), 'R');
rng(1);

bands   = fieldnames(cfg.bands);
metrics = {'dwpli','plv','fm_power'};
subs    = unique(R.subject);
results = struct();

fprintf('\n band   metric     n   low      high     dz     t       p\n');
fprintf('-------------------------------------------------------------\n');
for b = 1:numel(bands)
  bn = bands{b};
  for m = 1:numel(metrics)
    mt = metrics{m};
    lo = nan(numel(subs),1); hi = nan(numel(subs),1);
    for i = 1:numel(subs)
      rl = R(strcmp(R.subject,subs{i}) & strcmp(R.band,bn) & strcmp(R.load,'low'), :);
      rh = R(strcmp(R.subject,subs{i}) & strcmp(R.band,bn) & strcmp(R.load,'high'), :);
      if ~isempty(rl), lo(i) = rl.(mt)(1); end
      if ~isempty(rh), hi(i) = rh.(mt)(1); end
    end
    if strcmp(mt,'fm_power') && cfg.log_power
        lo = log(lo); hi = log(hi);
    end
    ok = ~isnan(lo) & ~isnan(hi);
    d  = hi(ok) - lo(ok); n = sum(ok);
    if n < 3, fprintf('%6s %-9s  n=%d too few\n', bn, mt, n); continue; end

    tobs = mean(d) / (std(d)/sqrt(n));
    np = cfg.stat.n_perm; tnull = zeros(np,1);
    for pp = 1:np
      dp = d .* sign(rand(n,1)-0.5);
      tnull(pp) = mean(dp) / (std(dp)/sqrt(n));
    end
    switch cfg.stat.tail
      case 0,  pval = mean(abs(tnull) >= abs(tobs));
      case 1,  pval = mean(tnull >= tobs);
      case -1, pval = mean(tnull <= tobs);
    end

    results.(bn).(mt) = struct('n',n,'mean_low',mean(lo(ok)),'mean_high',mean(hi(ok)), ...
        'mean_diff',mean(d),'cohen_dz',mean(d)/std(d),'t',tobs,'p',pval);
    fprintf('%6s %-9s %3d  %7.4f  %7.4f  %5.2f  %6.2f  %.4f\n', ...
        bn, mt, n, mean(lo(ok)), mean(hi(ok)), mean(d)/std(d), tobs, pval);
  end
end

save(fullfile(cfg.deriv_dir,'stats_results.mat'), 'results');
fprintf('\nSaved stats_results.mat\n');
