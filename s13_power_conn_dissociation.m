% s13_power_conn_dissociation.m -- Reviewer-requested check (C4): formal
% statistical test of the power-vs-connectivity "dissociation" claim, rather
% than inferring it from differential significance across two separate tests
% (comparing a significant p to a non-significant p is not itself a test of
% a difference; Nieuwenhuis et al., 2011).
%
% For the beta band and the same N=22 subjects, computes:
%   z_power_i = (fm_power_diff_i) / SD(fm_power_diff)   -- per-subject dz units
%   z_conn_i  = (dwpli_diff_i)    / SD(dwpli_diff)
% observed statistic = mean(z_power) - mean(z_conn) = dz_power - dz_conn
%
% Null distribution: the SAME random sign flip is applied to both z_power_i
% and z_conn_i for a given subject on each permutation (preserves whatever
% within-subject power-connectivity relationship exists under H0 of no
% differential load-sensitivity), 5000 reps -- same n_perm as elsewhere.

clear; clc;
cfg = config();
load(fullfile(cfg.deriv_dir,'conn_results.mat'), 'R');
rng(1);

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

dz_pw = mean(d_pw)/std(d_pw);
dz_cn = mean(d_cn)/std(d_cn);
z_pw = d_pw / std(d_pw);   % per-subject values, in dz units
z_cn = d_cn / std(d_cn);

obs = mean(z_pw) - mean(z_cn);   % = dz_pw - dz_cn

np_ = cfg.stat.n_perm; null = zeros(np_,1);
for pp = 1:np_
    s = sign(rand(n,1)-0.5);         % ONE sign per subject, applied to BOTH measures jointly
    null(pp) = mean(z_pw.*s) - mean(z_cn.*s);
end
pval = mean(abs(null) >= abs(obs));

fprintf('\nPower-vs-connectivity dissociation test (beta band), N = %d\n', n);
fprintf('  dz (power)         = %.3f\n', dz_pw);
fprintf('  dz (connectivity)  = %.3f\n', dz_cn);
fprintf('  Observed contrast (dz_power - dz_conn) = %.3f\n', obs);
fprintf('  Permutation p (joint sign-flip, %d perms) = %.4f\n', np_, pval);

save(fullfile(cfg.deriv_dir,'power_conn_dissociation_test.mat'), ...
    'dz_pw','dz_cn','obs','pval','n');
