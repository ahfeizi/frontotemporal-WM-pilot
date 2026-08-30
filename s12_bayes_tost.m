% s12_bayes_tost.m -- Reviewer-requested checks (C1, C5): sub-band-resolved
% Bayes factors and TOST equivalence tests for frontotemporal connectivity,
% computed IDENTICALLY for the main bands (theta, beta) and the beta
% sub-bands (beta_low, beta_high) from s09, so the "robust to how beta was
% defined" claim in Section 3.4 can be checked at each sub-band individually
% -- not just on the pooled 13-30 Hz band as in the current Table 4.
%
% BF10: JZS one-sample Bayes factor (Rouder et al., 2009), Cauchy(0,r) prior
%       on effect size, r = 0.707 (same default used to obtain BF01 ~ 4.5).
% TOST: two one-sided tests (Lakens, 2017) against equivalence bound
%       dz = 0.4 (same bound already reported in Section 3.6 / Table 4).
%
% Requires: conn_results.mat (from s03) and conn_results_betasplit.mat (from s09)
% to already exist in cfg.deriv_dir.

clear; clc;
cfg = config();
rng(1);
Delta   = 0.4;    % equivalence bound, Cohen's dz
r_prior = 0.707;  % Cauchy prior scale

S1 = load(fullfile(cfg.deriv_dir,'conn_results.mat'));           R_main = S1.R;
S2 = load(fullfile(cfg.deriv_dir,'conn_results_betasplit.mat')); R_sub  = S2.R;

metrics = {'dwpli','plv'};
rows = {};

% ---- main bands: theta, beta (reproduces/extends the current Table 4) ----
for b = {'theta','beta'}
    bn = b{1};
    for m = 1:numel(metrics)
        mt = metrics{m};
        d = paired_diff(R_main, bn, mt);
        rows(end+1,:) = eval_band(bn, mt, d, Delta, r_prior, cfg.stat.n_perm);
    end
end

% ---- beta sub-bands (the analysis the reviewer says is missing) ----
for b = {'beta_low','beta_high'}
    bn = b{1};
    for m = 1:numel(metrics)
        mt = metrics{m};
        d = paired_diff(R_sub, bn, mt);
        rows(end+1,:) = eval_band(bn, mt, d, Delta, r_prior, cfg.stat.n_perm);
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'band','metric','n','dz','t','p_perm','BF10','BF01','p_TOST','equivalent_at_0p4'});
disp(T);
writetable(T, fullfile(cfg.deriv_dir,'bayes_tost_expanded.csv'));
fprintf('\nSaved bayes_tost_expanded.csv -- expanded Table 4 (main bands + beta sub-bands)\n');

% ======================= local functions =======================
function d = paired_diff(R, bn, mt)
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
end

function row = eval_band(bn, mt, d, Delta, r_prior, n_perm)
    n  = numel(d);
    dz = mean(d)/std(d);
    se = std(d)/sqrt(n);
    tobs = mean(d)/se;
    v = n - 1;

    % Sign-flip permutation p, identical convention to s04_stats.m / s09
    tnull = zeros(n_perm,1);
    for pp = 1:n_perm
        dp = d .* sign(rand(n,1)-0.5);
        tnull(pp) = mean(dp)/(std(dp)/sqrt(n));
    end
    p_perm = mean(abs(tnull) >= abs(tobs));

    BF10 = jzs_bf10(tobs, n, r_prior);
    BF01 = 1/BF10;

    bound_raw = Delta * std(d);             % equivalence bound on the raw-difference scale
    t_lo = (mean(d) - (-bound_raw)) / se;   % H0: true diff <= -bound
    t_hi = (mean(d) -   bound_raw ) / se;   % H0: true diff >=  bound
    p_lo = 1 - tcdf(t_lo, v);
    p_hi = tcdf(t_hi, v);
    p_tost = max(p_lo, p_hi);
    equiv = p_tost < 0.05;

    row = {bn, mt, n, dz, tobs, p_perm, BF10, BF01, p_tost, equiv};
end

function bf10 = jzs_bf10(t, n, r)
    % Rouder et al. (2009), one-sample JZS Bayes factor.
    % g ~ Scale-Inv-chi2(1, r^2), i.e. density f(g) = (r/sqrt(2*pi)) * g^(-3/2) * exp(-r^2/(2*g)).
    % num   = null-model marginal likelihood (unnormalized; proportional to the
    %         central-t density with v df evaluated at the observed t).
    % denom = alternative-model marginal likelihood (unnormalized; Cauchy(0,r)
    %         effect-size prior via the g-mixture), obtained by integrating over g.
    % BF01 = num/denom (null over alternative); BF10 = 1/BF01 (alternative over null).
    v = n - 1;
    num = (1 + t^2/v)^(-(v+1)/2);
    integrand = @(g) (1+n.*g).^(-0.5) .* ...
        (1 + t^2 ./ ((1+n.*g).*v)).^(-(v+1)/2) .* ...
        (r/sqrt(2*pi)) .* g.^(-1.5) .* exp(-(r^2)./(2.*g));
    denom = integral(integrand, 0, Inf);
    bf01 = num / denom;
    bf10 = 1 / bf01;
end
